"""MRML editing session, live measurements and non-overwriting revisions. BSD-3-Clause."""
import datetime
import json
import logging
import os
from pathlib import Path
import uuid
import numpy as np
import qt
import vtk
import slicer
from vtk.util.numpy_support import vtk_to_numpy
from . import geometry, stateio

MODULE = "DentalArchEditor"

class EditorSession:
    def __init__(self, changed=None):
        self.changed = changed
        self.parameter_node = None
        self.md_nodes = {}
        self.handle_node = self.curve_node = None
        self.state = self.baseline = None
        self.result = {"valid": False, "warnings": ["Import an annotated scene or load the synthetic demo."]}
        self._observers = []
        self._scene_observers = []
        self._busy = False
        self._dragging = False
        self._structural_invalid = False
        self._locators = {}
        self._history = []
        self._cursor = 0
        self._timer = qt.QTimer()
        self._timer.setSingleShot(True)
        self._timer.connect("timeout()", self._deferred_update)
        self._scene_observers.append(slicer.mrmlScene.AddObserver(slicer.mrmlScene.EndImportEvent, self._scene_loaded))
        self._scene_observers.append(slicer.mrmlScene.AddObserver(slicer.mrmlScene.StartCloseEvent, self._scene_closing))
        self._scene_observers.append(slicer.mrmlScene.AddObserver(slicer.mrmlScene.NodeRemovedEvent, self._node_removed))
        self.restore()

    def cleanup(self):
        self._timer.stop()
        self._detach()
        for tag in self._scene_observers:
            slicer.mrmlScene.RemoveObserver(tag)
        self._scene_observers.clear()
        self.changed = None

    def _detach(self):
        for node, tag in self._observers:
            node.RemoveObserver(tag)
        self._observers.clear()
        self._locators.clear()

    def _scene_closing(self, *_):
        self._timer.stop()
        self._detach()
        self.parameter_node = None
        self.md_nodes = {}
        self.handle_node = self.curve_node = None
        self.state = self.baseline = None
        self._history = []
        self._invalid("Scene closed. Open a revision or import a baseline.")

    def _scene_loaded(self, *_):
        self.restore()

    def _node_removed(self, *_):
        if self.parameter_node and not self._busy:
            self._timer.start(40)

    def _emit(self):
        if self.changed:
            self.changed(self.result)

    def _invalid(self, message):
        self.result = {"valid": False, "R": None, "A": None, "C": None, "D": None, "widths": {}, "warnings": [str(message)]}
        self._emit()
        return self.result

    def _unique_named_node(self, names):
        matches = [node for node in slicer.util.getNodesByClass("vtkMRMLNode") if node.GetName() in names]
        if len(matches) > 1:
            raise ValueError("Ambiguous scene mapping: " + ", ".join(names))
        return matches[0] if matches else None

    def _reference(self, role, node):
        self.parameter_node.SetNodeReferenceID(role, node.GetID() if node else None)

    def _new_markup(self, name, node_type, points, color, role, labels=None, size=.4):
        node = slicer.mrmlScene.AddNewNodeByClass(node_type, name)
        node.CreateDefaultDisplayNodes()
        if node.IsA("vtkMRMLMarkupsCurveNode"):
            node.SetCurveTypeToLinear()
            if hasattr(node, "SetNumberOfPointsPerInterpolatingSegment"):
                node.SetNumberOfPointsPerInterpolatingSegment(1)
        for i, point in enumerate(points):
            index = node.AddControlPoint(vtk.vtkVector3d(*[float(x) for x in point]))
            node.SetNthControlPointLabel(index, labels[i] if labels else str(i))
        if hasattr(node, "SetFixedNumberOfControlPoints"):
            node.SetFixedNumberOfControlPoints(True)
        display = node.GetDisplayNode()
        display.SetSelectedColor(*color)
        display.SetColor(*color)
        display.SetPointLabelsVisibility(False)
        display.SetPropertiesLabelVisibility(False)
        display.SetUseGlyphScale(False)
        display.SetGlyphSize(size)
        if hasattr(display, "SetCurveLineSizeMode"):
            display.SetCurveLineSizeMode(slicer.vtkMRMLMarkupsDisplayNode.UseLineDiameter)
            display.SetLineDiameter(.20 if role == "DisplayArch" else .10)
        if hasattr(display, "SetOccludedVisibility"):
            display.SetOccludedVisibility(True)
            display.SetOccludedOpacity(.55)
        for i in range(node.GetNumberOfMeasurements()):
            node.GetNthMeasurement(i).SetEnabled(False)
        node.SetAttribute("DentalArchEditor.Role", role)
        node.SetAttribute("Units", "native coordinate units; no millimetre calibration")
        self._reference(role, node)
        return node

    def import_baseline(self, md_json_path, arch_json_path):
        state, provenance = stateio.import_baseline(md_json_path, arch_json_path)
        result = geometry.evaluate(state)
        if not result["valid"]:
            raise ValueError("; ".join(result.get("errors", []) + result.get("warnings", [])))
        if abs(result["A"] - provenance["expected_A"]) > 1e-6:
            raise ValueError("Imported spline length does not reproduce the recorded baseline.")
        protected = [str(Path(md_json_path).resolve()), str(Path(arch_json_path).resolve())]
        scene_url = slicer.mrmlScene.GetURL()
        if scene_url:
            protected.append(str(Path(scene_url).resolve()))
        self.initialize(state, provenance, protected)
        return self.result

    def initialize(self, state, provenance=None, protected_paths=None):
        if self.parameter_node:
            raise ValueError("An editing session already exists. Open a fresh scene before importing another baseline.")
        result = geometry.evaluate(state)
        if not result["valid"]:
            raise ValueError("; ".join(result.get("errors", []) + result.get("warnings", [])))
        # Validate all existing mapping candidates before any mutation.
        mapping = {}
        for tooth, points in state["teeth"].items():
            baseline = self._unique_named_node([tooth + "_MD", tooth + "_MD_CANDIDATE"])
            crown = self._unique_named_node([tooth + "_crown", tooth + "_isolated_crown"])
            if baseline:
                if baseline.GetParentTransformNode() or baseline.GetNumberOfControlPoints() != 2:
                    raise ValueError("Baseline landmark node has a transform or incorrect point count: " + tooth)
                xyz = np.array([baseline.GetNthControlPointPositionVector(i) for i in range(2)])
                if not np.allclose(xyz, [points["M"], points["D"]], atol=1e-6, rtol=0):
                    raise ValueError("Scene and metadata landmark coordinates disagree for " + tooth)
            if crown and crown.GetParentTransformNode():
                raise ValueError("Use untransformed crown coordinates for " + tooth)
            mapping[tooth] = (baseline, crown)
        self._busy = True
        try:
            parameter = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLScriptedModuleNode", "DentalArchEditorSession")
            parameter.SetModuleName(MODULE)
            self.parameter_node = parameter
            self.baseline = stateio.clone(state)
            self.state = stateio.clone(state)
            parameter.SetParameter("BaselineJSON", stateio.dumps(self.baseline))
            parameter.SetParameter("ProvenanceJSON", stateio.dumps(provenance or {}))
            parameter.SetParameter("ProtectedPathsJSON", stateio.dumps(protected_paths or []))
            parameter.SetParameter("SnapEnabled", "1")
            parameter.SetParameter("EditMode", "points")
            parameter.SetParameter("SelectedTooth", next(iter(state["teeth"])))
            for node in slicer.util.getNodesByClass("vtkMRMLMarkupsNode"):
                node.SetLocked(True)
                node.GetDisplayNode().SetVisibility(False)
            self.md_nodes = {}
            for tooth, points in state["teeth"].items():
                original, crown = mapping[tooth]
                if original is None:
                    original = self._new_markup("Baseline_" + tooth, "vtkMRMLMarkupsLineNode", [points["M"], points["D"]], (.55, .55, .55), "BaselineMD_" + tooth)
                self._reference("BaselineMD_" + tooth, original)
                original.SetLocked(True)
                original.GetDisplayNode().SetColor(.55, .55, .55)
                original.GetDisplayNode().SetSelectedColor(.55, .55, .55)
                original.GetDisplayNode().SetOpacity(.55)
                original.GetDisplayNode().SetPointLabelsVisibility(False)
                original.GetDisplayNode().SetVisibility(False)
                self._reference("Crown_" + tooth, crown)
                node = self._new_markup("Edit_" + tooth, "vtkMRMLMarkupsLineNode", [points["M"], points["D"]], (.22, .96, .65), "MD_" + tooth, [tooth + "_M", tooth + "_D"])
                node.GetDisplayNode().SetSelectedColor(1, .61, .25)
                node.SetNthControlPointSelected(0, False)
                node.SetNthControlPointSelected(1, True)
                self.md_nodes[tooth] = node
            base_arch = self._new_markup("Baseline_reference_arch", "vtkMRMLMarkupsCurveNode", result["curve_native"], (.55, .55, .55), "BaselineArch", size=.02)
            base_arch.SetLocked(True)
            base_arch.GetDisplayNode().SetVisibility(False)
            self.curve_node = self._new_markup("Edited_reference_arch", "vtkMRMLMarkupsCurveNode", result["curve_native"], (.1, .85, 1), "DisplayArch", size=.015)
            self.curve_node.SetLocked(True)
            self.handle_node = self._new_markup("Arch_shape_handles", "vtkMRMLMarkupsFiducialNode", result["controls_native"], (1, .85, .25), "ArchHandles", size=.65)
            self._history = [stateio.clone(self.state)]
            self._cursor = 0
            # Avoid retaining the baseline MRB as the default scene-save target.
            slicer.mrmlScene.SetURL("")
        finally:
            self._busy = False
        self._install_observers()
        self.set_edit_mode("points")
        self.evaluate_now(commit=False)
        self._persist()

    def restore(self):
        self._timer.stop()
        self._detach()
        parameters = [n for n in slicer.util.getNodesByClass("vtkMRMLScriptedModuleNode") if n.GetModuleName() == MODULE]
        if not parameters:
            return False
        if len(parameters) != 1:
            self._invalid("Multiple editor sessions in one scene are unsupported.")
            return False
        self.parameter_node = parameters[0]
        try:
            self.baseline = json.loads(self.parameter_node.GetParameter("BaselineJSON"))
            self.state = json.loads(self.parameter_node.GetParameter("CurrentJSON"))
            self.md_nodes = {t: self.parameter_node.GetNodeReference("MD_" + t) for t in self.state["teeth"]}
            self.handle_node = self.parameter_node.GetNodeReference("ArchHandles")
            self.curve_node = self.parameter_node.GetNodeReference("DisplayArch")
            saved = self.parameter_node.GetParameter("HistoryJSON")
            history = json.loads(saved) if saved else {"states": [self.state], "cursor": 0}
            self._history, self._cursor = history["states"], history["cursor"]
            if not all(self.md_nodes.values()) or not self.handle_node or not self.curve_node:
                raise ValueError("Revision is missing required markup nodes.")
            self._install_observers()
            self.set_edit_mode(self.parameter_node.GetParameter("EditMode") or "points")
            self.evaluate_now(commit=False)
            return True
        except Exception as error:
            self._invalid(error)
            return False

    def _install_observers(self):
        self._detach()
        for node in list(self.md_nodes.values()) + [self.handle_node]:
            if not node:
                continue
            for event_name in ["PointModifiedEvent", "PointAddedEvent", "PointRemovedEvent", "PointPositionDefinedEvent", "PointPositionUndefinedEvent", "TransformModifiedEvent"]:
                event = getattr(node, event_name, None)
                if event is not None:
                    self._observers.append((node, node.AddObserver(event, self._point_changed)))
            for event_name, callback in [("PointStartInteractionEvent", self._start_drag), ("PointEndInteractionEvent", self._end_drag)]:
                if hasattr(node, event_name):
                    self._observers.append((node, node.AddObserver(getattr(node, event_name), callback)))
        for tooth in self.md_nodes:
            crown = self.parameter_node.GetNodeReference("Crown_" + tooth)
            if crown:
                self._observers.append((crown, crown.AddObserver(vtk.vtkCommand.ModifiedEvent, self._point_changed)))
                self._observers.append((crown, crown.AddObserver(crown.TransformModifiedEvent, self._point_changed)))

    def _start_drag(self, *_):
        if not self._busy:
            self._dragging = True

    def _end_drag(self, *_):
        if not self._busy:
            self._dragging = False
            self._timer.stop()
            self.evaluate_now(commit=True)

    def _point_changed(self, node, *_):
        if self._busy:
            return
        try:
            if self.parameter_node.GetParameter("SnapEnabled") == "1" and node in self.md_nodes.values():
                tooth = next(t for t, n in self.md_nodes.items() if n == node)
                self._snap_to_crown(tooth)
        except Exception as error:
            self._structural_invalid = True
            self._invalid(error)
            return
        # A short coalescing timer keeps the event loop responsive during dragging.
        if not self._timer.isActive():
            self._timer.start(40)

    def _deferred_update(self):
        self.evaluate_now(commit=not self._dragging)

    def _check_node(self, node, count):
        if not node or node.GetScene() != slicer.mrmlScene:
            raise ValueError("A required point node was removed.")
        if node.GetParentTransformNode():
            raise ValueError("Parent transforms on editable nodes are unsupported; use the imported native coordinate frame.")
        if node.GetNumberOfControlPoints() != count:
            raise ValueError("Point count changed. Undo or restore the selected tooth.")
        for i in range(count):
            if node.GetNthControlPointPositionStatus(i) != slicer.vtkMRMLMarkupsNode.PositionDefined:
                raise ValueError("A control point is undefined. Restore it before measuring.")

    def _locator(self, tooth):
        model = self.parameter_node.GetNodeReference("Crown_" + tooth)
        if not model:
            return None
        if model.GetParentTransformNode():
            raise ValueError("Crown " + tooth + " has a parent transform; native-frame editing requires untransformed crowns.")
        poly = model.GetPolyData()
        if not poly or poly.GetNumberOfCells() == 0:
            raise ValueError("The assigned crown model is empty: " + tooth)
        key = (model.GetID(), poly.GetMTime())
        if tooth not in self._locators or self._locators[tooth][0] != key:
            locator = vtk.vtkStaticCellLocator()
            locator.SetDataSet(poly)
            locator.BuildLocator()
            self._locators[tooth] = (key, locator)
        return self._locators[tooth][1]

    def _snap_to_crown(self, tooth):
        node = self.md_nodes[tooth]
        self._check_node(node, 2)
        locator = self._locator(tooth)
        if locator is None:
            return
        self._busy = True
        try:
            for i in range(2):
                p = node.GetNthControlPointPositionVector(i)
                closest = [0., 0., 0.]
                cell, sub, distance = vtk.reference(0), vtk.reference(0), vtk.reference(0.)
                locator.FindClosestPoint(p, closest, cell, sub, distance)
                node.SetNthControlPointPosition(i, *closest)
        finally:
            self._busy = False

    def _capture(self):
        state = stateio.clone(self.state)
        for tooth, node in self.md_nodes.items():
            self._check_node(node, 2)
            self._locator(tooth)  # Also rejects unsupported transforms when snap is disabled.
            state["teeth"][tooth]["M"] = list(node.GetNthControlPointPositionVector(0))
            state["teeth"][tooth]["D"] = list(node.GetNthControlPointPositionVector(1))
        self._check_node(self.handle_node, len(state["arch"]["control_points_xy"]))
        center, axes = np.array(state["frame"]["center"]), np.array(state["frame"]["axes"])
        native = np.array([self.handle_node.GetNthControlPointPositionVector(i) for i in range(self.handle_node.GetNumberOfControlPoints())])
        old_xy = np.array(state["arch"]["control_points_xy"])
        old_native = np.column_stack([old_xy, np.full(len(old_xy), state["frame"]["plane_z"])]) @ axes.T + center
        changed = np.linalg.norm(native - old_native, axis=1) > 1e-10
        projected = ((native - center) @ axes)[:, :2]
        old_xy[changed] = projected[changed]
        state["arch"]["control_points_xy"] = old_xy.tolist()
        geometry.synchronize_endpoints(state)
        return state

    def _extra_warnings(self, state):
        warnings = []
        center, axes = np.array(state["frame"]["center"]), np.array(state["frame"]["axes"])
        for tooth, points in state["teeth"].items():
            model = self.parameter_node.GetNodeReference("Crown_" + tooth)
            if not model:
                warnings.append(tooth + ": no crown assigned; surface snapping unavailable.")
                continue
            if points.get("kind") in ("premolar", "molar"):
                z = ((vtk_to_numpy(model.GetPolyData().GetPoints().GetData()) - center) @ axes)[:, 2]
                h = (((np.array([points["M"], points["D"]]) - center) @ axes)[:, 2] - z.min()) / max(float(np.ptp(z)), 1e-12)
                if min(h) < 2 / 3 - 1e-5:
                    warnings.append(tooth + ": manual point is below the geometric upper-third band.")
            if points.get("kind") == "premolar":
                old = np.array(self.baseline["teeth"][tooth]["M"]) - self.baseline["teeth"][tooth]["D"]
                new = np.array(points["M"]) - points["D"]
                den = np.linalg.norm(old) * np.linalg.norm(new)
                if den and np.degrees(np.arccos(np.clip(old @ new / den, -1, 1))) > 5:
                    warnings.append(tooth + ": manual width direction differs from baseline by more than 5 degrees.")
        return warnings

    def evaluate_now(self, commit=True):
        if not self.parameter_node or self._busy:
            return self.result
        self._timer.stop()
        try:
            state = self._capture()
            self._structural_invalid = False
            result = geometry.evaluate(state)
            self.state = state
            result["warnings"] = list(result.get("warnings", [])) + self._extra_warnings(state)
            self.result = result
            self._busy = True
            try:
                if result.get("controls_native") is not None:
                    for i, p in enumerate(result["controls_native"]):
                        self.handle_node.SetNthControlPointPosition(i, *[float(x) for x in p])
                if result.get("curve_native") is not None and len(result["curve_native"]):
                    slicer.util.updateMarkupsControlPointsFromArray(self.curve_node, np.asarray(result["curve_native"]))
                self.curve_node.GetDisplayNode().SetSelectedColor(*((.1, .85, 1) if result["valid"] else (1, .25, .25)))
            finally:
                self._busy = False
            if commit:
                self._commit()
            self._persist()
            self._emit()
            return result
        except Exception as error:
            self._structural_invalid = True
            return self._invalid(error)

    def _commit(self):
        if self._history:
            prior = self._history[self._cursor]
            unchanged = all(np.allclose(prior["teeth"][t][r], self.state["teeth"][t][r], rtol=0, atol=1e-10) for t in self.state["teeth"] for r in ("M", "D"))
            unchanged &= np.allclose(prior["arch"]["control_points_xy"], self.state["arch"]["control_points_xy"], rtol=0, atol=1e-10)
            if unchanged:
                return
        self._history = self._history[:self._cursor + 1] + [stateio.clone(self.state)]
        self._history = self._history[-50:]
        self._cursor = len(self._history) - 1

    def _persist(self):
        if self.parameter_node:
            self.parameter_node.SetParameter("CurrentJSON", stateio.dumps(self.state))
            self.parameter_node.SetParameter("HistoryJSON", stateio.dumps({"states": self._history, "cursor": self._cursor}))
            self.parameter_node.SetParameter("EditorVersion", stateio.VERSION)

    def _apply(self, state):
        self._timer.stop()
        self._busy = True
        try:
            self.state = stateio.clone(state)
            for tooth, node in self.md_nodes.items():
                points = np.array([state["teeth"][tooth]["M"], state["teeth"][tooth]["D"]])
                slicer.util.updateMarkupsControlPointsFromArray(node, points)
            xy = np.array(state["arch"]["control_points_xy"])
            q = np.column_stack([xy, np.full(len(xy), state["frame"]["plane_z"])])
            points = q @ np.array(state["frame"]["axes"]).T + state["frame"]["center"]
            slicer.util.updateMarkupsControlPointsFromArray(self.handle_node, points)
        finally:
            self._busy = False
        self.evaluate_now(commit=False)
        self.set_edit_mode(self.parameter_node.GetParameter("EditMode"))

    def set_point(self, tooth, role, xyz, snap=True):
        if role not in ("M", "D"):
            raise ValueError("Point role must be M or D.")
        self._busy = True
        try:
            self.md_nodes[tooth].SetNthControlPointPosition(0 if role == "M" else 1, *[float(x) for x in xyz])
        finally:
            self._busy = False
        if snap:
            self._snap_to_crown(tooth)
        return self.evaluate_now(commit=True)

    def set_arch_control(self, index, xy):
        if not 0 < index < self.handle_node.GetNumberOfControlPoints() - 1:
            raise ValueError("Boundary handles follow the M1 mesial points and cannot be edited independently.")
        q = np.array([*xy, self.state["frame"]["plane_z"]], float)
        p = q @ np.array(self.state["frame"]["axes"]).T + self.state["frame"]["center"]
        self._busy = True
        try:
            self.handle_node.SetNthControlPointPosition(index, *p)
        finally:
            self._busy = False
        return self.evaluate_now(commit=True)

    def undo(self):
        if self._structural_invalid:
            self._apply(self._history[self._cursor])
        elif self._cursor > 0:
            self._cursor -= 1
            self._apply(self._history[self._cursor])
        return self.result

    def redo(self):
        if self._cursor + 1 < len(self._history):
            self._cursor += 1
            self._apply(self._history[self._cursor])
        return self.result

    def reset_tooth(self, tooth):
        state = stateio.clone(self.state)
        state["teeth"][tooth] = stateio.clone(self.baseline["teeth"][tooth])
        self._apply(state)
        self._commit()
        self._persist()
        self._emit()

    def refit(self):
        if not self.evaluate_now()["valid"]:
            raise ValueError("Correct invalid input before refitting.")
        candidate = geometry.refit_state(self.state)
        if not geometry.evaluate(candidate)["valid"]:
            raise ValueError("Refit produced an invalid arch; working shape retained.")
        self._apply(candidate)
        self._commit()
        self._persist()
        self._emit()
        return self.result

    def set_edit_mode(self, mode, tooth=None):
        if not self.parameter_node:
            return
        mode = "arch" if mode == "arch" else "points"
        tooth = tooth or self.parameter_node.GetParameter("SelectedTooth") or next(iter(self.md_nodes))
        self.parameter_node.SetParameter("SelectedTooth", tooth)
        self.parameter_node.SetParameter("EditMode", mode)
        for name, node in self.md_nodes.items():
            node.SetLocked(mode != "points" or name != tooth)
            node.GetDisplayNode().SetPointLabelsVisibility(mode == "points" and name == tooth)
        self.handle_node.SetLocked(mode != "arch")
        self.handle_node.GetDisplayNode().SetVisibility(mode == "arch")
        self.handle_node.GetDisplayNode().SetPointLabelsVisibility(mode == "arch")
        for i in range(self.handle_node.GetNumberOfControlPoints()):
            self.handle_node.SetNthControlPointLocked(i, i in (0, self.handle_node.GetNumberOfControlPoints() - 1))
        self.curve_node.SetLocked(True)

    def set_crown(self, tooth, node):
        if node and node.GetParentTransformNode():
            raise ValueError("Assign a crown without a parent transform.")
        self._reference("Crown_" + tooth, node)
        self._install_observers()
        self.evaluate_now(commit=False)

    def show_baseline(self, visible):
        for role in ["BaselineArch"] + ["BaselineMD_" + t for t in self.md_nodes]:
            node = self.parameter_node.GetNodeReference(role)
            if node:
                node.SetLocked(True)
                node.GetDisplayNode().SetVisibility(visible)

    def save_revision(self, parent_dir):
        if not self.evaluate_now(commit=True)["valid"]:
            raise ValueError("Cannot save a measured revision while inputs are invalid.")
        root = Path(parent_dir).expanduser().resolve()
        if root.suffix.lower() == ".mrb":
            raise ValueError("Choose a directory; revision saving always creates a new file.")
        root.mkdir(parents=True, exist_ok=True)
        revision_id = datetime.datetime.now(datetime.timezone.utc).strftime("Revision_%Y%m%dT%H%M%SZ_") + uuid.uuid4().hex[:8]
        staging = root / ("." + revision_id + ".tmp")
        final = root / revision_id
        staging.mkdir(exist_ok=False)
        provenance = json.loads(self.parameter_node.GetParameter("ProvenanceJSON") or "{}")
        provenance.update(revision_id=revision_id, saved_utc=datetime.datetime.now(datetime.timezone.utc).isoformat())
        previous_url = slicer.mrmlScene.GetURL()
        try:
            self.parameter_node.SetParameter("LastSavedRevision", str(final))
            self._persist()
            stateio.write_exports(staging, self.state, self.baseline, self.result, provenance)
            if not slicer.util.saveScene(str(staging / "DentalArchRevision.mrb")):
                raise IOError("Slicer could not save the revision scene.")
            os.rename(staging, final)
            slicer.mrmlScene.SetURL(str(final / "DentalArchRevision.mrb"))
            return str(final)
        except Exception:
            slicer.mrmlScene.SetURL(previous_url)
            # Keep incomplete exports in the uniquely named staging folder for recovery.
            raise
