"""Build a portable, editor-ready Slicer scene from pipeline results. BSD-3-Clause.

Run only inside 3D Slicer. Native STL XYZ values are preserved numerically;
an anatomical RAS/LPS orientation or physical millimetre scale is not inferred.
"""
import json
from pathlib import Path
import sys

import numpy as np
import vtk
import slicer
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray, vtk_to_numpy

from case_io import load_case, tooth_region


def polydata(vertices, triangles):
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(np.ascontiguousarray(vertices, dtype=np.float64), deep=True))
    faces = np.column_stack((np.full(len(triangles), 3, dtype=np.int64), triangles)).astype(np.int64)
    cells = vtk.vtkCellArray()
    cells.SetCells(len(triangles), numpy_to_vtkIdTypeArray(faces.ravel(), deep=True))
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetPolys(cells)
    return mesh


def run(source, config_path, output_dir):
    """Create DentalArchScene.mrb in a fresh Slicer scene; refuse overwriting it."""
    out = Path(output_dir).resolve()
    target = out / "DentalArchScene.mrb"
    if target.exists():
        raise FileExistsError("The algorithm scene already exists; use a new output directory.")
    # This is a batch operation. Never clear a user's loaded model or editor session.
    user_models = [n for n in slicer.util.getNodesByClass("vtkMRMLModelNode") if not n.GetHideFromEditors()]
    if user_models or slicer.util.getNodesByClass("vtkMRMLMarkupsNode"):
        raise RuntimeError("Use a fresh Slicer process for scene generation.")
    def read(relative):
        return json.loads((out / relative).read_text(encoding="utf-8-sig"))
    md_path = out / "mesiodistal/mesiodistal_results.json"
    arch_path = out / "arch/reference_arch_results.json"
    md = read("mesiodistal/mesiodistal_results.json")
    arch = read("arch/reference_arch_results.json")
    inc = read("incisal/incisal_results.json")
    landmarks = read("landmarks/landmark_candidates.json")
    crowding = read("crowding_results.json")
    mesh = load_case(source, config_path)
    digest = mesh["config"]["sha256"]
    for label, record, key in (("mesiodistal", md, "source_sha256"),
                               ("landmarks", landmarks, "source_sha256"),
                               ("incisal", inc, "source_sha256"),
                               ("arch", arch, "source_STL_sha256"),
                               ("crowding", crowding, "source_sha256")):
        if record.get(key) != digest:
            raise ValueError("Scene source and " + label + " fingerprints disagree.")

    manager = slicer.app.layoutManager()
    if manager:
        manager.setLayout(slicer.vtkMRMLLayoutNode.SlicerLayoutOneUp3DView)
    if not slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLViewNode"):
        slicer.mrmlScene.AddNewNodeByClass("vtkMRMLViewNode")
    hierarchy = slicer.vtkMRMLSubjectHierarchyNode.GetSubjectHierarchyNode(slicer.mrmlScene)
    folders = {name: hierarchy.CreateFolderItem(hierarchy.GetSceneItemID(), name) for name in
               ("01_Crowns", "02_Cusps", "03_Central_fossae", "04_Incisal_edges",
                "05_Baseline_MD", "06_Groove_references", "07_Baseline_arch", "08_Algorithm_results", "09_Editable_measurements")}
    def group(node, folder):
        hierarchy.SetItemParent(hierarchy.GetItemByDataNode(node), folders[folder])
    original = slicer.modules.models.logic().AddModel(polydata(mesh["v"], mesh["f"]))
    original.SetName("Original_STL_native_coordinates")
    original.GetDisplayNode().SetVisibility(False)
    original.SetAttribute("SourceSHA256", mesh["config"]["sha256"])
    original.SetAttribute("Coordinates", "Raw native STL XYZ; no RAS/LPS reflection applied")
    crown_counts = {}
    for region in mesh["config"]["regions"]:
        ids, face_ids = tooth_region(mesh, region)
        remap = np.full(len(mesh["v"]), -1, dtype=np.int64)
        remap[ids] = np.arange(len(ids))
        crown = slicer.modules.models.logic().AddModel(polydata(mesh["v"][ids], remap[mesh["f"][face_ids]]))
        crown.SetName(region["id"] + "_crown")
        crown.GetDisplayNode().SetColor(.77, .80, .84)
        group(crown, "01_Crowns")
        crown_counts[crown.GetName()] = len(face_ids)

    expected = {}
    def mark(name, node_class, points, labels, color, folder, visible=False, glyph=.35):
        node = slicer.mrmlScene.AddNewNodeByClass(node_class, name)
        node.CreateDefaultDisplayNodes()
        if node.IsA("vtkMRMLMarkupsCurveNode"):
            node.SetCurveTypeToLinear()
            node.SetNumberOfPointsPerInterpolatingSegment(1)
        for i, point in enumerate(points):
            index = node.AddControlPoint(vtk.vtkVector3d(*[float(x) for x in point]))
            node.SetNthControlPointLabel(index, labels[i] if labels else name + "_" + str(i))
        display = node.GetDisplayNode()
        display.SetColor(*color)
        display.SetSelectedColor(*color)
        display.SetPointLabelsVisibility(False)
        display.SetPropertiesLabelVisibility(False)
        display.SetUseGlyphScale(False)
        display.SetGlyphSize(glyph)
        display.SetVisibility(visible)
        for i in range(node.GetNumberOfMeasurements()):
            node.GetNthMeasurement(i).SetEnabled(False)
        node.SetLocked(True)
        node.SetAttribute("ReviewRequired", "True")
        node.SetAttribute("Units", "native coordinate units; millimetre calibration not assumed")
        group(node, folder)
        expected[name] = np.asarray(points, dtype=float).reshape(-1, 3)
        return node

    for kind, name, color, folder in (
        ("cusp_candidate", "Cusp_candidates", (1, .36, .26), "02_Cusps"),
        ("central_fossa_candidate", "Central_fossa_candidates", (.36, .65, 1), "03_Central_fossae"),
    ):
        points = [p for p in landmarks["points"] if p["kind"] == kind]
        node = mark(name, "vtkMRMLMarkupsFiducialNode", [p["native_xyz"] for p in points],
                    [p["label"] for p in points], color, folder)
        for i, point in enumerate(points):
            node.SetNthControlPointDescription(i, json.dumps({"triangle": point["triangle_index"], "review_flags": point.get("review_flags", [])}))
    incisal_points, incisal_labels = [], []
    for curve in inc["curves"]:
        mark(curve["region"] + "_incisal_edge", "vtkMRMLMarkupsCurveNode", curve["points_native"],
             None, (1, .86, .27), "04_Incisal_edges", glyph=.04)
        for key, label in (("endpoint_A_native", "A"), ("midpoint_native", "Mid"), ("endpoint_B_native", "B")):
            incisal_points.append(curve[key])
            incisal_labels.append(curve["region"] + "_" + label)
    mark("Incisal_endpoints_and_midpoints", "vtkMRMLMarkupsFiducialNode", incisal_points,
         incisal_labels, (1, .86, .27), "04_Incisal_edges")
    for tooth in md["teeth"]:
        node = mark(tooth["id"] + "_MD", "vtkMRMLMarkupsLineNode",
                    [tooth["mesial_native"], tooth["distal_native"]], [tooth["id"] + "_M", tooth["id"] + "_D"],
                    (.25, .96, .62), "05_Baseline_MD")
        node.SetAttribute("InitialWidthNativeUnits", str(tooth["width"]))
        node.SetAttribute("ReviewFlags", "; ".join(tooth["review_flags"]))
        if tooth["kind"] == "premolar":
            mark(tooth["id"] + "_groove", "vtkMRMLMarkupsCurveNode",
                 tooth["axis_diagnostics"]["groove"]["points_native"], None, (.89, .46, .95), "06_Groove_references", glyph=.04)
    mark("Algorithm_reference_arch", "vtkMRMLMarkupsCurveNode", arch["curve_samples_native"],
         None, (.14, .87, .96), "07_Baseline_arch", glyph=.03)
    mark("Reference_M1_mesial_boundaries", "vtkMRMLMarkupsFiducialNode", arch["endpoint_projected_native"],
         arch["endpoint_ids"], (.14, .87, .96), "07_Baseline_arch")

    table = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTableNode", "Algorithm_measurements_native_units")
    rows = [(t["id"], str(t["width"]), "Included in R" if t["id"] in crowding["range_tooth_ids"] else "Outside R") for t in md["teeth"]]
    rows += [("A", str(crowding["available_space_A"]), "Available planar reference arc"),
             ("R", str(crowding["required_space_R"]), "Required ten-tooth width sum"),
             ("C=R-A", str(crowding["signed_crowding_R_minus_A"]), "Positive shortage; negative surplus")]
    for index, label in enumerate(("Quantity", "Native_coordinate_units", "Meaning")):
        column = vtk.vtkStringArray()
        column.SetName(label)
        for row in rows:
            column.InsertNextValue(row[index])
        table.AddColumn(column)
    group(table, "08_Algorithm_results")
    note = slicer.mrmlScene.AddNewNodeByClass("vtkMRMLTextNode", "READ_ME")
    note.SetText("Algorithm candidates require review. Select DentalArchEditor to edit contact points and arch handles; the panel recalculates C=R-A. Baseline annotations are locked and initially hidden. Enable them in Data to inspect cusps, fossae and incisal curves. Save a revision from the module to create a new directory without overwriting this scene. Coordinates are native STL values; no millimetre calibration or anatomical RAS/LPS orientation is implied. The saved algorithm table is a baseline snapshot, not the live panel.")
    group(note, "08_Algorithm_results")

    sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "DentalArchEditor"))
    from DentalArchEditorLib.session import EditorSession
    session = EditorSession()
    session.import_baseline(md_path, arch_path)
    for node in [*session.md_nodes.values(), session.curve_node, session.handle_node]:
        group(node, "09_Editable_measurements")
    camera_node = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLCameraNode") or slicer.mrmlScene.AddNewNodeByClass("vtkMRMLCameraNode")
    view = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLViewNode")
    camera_node.SetActiveTag(view.GetID())
    camera = camera_node.GetCamera()
    focus = np.asarray(mesh["center"])
    extent = float(np.ptp(mesh["q"], axis=0).max())
    camera.SetFocalPoint(*focus)
    camera.SetPosition(*(focus + 3 * extent * mesh["axes"][:, 2]))
    camera.SetViewUp(*mesh["axes"][:, 1])
    camera.ParallelProjectionOn()
    camera.SetParallelScale(.60 * extent)
    view.SetBackgroundColor(.055, .075, .105)
    view.SetBackgroundColor2(.055, .075, .105)
    view.SetBoxVisible(False)
    view.SetAxisLabelsVisible(False)
    if not slicer.util.saveScene(str(target)):
        raise RuntimeError("Slicer failed to save the scene bundle.")
    session.cleanup()
    slicer.mrmlScene.Clear(0)
    if not slicer.util.loadScene(str(target)):
        raise RuntimeError("Saved scene could not be reloaded.")
    errors = []
    for name, points in expected.items():
        node = slicer.mrmlScene.GetFirstNodeByName(name)
        actual = np.array([node.GetNthControlPointPositionVector(i) for i in range(node.GetNumberOfControlPoints())]).reshape(-1, 3)
        if actual.shape != points.shape or not np.allclose(actual, points, rtol=0, atol=1e-8):
            raise RuntimeError("Markup round-trip mismatch: " + name)
        errors.append(float(np.linalg.norm(actual - points, axis=1).max()) if len(points) else 0.)
    for name, count in crown_counts.items():
        if slicer.mrmlScene.GetFirstNodeByName(name).GetPolyData().GetNumberOfCells() != count:
            raise RuntimeError("Crown triangle count changed during saving.")
    reloaded = slicer.mrmlScene.GetFirstNodeByName("Original_STL_native_coordinates").GetPolyData()
    raw_error = float(np.linalg.norm(vtk_to_numpy(reloaded.GetPoints().GetData()) - mesh["v"], axis=1).max())
    if raw_error > 1e-8 or reloaded.GetNumberOfCells() != len(mesh["f"]):
        raise RuntimeError("Raw source mesh changed during saving.")
    restored = EditorSession()
    try:
        if not restored.result["valid"] or abs(restored.result["C"] - crowding["signed_crowding_R_minus_A"]) > 1e-6:
            raise RuntimeError("Editor state did not reproduce algorithm measurements after scene reload.")
        validation = {"saved_and_reloaded": True, "editor_ready": True,
                      "max_coordinate_roundtrip_error": max(errors), "original_mesh_coordinate_error": raw_error,
                      "crown_count": len(crown_counts), "markup_nodes_checked": len(expected),
                      "A": restored.result["A"], "R": restored.result["R"], "C": restored.result["C"]}
        (out / "scene_validation.json").write_text(json.dumps(validation, indent=2), encoding="utf8")
    finally:
        restored.cleanup()
    return target
