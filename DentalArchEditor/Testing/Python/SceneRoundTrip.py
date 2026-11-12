"""Synthetic MRML integration test. Run with Slicer --python-script this_file.

Uses no clinical inputs and writes only to a TemporaryDirectory. BSD-3-Clause.
"""
from pathlib import Path
import json
import sys
import tempfile
import time
import traceback
import numpy as np
import slicer
sys.path.insert(0, str(Path(__file__).resolve().parents[2]))
from DentalArchEditorLib.session import EditorSession
from DentalArchEditorLib import demo, stateio

def run():
    slicer.mrmlScene.Clear(0)
    logic = EditorSession()
    checks = []
    with tempfile.TemporaryDirectory(prefix="DentalArchEditor-test-") as temporary:
        try:
            demo.load_demo(logic)
            assert logic.result["valid"]
            baseline = stateio.clone(logic.state)
            initial_a, initial_r = logic.result["A"], logic.result["R"]
            tooth = baseline["included_teeth"][2]
            p = np.array(logic.state["teeth"][tooth]["M"])
            logic.set_point(tooth, "M", p + [.15, .1, .1])
            assert abs(logic.result["A"] - initial_a) < 1e-9
            assert abs(logic.result["R"] - initial_r) > 1e-5
            checks.append("point-edit owned-surface snapping and manual arch retention")
            before_control = stateio.clone(logic.state)
            c = np.array(logic.state["arch"]["control_points_xy"][3])
            logic.set_arch_control(3, c + [0, .5])
            assert abs(logic.result["A"] - initial_a) > .01
            logic.evaluate_now(commit=True)
            logic.undo()
            assert np.allclose(logic.state["arch"]["control_points_xy"], before_control["arch"]["control_points_xy"], atol=1e-10, rtol=0)
            logic.redo()
            checks.append("control-edit, numerical-noop history, undo and redo")
            first = Path(logic.save_revision(temporary))
            second = Path(logic.save_revision(temporary))
            assert first != second and (first / "DentalArchRevision.mrb").is_file()
            saved = stateio.clone(logic.state)
            saved_c = logic.result["C"]
            slicer.mrmlScene.Clear(0)
            assert slicer.util.loadScene(str(second / "DentalArchRevision.mrb"))
            assert logic.parameter_node and logic.result["valid"]
            assert abs(logic.result["C"] - saved_c) < 1e-8
            assert np.allclose(logic.state["arch"]["control_points_xy"], saved["arch"]["control_points_xy"], atol=1e-10)
            assert np.allclose(logic.baseline["arch"]["control_points_xy"], baseline["arch"]["control_points_xy"], atol=1e-10)
            checks.append("unique revision files, immutable baseline and self-contained MRB reload")
            old_r = logic.result["R"]
            node = logic.md_nodes[tooth]
            p = np.array(node.GetNthControlPointPositionVector(0))
            node.SetNthControlPointPosition(0, *(p + [.3, .2, .15]))
            time.sleep(.08)
            slicer.app.processEvents()
            assert logic.result["valid"] and abs(logic.result["R"] - old_r) > 1e-5
            checks.append("live observer reattachment after reload")
            node.UnsetNthControlPointPosition(0)
            time.sleep(.08)
            slicer.app.processEvents()
            assert not logic.result["valid"] and logic.result["C"] is None
            logic.undo()
            assert logic.result["valid"]
            checks.append("undefined-point invalidation and undo recovery")
            try:
                logic.save_revision(first / "DentalArchRevision.mrb")
                raise AssertionError("Existing scene path was accepted")
            except ValueError:
                pass
            checks.append("original scene path rejected as save target")
            return {"passed": True, "checks": checks}
        finally:
            logic.cleanup()
            slicer.mrmlScene.Clear(0)

if __name__ == "__main__":
    try:
        print(json.dumps(run(), indent=2), flush=True)
        slicer.util.exit(0)
    except Exception:
        traceback.print_exc()
        slicer.util.exit(1)
