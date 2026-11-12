"""Procedurally generated non-anatomical teaching example. No patient data."""
import numpy as np
import vtk
import slicer

def load_demo(logic):
    user_models = [node for node in slicer.util.getNodesByClass("vtkMRMLModelNode") if not node.GetHideFromEditors() and node.GetPolyData() and node.GetPolyData().GetNumberOfCells()]
    user_marks = [node for node in slicer.util.getNodesByClass("vtkMRMLMarkupsNode") if not node.GetHideFromEditors() and node.GetNumberOfControlPoints()]
    if user_models or user_marks:
        raise ValueError("Load the demo in an empty scene (File > Close Scene). Existing work is kept.")
    order = ["VL_P2", "VL_P1", "VL_C", "VL_I2", "VL_I1", "VR_I1", "VR_I2", "VR_C", "VR_P1", "VR_P2"]
    controls = [[-24, 0], [-23, 9], [-15, 23], [0, 30], [15, 23], [23, 9], [24, 0]]
    from scipy.interpolate import BSpline
    knots = [0, 0, 0, 0, .25, .5, .75, 1, 1, 1, 1]
    curve = BSpline(knots, np.array(controls), 3)
    teeth = {}
    for i, (tooth, t) in enumerate(zip(order, np.linspace(.055, .945, len(order)))):
        xy = curve(t)
        tangent = curve(t, 1)
        tangent /= np.linalg.norm(tangent)
        toward_midline = tangent if i < 5 else -tangent
        radius = 3.3
        teeth[tooth] = {"M": [*(xy + radius * toward_midline), 0.], "D": [*(xy - radius * toward_midline), 0.], "kind": "synthetic"}
        source = vtk.vtkSphereSource()
        source.SetCenter(*xy, 0)
        source.SetRadius(radius)
        source.SetThetaResolution(40)
        source.SetPhiResolution(30)
        source.Update()
        model = slicer.modules.models.logic().AddModel(source.GetOutput())
        model.SetName(tooth + "_crown")
        model.GetDisplayNode().SetColor(.78, .80, .84)
    for tooth, x in [("VL_M1", -24), ("VR_M1", 24)]:
        teeth[tooth] = {"M": [x, 0, 0], "D": [x, -7, 0], "kind": "synthetic"}
        source = vtk.vtkSphereSource()
        source.SetCenter(x, -3.5, 0)
        source.SetRadius(3.5)
        source.SetThetaResolution(40)
        source.SetPhiResolution(30)
        source.Update()
        model = slicer.modules.models.logic().AddModel(source.GetOutput())
        model.SetName(tooth + "_crown")
        model.GetDisplayNode().SetColor(.78, .80, .84)
    state = {"schema_version": 1, "frame": {"center": [0, 0, 0], "axes": np.eye(3).tolist(), "plane_z": 0.},
             "arch": {"control_points_xy": controls, "knots": knots, "degree": 3}, "teeth": teeth,
             "included_teeth": order, "boundaries": [{"tooth": "VL_M1", "point": "M"}, {"tooth": "VR_M1", "point": "M"}],
             "units": "synthetic coordinate units", "fit": {"smooth": 1e-4, "huber_fraction": .015}}
    logic.initialize(state, {"source": "Procedurally generated spheres; not anatomical teeth or clinical data"})
    view = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLViewNode")
    if view:
        view.SetBackgroundColor(.055, .075, .105)
        view.SetBackgroundColor2(.055, .075, .105)
        view.SetBoxVisible(False)
        view.SetAxisLabelsVisible(False)
    camera_node = slicer.mrmlScene.GetFirstNodeByClass("vtkMRMLCameraNode")
    if camera_node:
        camera = camera_node.GetCamera()
        camera.SetFocalPoint(0, 11, 0)
        camera.SetPosition(0, 11, 100)
        camera.SetViewUp(0, 1, 0)
        camera.ParallelProjectionOn()
        camera.SetParallelScale(33)
    return logic.result
