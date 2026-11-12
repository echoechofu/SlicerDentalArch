"""Inspect binary-STL face components and write an UNREVIEWED mapping template.

PCA is a geometric proposal, not anatomical orientation. Component IDs identify
topology, not tooth classes. Review both before enabling the algorithm pipeline.
"""
from pathlib import Path
import argparse
import csv
import json

import numpy as np
import vtk
from vtk.util.numpy_support import numpy_to_vtk, numpy_to_vtkIdTypeArray

from case_io import EXPECTED_REGIONS, read_binary_stl, face_components, proposed_alignment, topology


def preview(q, faces, groups, components, path, bottom=False):
    points = vtk.vtkPoints()
    points.SetData(numpy_to_vtk(q, deep=True))
    cells = vtk.vtkCellArray()
    cells.SetData(numpy_to_vtkIdTypeArray(np.arange(0, 3 * len(faces) + 1, 3, dtype=np.int64), deep=True),
                  numpy_to_vtkIdTypeArray(faces.astype(np.int64).ravel(), deep=True))
    mesh = vtk.vtkPolyData()
    mesh.SetPoints(points)
    mesh.SetPolys(cells)
    # Deterministic colors distinguish topological groups only.
    palette = np.random.default_rng(10).integers(80, 230, (int(groups.max()) + 1, 3), dtype=np.uint8)
    colors = numpy_to_vtk(palette[groups], deep=True, array_type=vtk.VTK_UNSIGNED_CHAR)
    colors.SetName('Face component')
    mesh.GetCellData().SetScalars(colors)
    mapper = vtk.vtkPolyDataMapper()
    mapper.SetInputData(mesh)
    mapper.SetScalarModeToUseCellData()
    actor = vtk.vtkActor()
    actor.SetMapper(mapper)
    actor.GetProperty().SetAmbient(.35)
    renderer = vtk.vtkRenderer()
    renderer.SetBackground(.06, .075, .10)
    renderer.AddActor(actor)
    sign = -1 if bottom else 1
    # Label the 32 largest groups; every group remains available in the tables.
    for component in sorted(components, key=lambda c: c['face_count'], reverse=True)[:32]:
        center = np.array(component['center_aligned'])
        center[2] = q[:, 2].max() + 1 if sign > 0 else q[:, 2].min() - 1
        label = vtk.vtkBillboardTextActor3D()
        label.SetInput(str(component['group']))
        label.SetPosition(center)
        label.GetTextProperty().SetFontSize(23)
        label.GetTextProperty().SetColor(1, 1, 1)
        label.GetTextProperty().SetBold(True)
        label.GetTextProperty().SetBackgroundColor(0, 0, 0)
        label.GetTextProperty().SetBackgroundOpacity(.65)
        renderer.AddActor(label)
    caption = vtk.vtkTextActor()
    caption.SetInput(('PCA -Z view' if bottom else 'PCA +Z view') + ' | UNREVIEWED | numbers = face groups')
    caption.SetPosition(18, 18)
    caption.GetTextProperty().SetFontSize(21)
    renderer.AddActor2D(caption)
    camera = renderer.GetActiveCamera()
    center = (q.min(0) + q.max(0)) / 2
    camera.SetFocalPoint(center)
    camera.SetPosition(center + np.array([0., 0., sign * max(np.ptp(q, axis=0).max(), 1.) * 3]))
    camera.SetViewUp(0, 1, 0)
    camera.ParallelProjectionOn()
    renderer.ResetCamera()
    window = vtk.vtkRenderWindow()
    window.SetOffScreenRendering(1)
    window.SetSize(1500, 1300)
    window.AddRenderer(renderer)
    window.Render()
    snapshot = vtk.vtkWindowToImageFilter()
    snapshot.SetInput(window)
    snapshot.ReadFrontBufferOff()
    snapshot.Update()
    writer = vtk.vtkPNGWriter()
    writer.SetFileName(str(path))
    writer.SetInputConnection(snapshot.GetOutputPort())
    writer.Write()
    window.Finalize()


def run(source, output_dir):
    out = Path(output_dir)
    out.mkdir(parents=True, exist_ok=True)
    if (out / 'config.template.json').exists():
        raise FileExistsError('Inspection template already exists. Choose a new output directory to preserve any review edits.')
    vertices, faces, digest = read_binary_stl(source)
    groups = face_components(faces)
    center, axes = proposed_alignment(vertices)
    q = (vertices - center) @ axes
    components = []
    for group in range(int(groups.max()) + 1):
        selected = faces[groups == group]
        ids = np.unique(selected)
        components.append(dict(group=group, face_count=len(selected), vertex_count=len(ids),
                               center_native=vertices[ids].mean(0).tolist(), center_aligned=q[ids].mean(0).tolist(),
                               bounds_aligned=[q[ids].min(0).tolist(), q[ids].max(0).tolist()],
                               topology=topology(vertices, selected)))
    report = dict(source_sha256=digest, vertex_count=len(vertices), face_count=len(faces),
                  component_count=len(components), components=components,
                  component_definition='Faces connected across edges incident to exactly two triangles; IDs depend on source triangle order.',
                  limitations=['A component is not an automatically recognized tooth.',
                               'Fused gum/tooth scans need external segmentation and a new review.',
                               'PCA axis order and signs require review; +qZ must face the occlusal surfaces and +qY the anterior arch.',
                               'Preview labels show the 32 largest groups; all groups are in this report.'])
    (out / 'components.json').write_text(json.dumps(report, indent=2), encoding='utf8')
    with (out / 'components.csv').open('w', newline='', encoding='utf-8-sig') as stream:
        writer = csv.writer(stream)
        writer.writerow(['group', 'face_count', 'vertex_count', 'qX', 'qY', 'qZ',
                         'closed', 'outward', 'boundary_edges', 'nonmanifold_edges',
                         'inconsistent_winding_edges', 'degenerate_triangles'])
        for item in components:
            t = item['topology']
            writer.writerow([item['group'], item['face_count'], item['vertex_count'], *item['center_aligned'],
                             *[t[k] for k in ('closed', 'outward', 'boundary_edges', 'nonmanifold_edges',
                                             'inconsistent_winding_edges', 'degenerate_triangles')]])
    config = dict(schema_version=1, sha256=digest,
                  alignment=dict(center=center.tolist(), axes_columns=axes.tolist()),
                  regions=[dict(id=name, group=None, kind=kind) for name, kind in EXPECTED_REGIONS],
                  orientation_reviewed=False, tooth_mapping_reviewed=False,
                  groove_anchors={},
                  notes='UNREVIEWED template. +qZ must face occlusally, +qY anteriorly; VL/VR are the left/right sides of the reviewed +qZ preview, not patient laterality. Review each group and preserve a right-handed orthonormal frame. The displayed group centroids use this proposed frame.')
    (out / 'config.template.json').write_text(json.dumps(config, indent=2), encoding='utf8')
    preview(q, faces, groups, components, out / 'inspection_top.png')
    preview(q, faces, groups, components, out / 'inspection_bottom.png', bottom=True)
    print(json.dumps(dict(components=len(components), template=str(out / 'config.template.json'), reviewed=False)), flush=True)
    return report


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    args = parser.parse_args()
    run(args.input, args.output_dir)
