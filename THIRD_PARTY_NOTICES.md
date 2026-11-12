# Third-party software and method references

SlicerDentalArch project code, including DentalArchEditor and the algorithm pipeline,
is distributed under the [BSD-3-Clause license](LICENSE).
It uses the host application's APIs and runtime libraries. The source archive
does not bundle Slicer, NumPy, SciPy, VTK, Qt, CTK, or PythonQt binaries or source.
Their respective licenses remain in force; the extension's license does not
relicense them.

| Runtime dependency | Purpose | Upstream license information |
| --- | --- | --- |
| 3D Slicer | MRML scenes, markups, module lifecycle, Python environment | [Slicer Contribution and Software License Agreement](https://github.com/Slicer/Slicer/blob/main/License.txt), a BSD-style agreement with additional terms |
| NumPy | Vector and matrix calculations | [NumPy BSD license](https://numpy.org/doc/stable/license.html) |
| SciPy | B-splines, integration and numerical optimization | [SciPy BSD-3-Clause license](https://github.com/scipy/scipy/blob/main/LICENSE.txt) |
| VTK | Mesh geometry and surface queries | [VTK license](https://vtk.org/about/#license), BSD-style |
| Qt / PythonQt / CTK | User interface and application integration supplied by Slicer | Use the notices shipped with the installed Slicer distribution; these projects have separate license terms. |

The small CMake files use the documented Slicer extension/module build API.
The module and numerical helpers are project code, not a copy of Slicer core
implementation. If future contributions copy or bundle upstream implementation
or assets, retain their notices and update this file.

## Method background

The project implements its own explicit measurement protocol: an intrinsic
mesiodistal width sum compared with a planar, smooth reference curve. The
following references provide background, not clinical validation of this code:

- [BeGole, 1979: cubic spline dental arch analysis](https://pubmed.ncbi.nlm.nih.gov/487781/).
- [AlHarbi et al., 2008: mathematical analyses of dental arch curvature](https://pubmed.ncbi.nlm.nih.gov/18251619/).
- [Lang et al., 2025: digital three-dimensional cast analysis](https://link.springer.com/article/10.1186/s13005-025-00515-8). Its segmented perimeter measurement is not identical to this editor's continuous reference-curve arc length.

No patient mesh, scan, landmark coordinate set, case screenshot, or other clinical
dataset is licensed or distributed as part of this source package. The demo is
generated mathematically from simple curves and spheres.
