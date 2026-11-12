# DentalArchEditor for 3D Slicer

[中文说明](editor.zh-CN.md) · [Project overview](../README.md)

DentalArchEditor is a scripted Slicer module for correcting mesiodistal crown-width
landmarks and a planar reference arch while viewing the resulting space analysis.
It keeps the imported algorithm baseline separate from the editable revision and
saves revisions to new files.

This is an experimental measurement editor. It consumes an annotated scene and
measurement metadata; it does not automatically identify or segment teeth in an
arbitrary STL. It has not been validated for clinical diagnostic accuracy. All
lengths remain in the input coordinate units unless those units are independently
established.

The repository is distributed under BSD-3-Clause. It is not currently a
listing in the Slicer Extensions Index. No patient data is included.

## Install locally

Target environment: **3D Slicer 5.12.4**. The module uses the NumPy, SciPy, VTK and
Qt/PythonQt environment supplied by Slicer. A system Python installation alone
cannot run its graphical interface.

1. Extract the source archive into a stable local directory.
2. In Slicer, open **Edit → Application Settings → Modules**.
3. Add the directory containing `DentalArchEditor.py` to **Additional module paths**.
   This is the `DentalArchEditor` subdirectory, not the repository's top-level directory.
4. Restart Slicer and select **Dentistry → Dental Arch Editor / 牙弓修订**
   from the module selector.

Python-only modules can be loaded this way without compiling Slicer or the
extension. See [Slicer's extension documentation](https://slicer.readthedocs.io/en/latest/developer_guide/extensions.html).

## Try the synthetic demo

Use **加载合成示例 / Synthetic demo** to create a mathematical arch and simple spherical
crown models. The demo lets a new user exercise point editing, arch editing,
undo/redo, live measurements, and revision saving without any clinical files.
Its shapes and measurements are artificial and are not anatomical examples.

## Open a generated scene

The [algorithm pipeline](algorithms.md) creates `DentalArchScene.mrb` with the
editor state already initialized. Open that MRB and this module to start reviewing
it; do not import its JSON files again. MRB scenes contain data, not running Python
callbacks, so the module must be installed and opened to resume live calculation.

## Import an older annotated scene

Open the annotated Slicer scene, then use **Import baseline** and select its two
associated metadata files:

- Mesiodistal results JSON: tooth IDs, original M/D coordinates, coordinate frame,
  and the source identity.
- Reference-arch results JSON: the original B-spline controls, knots, degree,
  reference-plane height, measurement range and boundary IDs.

The initial importer understands the measurement-result schemas used by this
project. For other pipelines, provide matching metadata or an adapter; an arbitrary
point list is not enough to establish tooth identities or the measurement range.
The existing 501-point display polyline is not treated as 501 independent arch
editing handles. Import preserves the analytic spline, rather than refitting it
silently.

The current tooth can be associated with its own crown model through the module's
crown selector. Confirm this association before using surface snapping. The labels
`VL` and `VR` describe the input convention; the module does not infer patient
laterality or FDI numbering from those labels.

For an adapter, the importer expects these fields:

| File | Required measurement fields |
| --- | --- |
| MD JSON | `source_sha256`, `pca_center`, `pca_axes_columns`, and a `teeth` list containing `id`, `mesial_native`, `distal_native` and optional `kind` |
| Arch JSON | `source_STL_sha256`, the same frame fields, `reference_plane_aligned_z`, `control_points_xy`, `knots`, `degree`, `input_tooth_ids`, `endpoint_ids`, and `projected_arch_length` |

Both source identifiers must agree. When `input_v2_json_sha256` is present, it must
match the exact MD JSON bytes. An endpoint ID has the form `toothID_M` or `toothID_D`.
Frame axes are columns of a right-handed orthonormal matrix: native coordinates
map to the local frame as `(point - pca_center) @ pca_axes_columns`.
The current importer discovers baseline line nodes named `toothID_MD` or
`toothID_MD_CANDIDATE` and crown models named `toothID_crown` or
`toothID_isolated_crown`. Duplicate matches are rejected. Newly saved revisions
use persistent node references rather than relying on those display names.

## Edit and measure

1. Select a tooth and choose **Edit M-D points**. Drag its M or D point in the 3D view.
   With snapping enabled, the point is projected to that tooth's associated crown.
   The current distance between its two points becomes its edited width.
2. Choose **Edit arch handles** to move the internal B-spline control handles in the fixed
   reference plane. Control handles are not necessarily on the visible curve.
3. Read required space `R`, available space `A`, signed crowding `C = R - A`, and
   the baseline/current comparison in the side panel.
4. Use **Undo**, **Redo**, **Reset tooth**, or the baseline comparison
   display to review changes.
5. Use **Refit arch** explicitly when the reference curve should be recomputed
   from the current measurement points.

In manual editing, changing an internal tooth's M/D points changes `R` and keeps
the internal arch controls unchanged. The two curve endpoints remain bound to the
projected first-molar mesial landmarks: moving a boundary point also changes `A`.
The reference frame stays fixed during editing. Refit is a separate operation so
that a manually adjusted arch is not silently replaced by moving another tooth's
measurement point.

The supported curve is a planar cubic B-spline. Its progress along the
boundary-to-boundary direction must be monotone; the numerical core checks each
polynomial span. This prevents self-intersections but deliberately rejects some
non-self-intersecting curves that double back transversely. Reverse bending is a
review warning. Invalid inputs clear the current totals instead of retaining a
previous value as if it were current. Parent transforms are not an alternate way
to change the measurement frame; keep imported landmarks and crown models in
their shared native coordinate system.

Posterior upper-third and direction rules from an automatic initialization are
review aids. Manually edited points are not silently moved back to the algorithm's
maximum-chord solution. Inspect rule deviations and any geometry/status messages
before accepting a revision.

## Measurement definition

For the default upper-arch protocol:

```text
R = sum of ||M_i - D_i|| over the 10 teeth from P2 to contralateral P2
A = integral of the planar reference spline's speed between the two M1 mesial boundaries
C = R - A       positive: net shortage; negative: net surplus
D = A - R       positive: net surplus; negative: net shortage
```

First molars define the boundaries; their widths are not included in `R`. Second
molars are outside this default space range. A small net surplus does not exclude
local rotation, overlap or displacement. These landmarks represent intrinsic
crown-width measurement locations, not necessarily existing physical contacts
between neighboring teeth. The reference arch is a measurement curve, not a
treatment target.

## Save a revision and reopen it

Use the module's **Save new revision** action. Choose a destination directory and
keep each revision's scene, JSON state and CSV measurements together. Existing
revision targets and the imported original file are not overwritten by this action.

The revision contains both the baseline and the edited state, including point
coordinates, analytic spline controls, reference frame, measurement range and
result values. The scene contains the models needed for viewing. Open the installed
Dental Arch Editor module and load a revision scene to resume interactive calculation. An MRB stores
data, not running Python callbacks; the module reconnects its observers after loading.

The module's save protection applies to its revision-saving action. Slicer's
general-purpose **Save** dialog remains a separate application feature; do not
use it to intentionally replace your original source file.

Saving a revision stores the current calculation. Editing a CSV outside Slicer
does not update the scene. Keep unmodified algorithm outputs as the source record.

## Development and testing

The numerical core is independent of Slicer and uses synthetic test cases:

```sh
python -m pip install numpy scipy vtk
python -m unittest discover -s tests -p "test_*.py" -v
```

The GitHub Actions workflow runs this numerical suite. Passing it does not validate
mouse interaction, scene lifecycle or clinical measurement accuracy. Slicer
integration checks must additionally exercise dragging, surface ownership,
undo/redo, endpoint binding, explicit refit, save/reload, and baseline preservation.

A distributable, synthetic MRML integration test is also included:

```sh
Slicer --no-main-window --no-splash --python-script DentalArchEditor/Testing/Python/SceneRoundTrip.py
```

It verifies editing, undo, unique saves, source preservation and automatic observer
restoration after loading a saved scene. It creates only temporary synthetic data.
See [VALIDATION.md](../VALIDATION.md) for the locally checked release scope.

The CMake files support Slicer's extension build process. Direct module-path
installation is the primary local workflow; package builds and cross-platform
distribution must be verified in their target environments. Before an Extensions
Index submission, verify the extension metadata, add any required public icon
and documentation URLs, and follow the
[official distribution process](https://slicer.readthedocs.io/en/latest/developer_guide/extensions.html).

See [CONTRIBUTING.md](../CONTRIBUTING.md) for contribution and data-handling rules and
[THIRD_PARTY_NOTICES.md](../THIRD_PARTY_NOTICES.md) for runtime licenses and method
references. The project license is [BSD-3-Clause](../LICENSE).

