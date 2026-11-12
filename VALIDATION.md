# Validation record

Checked locally on Windows with **3D Slicer 5.12.4**. These checks establish
software behavior and compatibility, not anatomical localization accuracy or
clinical diagnostic validation.

## Public automated checks

- **22 tests passed** across the numerical editor core and pipeline configuration
  tests. They cover curve/measurement behavior and input rejection, including
  source-bound configuration, required review flags, tooth maps and mesh checks.
- The synthetic `DentalArchEditor/Testing/Python/SceneRoundTrip.py` integration
  test passed: owned-surface edits, preservation of manually edited curves,
  control-point changes, numerical-noop undo handling, separate revision files,
  immutable baseline data, self-contained MRB reload, automatic event
  reconnection, undefined-point invalidation and recovery, and rejection of an
  existing MRB filename as a revision directory.
- The module was discovered through Additional module paths. Its Qt panel,
  synthetic demo, mode selector, baseline comparison, live measurement table,
  and undo/redo buttons were exercised.

Run the public checks using the commands in the [README](README.md). The
procedural demo and synthetic tests contain no clinical data and do not test
anatomical landmark accuracy.

## Private compatibility checks

Two private local cases were rerun with the generalized reviewed-configuration
pipeline. Their numerical results matched the earlier algorithm baseline.

A fresh full batch run created `DentalArchScene.mrb` and successfully reloaded
its editor state. The scene included 14 crown models and 27 checked algorithm
annotation nodes, plus editor nodes.
The source mesh was unchanged, and exported/reloaded coordinates agreed to
floating-point precision (maximum observed difference approximately `7.3e-15`
in native coordinate units).

Existing annotated-scene compatibility checks also covered initial values,
mapped-crown snapping, boundary linkage, node renaming, explicit refit, invalid
curve and transform handling, revision save/reload, and source-file hashes.
Private meshes, configurations, scenes, coordinates and reports are excluded
from the public repository.

These small development checks are not an independent validation dataset,
a multi-operator reproducibility study, or a measure of diagnostic performance.

## Not yet established

- Landmark accuracy against independently annotated anatomical ground truth.
- Robustness on different scanners, coordinate scales, tooth shapes, missing
  teeth, lower arches, or fused unsegmented scans.
- macOS/Linux GUI behavior or a CMake-built extension package.
- Remote GitHub Actions results are pending the first source push; the 22-test
  result above is from the local run.
- Installation through the Slicer Extensions Index.
- Clinical accuracy or clinical severity classification.

See the [algorithm guide](docs/algorithms.md) for the exact supported input,
measurement protocol and fixed-scale detection assumptions. See the
[editor guide](docs/editor.md) for editing and revision behavior.

