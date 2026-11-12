# Contributing

This project welcomes reproducible fixes and improvements to landmark algorithms,
measurement editing, geometry, accessibility, documentation, and interoperability.

## Scope

The algorithm pipeline consumes reviewed orientation and tooth-to-component maps;
the module edits its scene and explicit measurement protocol. Keep algorithm
changes in `algorithms/` separate from the editor's measurement state and UI.
Do not describe component discovery as segmentation of arbitrary fused scans.
Keep the numerical core importable without Slicer; place MRML and widget code in
the Slicer layer. Avoid machine-specific paths and tooth mappings tied to one case.

## Reproduce changes

Describe the observed behavior, the triggering operation, and the expected result.
Add a small synthetic regression case when a change fixes geometry, calculation,
state restoration or file handling. Run the numerical tests documented in the
README, then verify affected Slicer interactions in the target version.

Relevant invariants include:

- Reviewed configurations are bound to the exact source STL and reject incomplete maps.
- Algorithm runs preserve the source and refuse to overwrite earlier output.
- Import preserves the baseline measurements and curve.
- Each M/D pair stays associated with the intended tooth and crown.
- Manual point edits do not silently refit an already edited arch.
- The reference plane and measurement range remain explicit and consistent.
- Saving a revision does not overwrite the imported source or another revision.
- Reload restores values and reconnects live calculation.
- Invalid or incomplete state is reported rather than shown as a valid current result.

## Data and reports

Public changes, tests, issues and screenshots must use synthetic data or data
explicitly cleared for public redistribution. The default demo is generated from
spheres and mathematical curves. Do not include patient names, scan identifiers,
source hashes, local clinical paths, real coordinate sets, or clinical files in a
public report. JSON, CSV, logs and screenshots can expose data even without an STL.
Keep private validation outside the repository. Review the staged diff and history;
the supplied `.gitignore` only reduces accidental additions.

## Dependencies and licenses

Contributions to project code are accepted under the repository's BSD-3-Clause
license. Preserve existing copyright notices. Document any new runtime dependency,
and include required notices for any copied code or asset. Do not copy unrelated
implementation into the project without checking its license and attribution.

Slicer supplies the GUI and scientific runtime. Avoid downloading or installing
dependencies automatically during routine module use. Keep tests independent of
network access and external services.

## Release preparation

Before a public release, verify the selected source files contain no private
validation data, run the documented tests, check a clean local installation, and
state exactly which Slicer version and operating systems were exercised. Review
the LICENSE copyright attribution and any dependency changes with maintainers.
Populate repository and extension metadata only after their URLs exist.

An Extensions Index submission is a separate distribution step; see the
[Slicer documentation](https://slicer.readthedocs.io/en/latest/developer_guide/extensions.html).
Do not describe the extension as indexed or clinically validated until that status
has actually been established.
