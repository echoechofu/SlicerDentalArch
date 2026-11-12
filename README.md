# SlicerDentalArch

[中文说明](README.zh-CN.md)

SlicerDentalArch combines dental landmark algorithms with a **3D Slicer scene and
interactive measurement editor**. It proposes cusps, molar central fossae, four
incisal edges, mesiodistal width points and a reference arch. DentalArchEditor then
lets a user adjust width points and the arch, see space analysis update, and save a
separate revision while preserving the algorithm baseline.

```text
STL → inspect components and review tooth mapping / orientation
    → landmark and reference-arch algorithms → DentalArchScene.mrb
    → DentalArchEditor → manual correction → new revision MRB + JSON + CSV
```

This is an experimental, semi-automatic workflow. Its input must already contain
separable crown mesh components; it does not segment arbitrary fused dental scans.
Orientation and tooth identities require human review. Software checks do not
establish anatomical or clinical accuracy. Measurements use the original STL
coordinate units; an STL alone does not establish millimetres.

Project code is BSD-3-Clause. No patient data is included. This repository is a
source distribution, not an entry in the Slicer Extensions Index.

## Requirements

- Target environment: **3D Slicer 5.12.4**, tested locally on Windows.
- Slicer's Python runtime supplies NumPy, SciPy, VTK and Qt/PythonQt.
- A binary STL with 14 individually separable upper crown components, plus any
  unrelated base components. Each selected crown must be closed, consistently
  outward oriented and nondegenerate. See [input requirements](docs/algorithms.md).
- Coordinate scale similar to a dental model represented in millimetre-sized
  coordinates. Detection uses fixed distance parameters: inspect bounds and
  parameters before applying it to a different scale. No automatic rescaling is
  performed.

## From STL to an editable Slicer scene

### 1. Install the editor

Clone or download this repository. In Slicer, open **Edit → Application Settings →
Modules → Additional module paths**, add the repository's `DentalArchEditor`
directory, then restart Slicer. Select **Dentistry → Dental Arch Editor / 牙弓修订**.

The [editor guide](docs/editor.md) explains installation, controls, saved revisions
and importing an existing annotated scene.

### 2. Inspect the mesh and prepare a reviewed configuration

Run with Slicer's `PythonSlicer` executable, replacing the executable and input
paths with your own. The following commands are PowerShell examples run from the
repository root:

```powershell
& "C:\path\to\Slicer\bin\PythonSlicer.exe" algorithms/inspect_stl.py --input model.stl --output-dir work/inspect
```

Inspect the component reports and labelled overview images. Copy
`work/inspect/config.template.json` to a private `reviewed.json`, map the 14 tooth
labels to the correct component IDs, and review the alignment frame. Set the two
review flags to `true` only after checking them. The initial PCA frame is a
proposal, not an anatomical orientation decision.

See the [algorithm guide](docs/algorithms.md) for the complete configuration
contract, limitations and troubleshooting. Store clinical inputs and outputs
outside the public source tree.

### 3. Run the algorithms and create the scene

```powershell
& "C:\path\to\Slicer\Slicer.exe" --ignore-slicerrc --disable-settings --no-main-window --no-splash --python-script algorithms/run_in_slicer.py -- --input model.stl --config reviewed.json --output-dir output/session
```

The main deliverable is **`output/session/DentalArchScene.mrb`**, a self-contained
Slicer scene with the crown model, landmarks, reference arch, preserved baseline
and editable measurement state. JSON/CSV intermediate results and a pipeline
manifest are also written. Use a fresh output directory for each run. Optional
`--skip-sensitivity` skips additional incisal, width and reference-arch sensitivity runs.

To calculate JSON/CSV results without creating a scene:

```powershell
& "C:\path\to\Slicer\bin\PythonSlicer.exe" algorithms/run_all.py --input model.stl --config reviewed.json --output-dir output/session
```

### 4. Review, correct and save

Open the generated MRB and select DentalArchEditor. It is already initialized;
there is no need to import the algorithm JSON again. Drag the M/D points or the
reference-arch controls and inspect the live required length, available length and
their difference. Internal tooth-point edits do not silently refit the arch.

Use the editor's **Save new revision** operation to create a new revision directory
with an MRB, JSON and CSV. Both the baseline and working result are retained.
Reopen the revision with the module installed to reconnect live calculations.
Slicer's general Save dialog is separate and can overwrite a chosen file.

## What is measured

The default analysis spans from the projected mesial point of one first molar to
the other. The width sum covers **10 teeth, second premolar to second premolar**.

| Quantity | Definition |
| --- | --- |
| Required length `R` | Sum of the 10 intrinsic 3D M/D chord lengths |
| Available length `A` | Analytic arc length of the planar cubic B-spline reference arch |
| Crowding `C = R − A` | Positive: shortage; negative: surplus |
| Space difference `D = A − R` | The opposite sign convention |

The reference arch smooths local displacement while following the current overall
arch form; it is not a treatment target. The first-molar widths are excluded and
second molars lie outside the measurement span. A net surplus does not rule out
local crowding or rotation. No clinical severity category is assigned.

## Demo, testing and development

**加载合成示例 / Synthetic demo** in the module generates a simple mathematical scene
without patient data. It demonstrates editing and saving; it is not an accuracy
benchmark for tooth landmark detection.

```powershell
python -m pip install numpy scipy vtk
python -m unittest discover -s tests -v
& "C:\path\to\Slicer\Slicer.exe" --no-main-window --no-splash --python-script DentalArchEditor/Testing/Python/SceneRoundTrip.py
```

See [validation and its limits](VALIDATION.md), [contributing](CONTRIBUTING.md),
[third-party notices](THIRD_PARTY_NOTICES.md) and the [license](LICENSE).

| Directory | Contents |
| --- | --- |
| `algorithms/` | Inspection, reviewed configuration, landmark/width/arch computation, scene creation |
| `DentalArchEditor/` | Slicer module, numerical editor core and synthetic demonstration |
| `docs/` | Algorithm and editor usage guides in English and Chinese |
| `tests/` | Numerical and Slicer integration tests |


