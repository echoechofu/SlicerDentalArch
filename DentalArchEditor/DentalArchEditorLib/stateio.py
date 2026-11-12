"""Portable, case-independent state import and export. BSD-3-Clause."""
import copy
import csv
import hashlib
import json
from pathlib import Path
import numpy as np

VERSION = "0.1.0"

def dumps(value):
    return json.dumps(value, ensure_ascii=False, allow_nan=False, separators=(",", ":"), default=lambda x: x.tolist() if hasattr(x, "tolist") else float(x))

def import_baseline(md_path, arch_path):
    md_raw = Path(md_path).read_bytes()
    md = json.loads(md_raw)
    arch = json.loads(Path(arch_path).read_text(encoding="utf-8-sig"))
    if md.get("source_sha256") != arch.get("source_STL_sha256"):
        raise ValueError("Measurement and arch files refer to different sources.")
    if arch.get("input_v2_json_sha256") and hashlib.sha256(md_raw).hexdigest() != arch["input_v2_json_sha256"]:
        raise ValueError("The arch was fitted from a different version of the measurement file.")
    if not np.allclose(md["pca_axes_columns"], arch["pca_axes_columns"], atol=1e-10) or not np.allclose(md["pca_center"], arch["pca_center"], atol=1e-10):
        raise ValueError("Source coordinate frames disagree.")
    teeth = {}
    for tooth in md["teeth"]:
        if tooth["id"] in teeth:
            raise ValueError("Duplicate tooth identifier.")
        teeth[tooth["id"]] = {"M": tooth["mesial_native"], "D": tooth["distal_native"], "kind": tooth.get("kind", "unknown")}
    ends = []
    for name in arch["endpoint_ids"]:
        tooth, point = name.rsplit("_", 1)
        ends.append({"tooth": tooth, "point": point})
    state = {
        "schema_version": 1,
        "frame": {"center": arch["pca_center"], "axes": arch["pca_axes_columns"], "plane_z": arch["reference_plane_aligned_z"]},
        "arch": {"control_points_xy": arch["control_points_xy"], "knots": arch["knots"], "degree": arch["degree"]},
        "teeth": teeth, "included_teeth": arch["input_tooth_ids"], "boundaries": ends,
        "units": "native coordinate units", "fit": {"smooth": arch.get("smoothing_lambda", 1e-4), "huber_fraction": .015},
    }
    return state, {"expected_A": float(arch["projected_arch_length"]), "source_fingerprint": md.get("source_sha256", ""), "algorithm_version": md.get("algorithm_version", "unknown")}

def write_exports(directory, state, baseline, result, provenance):
    directory = Path(directory)
    payload = {"format": "DentalArchEditorRevision", "schema_version": 1, "editor_version": VERSION,
               "baseline": baseline, "current": state,
               "measurements": {k: result[k] for k in ("R", "A", "C", "D", "widths", "warnings")},
               "provenance": provenance}
    (directory / "revision.json").write_text(json.dumps(json.loads(dumps(payload)), ensure_ascii=False, indent=2), encoding="utf8")
    with (directory / "measurements.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(["tooth", "included_in_R", "M_x", "M_y", "M_z", "D_x", "D_y", "D_z", "width_native_units"])
        for tooth, points in state["teeth"].items():
            writer.writerow([tooth, tooth in state["included_teeth"], *points["M"], *points["D"], result["widths"][tooth]])
        writer.writerow([])
        writer.writerow(["A", result["A"], "available reference arc length"])
        writer.writerow(["R", result["R"], "required sum of included tooth widths"])
        writer.writerow(["C=R-A", result["C"], "positive shortage, negative surplus"])
        writer.writerow(["D=A-R", result["D"], "positive surplus, negative shortage"])
    with (directory / "arch.csv").open("w", newline="", encoding="utf-8-sig") as stream:
        writer = csv.writer(stream)
        writer.writerow(["sample", "x_native", "y_native", "z_native"])
        for i, point in enumerate(result["curve_native"]):
            writer.writerow([i, *point])

def clone(state):
    return copy.deepcopy(state)
