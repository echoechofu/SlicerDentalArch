"""Batch entry point: reviewed STL configuration -> editor-ready MRB scene.

Slicer --no-main-window --no-splash --python-script algorithms/run_in_slicer.py -- \
  --input model.stl --config reviewed.json --output-dir output/session
"""
import argparse
import json
from pathlib import Path
import sys
import traceback

sys.path.insert(0, str(Path(__file__).resolve().parent))


def main(argv=None):
    import slicer
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--input", required=True, type=Path)
    parser.add_argument("--config", required=True, type=Path)
    parser.add_argument("--output-dir", required=True, type=Path)
    parser.add_argument("--skip-sensitivity", action="store_true")
    args = parser.parse_args(argv)
    out = args.output_dir.resolve()
    if out.exists() and any(out.iterdir()):
        raise FileExistsError("Use a new or empty output directory to preserve earlier algorithm results.")
    out.mkdir(parents=True, exist_ok=True)
    from run_all import run
    from build_scene import run as build_scene
    run(args.input, args.config, out, sensitivity=not args.skip_sensitivity)
    scene = build_scene(args.input, out / "resolved_config.json", out)
    (out / "pipeline_status.json").write_text(json.dumps({"success": True, "scene": scene.name}, indent=2), encoding="utf8")
    print("Created " + str(scene), flush=True)
    return 0


if __name__ == "__main__":
    import slicer
    arguments = sys.argv[sys.argv.index("--") + 1:] if "--" in sys.argv else sys.argv[1:]
    try:
        status = main(arguments)
    except SystemExit as error:
        status = int(error.code or 0)
    except Exception:
        traceback.print_exc()
        status = 1
    slicer.util.exit(status)
