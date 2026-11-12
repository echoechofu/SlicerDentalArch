"""Generate reviewable numeric results from a source-bound reviewed case config.

Use run_in_slicer.py to additionally create an editable Slicer MRB scene.
"""
from pathlib import Path
import argparse
import hashlib
import json
import time

from case_io import load_case, PREMOLARS
import detect_cusps
import detect_incisal
import mesiodistal_algorithm
import reference_arch_algorithm
import calculate_crowding


def run(source, config_path, out, sensitivity=True):
    source, out = Path(source).resolve(), Path(out).resolve()
    # Validate before creating output files. Never replace an existing run.
    mesh = load_case(source, config_path)
    config = mesh['config']
    del mesh
    known = ['pipeline_manifest.json', 'resolved_config.json', 'landmarks', 'incisal',
             'mesiodistal', 'arch', 'crowding_results.json']
    if any((out / name).exists() for name in known):
        raise FileExistsError('Output contains an earlier or incomplete run. Choose a new output directory; baseline results are never overwritten.')
    out.mkdir(parents=True, exist_ok=True)
    manifest = dict(schema_version=1, status='running', source_sha256=config['sha256'],
                    units='native STL coordinate units', sensitivity_enabled=bool(sensitivity),
                    source_unchanged=False, stages={}, review_required=True,
                    scope='Reviewed 14-crown upper arch; provisional geometric candidates, not automatic tooth segmentation or clinical ground truth.')
    manifest_path = out / 'pipeline_manifest.json'
    resolved_path = out / 'resolved_config.json'

    def save_manifest():
        manifest_path.write_text(json.dumps(manifest, indent=2), encoding='utf8')

    def stage(name, operation, result_path):
        start = time.monotonic()
        print('Running ' + name, flush=True)
        result = operation()
        manifest['stages'][name] = dict(status='complete', seconds=round(time.monotonic() - start, 3), result=result_path)
        save_manifest()
        return result

    save_manifest()
    try:
        landmarks = stage('landmarks', lambda: detect_cusps.run(out / 'landmarks', source, config, previews=False),
                          'landmarks/landmark_candidates.json')
        anchors = dict(config.get('groove_anchors', {}))
        provenance = {}
        for name in PREMOLARS:
            if name in anchors:
                provenance[name] = 'supplied_in_reviewed_configuration'
                continue
            cusps = [point for point in landmarks['points'] if point['region'] == name and point['kind'] == 'cusp_candidate']
            if len(cusps) != 2:
                raise ValueError(f'{name}: need two distinct cusp candidates to initialize the groove corridor. Review and supply groove_anchors in a new configuration/run.')
            anchors[name] = [point['aligned_xyz'][:2] for point in cusps]
            provenance[name] = 'two_detected_cusp_candidates_provisional_requires_anatomical_review'
        config['groove_anchors'] = anchors
        config['groove_anchor_provenance'] = provenance
        resolved_path.write_text(json.dumps(config, indent=2), encoding='utf8')
        stage('incisal', lambda: detect_incisal.run(source, out / 'incisal', resolved_path, sensitivity=sensitivity),
              'incisal/incisal_results.json')
        stage('mesiodistal', lambda: mesiodistal_algorithm.run(source, out / 'mesiodistal', sensitivity=sensitivity, config_path=resolved_path),
              'mesiodistal/mesiodistal_results.json')
        stage('arch', lambda: reference_arch_algorithm.run(out / 'mesiodistal/mesiodistal_results.json', out / 'arch', sensitivity=sensitivity),
              'arch/reference_arch_results.json')
        stage('crowding', lambda: calculate_crowding.run(out / 'mesiodistal/mesiodistal_results.json', out / 'arch/reference_arch_results.json', out),
              'crowding_results.json')
        if hashlib.sha256(source.read_bytes()).hexdigest() != config['sha256']:
            raise RuntimeError('Source file changed during computation; discard this run and inspect the source again.')
        manifest.update(status='complete', source_unchanged=True, resolved_config='resolved_config.json')
        save_manifest()
    except Exception as error:
        manifest.update(status='failed', error=str(error))
        save_manifest()
        raise
    return manifest


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--input', type=Path, required=True)
    parser.add_argument('--config', type=Path, required=True)
    parser.add_argument('--output-dir', type=Path, required=True)
    parser.add_argument('--skip-sensitivity', action='store_true')
    args = parser.parse_args()
    run(args.input, args.config, args.output_dir, sensitivity=not args.skip_sensitivity)
