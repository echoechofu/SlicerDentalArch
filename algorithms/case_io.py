"""Load a binary STL only with an explicitly reviewed, source-bound tooth map.

The pipeline consumes separate, closed crown components; it does not segment a
fused scan. Aligned coordinates are ``(native - center) @ axes_columns``.
"""
from pathlib import Path
import copy
import hashlib
import json
import re

import numpy as np
from scipy.sparse import coo_matrix
from scipy.sparse.csgraph import connected_components

TOOTH_TYPES = (('I1', 'central_incisor'), ('I2', 'lateral_incisor'),
               ('C', 'canine'), ('P1', 'premolar'), ('P2', 'premolar'),
               ('M1', 'molar'), ('M2', 'molar'))
EXPECTED_REGIONS = [(side + '_' + tooth, kind)
                    for side in ('VL', 'VR') for tooth, kind in TOOTH_TYPES]
PREMOLARS = [name for name, kind in EXPECTED_REGIONS if kind == 'premolar']


def read_binary_stl(source):
    """Read triangles without repairing, merging by tolerance, or changing units."""
    raw = Path(source).read_bytes()
    if len(raw) < 84:
        raise ValueError('Expected a nonempty binary STL with an 80-byte header.')
    count = int.from_bytes(raw[80:84], 'little')
    if count == 0 or len(raw) != 84 + 50 * count:
        raise ValueError('Only binary STL with an exact triangle record count is supported. Export an unmodified binary STL first.')
    dtype = np.dtype([('normal', '<f4', (3,)), ('points', '<f4', (3, 3)), ('attr', '<u2')])
    triangles = np.frombuffer(raw, dtype=dtype, offset=84, count=count)['points'].astype(float)
    if not np.isfinite(triangles).all():
        raise ValueError('STL contains non-finite vertex coordinates.')
    vertices, inverse = np.unique(triangles.reshape(-1, 3), axis=0, return_inverse=True)
    return vertices, inverse.reshape(-1, 3), hashlib.sha256(raw).hexdigest()


def face_components(faces):
    """Stable groups connected across edges shared by exactly two triangles.

    Group IDs are tied to the input triangle order and therefore to its SHA-256.
    Non-manifold junctions are not silently treated as crown boundaries: each
    selected group is subsequently required to be a closed manifold.
    """
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    owner = np.tile(np.arange(len(faces)), 3)
    _, inverse, count = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    order = np.argsort(inverse, kind='stable')
    starts = np.r_[0, np.cumsum(count)][:-1][count == 2]
    a, b = owner[order[starts]], owner[order[starts + 1]]
    graph = coo_matrix((np.ones(len(a)), (a, b)), shape=(len(faces), len(faces))).tocsr()
    return connected_components(graph, directed=False)[1]


def topology(vertices, faces):
    """Report topology and winding; never fix a mesh as a side effect."""
    points = vertices[np.unique(faces)]
    triangles = vertices[faces]
    cross = np.cross(triangles[:, 1] - triangles[:, 0], triangles[:, 2] - triangles[:, 0])
    scale = max(float(np.linalg.norm(np.ptp(points, axis=0))), np.finfo(float).tiny)
    zero_area = int(np.count_nonzero(np.linalg.norm(cross, axis=1) <= scale * scale * 1e-14))
    edges = np.concatenate([faces[:, [0, 1]], faces[:, [1, 2]], faces[:, [2, 0]]])
    _, inverse, counts = np.unique(np.sort(edges, axis=1), axis=0, return_inverse=True, return_counts=True)
    winding = np.bincount(inverse, weights=np.where(edges[:, 0] < edges[:, 1], 1, -1))
    origin = points.mean(0)
    shifted = triangles - origin
    volume = float(np.einsum('ij,ij->i', shifted[:, 0], np.cross(shifted[:, 1], shifted[:, 2])).sum() / 6)
    return dict(boundary_edges=int(np.count_nonzero(counts == 1)),
                nonmanifold_edges=int(np.count_nonzero(counts > 2)),
                inconsistent_winding_edges=int(np.count_nonzero((counts == 2) & (winding != 0))),
                degenerate_triangles=zero_area, signed_volume=volume,
                closed=bool(np.all(counts == 2)), outward=bool(volume > scale ** 3 * 1e-14))


def proposed_alignment(vertices):
    """An unreviewed PCA proposal, with arbitrary signs, for inspection only."""
    center = vertices.mean(0)
    _, axes = np.linalg.eigh(np.cov((vertices - center).T))
    axes = axes[:, ::-1]
    if np.linalg.det(axes) < 0:
        axes[:, 0] *= -1
    return center, axes


def validate_config(config):
    config = copy.deepcopy(config)
    if config.get('schema_version', 1) != 1:
        raise ValueError('Unsupported case configuration schema_version; expected 1.')
    for flag in ('orientation_reviewed', 'tooth_mapping_reviewed'):
        if config.get(flag) is not True:
            raise ValueError(f'{flag} must be true after human review; the inspection template is not runnable.')
    if not isinstance(config.get('sha256'), str) or not re.fullmatch('[0-9a-fA-F]{64}', config['sha256']):
        raise ValueError('sha256 must be the 64-character digest of the reviewed source STL.')
    config['sha256'] = config['sha256'].lower()
    alignment = config.get('alignment', {})
    center = np.asarray(alignment.get('center'), dtype=float)
    axes = np.asarray(alignment.get('axes_columns'), dtype=float)
    if center.shape != (3,) or axes.shape != (3, 3) or not np.isfinite(center).all() or not np.isfinite(axes).all():
        raise ValueError('alignment needs a finite center[3] and axes_columns[3][3].')
    if not np.allclose(axes.T @ axes, np.eye(3), atol=1e-10, rtol=0) or not np.isclose(np.linalg.det(axes), 1, atol=1e-10, rtol=0):
        raise ValueError('alignment.axes_columns must be orthonormal and right-handed; do not use scale or shear.')
    regions = config.get('regions', [])
    if not isinstance(regions, list) or any(not isinstance(r, dict) for r in regions):
        raise ValueError('regions must be a list of reviewed tooth-to-group mappings.')
    expected = dict(EXPECTED_REGIONS)
    if len(regions) != 14 or len({r.get('id') for r in regions}) != 14 or {r.get('id') for r in regions} != set(expected):
        raise ValueError('Exactly the 14 expected VL/VR I1,I2,C,P1,P2,M1,M2 regions are required; partial arches are unsupported.')
    for region in regions:
        if region.get('kind') != expected[region['id']]:
            raise ValueError(f"Wrong tooth kind for {region['id']}.")
        group = region.get('group')
        if isinstance(group, bool) or not isinstance(group, int) or group < 0:
            raise ValueError(f"Assign a nonnegative integer component group to {region['id']}.")
    if len({r['group'] for r in regions}) != 14:
        raise ValueError('Each tooth must map to a different independent crown component.')
    by_id = {r['id']: r for r in regions}
    config['regions'] = [by_id[name] for name, _ in EXPECTED_REGIONS]
    anchors = config.get('groove_anchors', {})
    if not isinstance(anchors, dict) or set(anchors) - set(PREMOLARS):
        raise ValueError('groove_anchors may only map premolar IDs to two aligned XY cusp anchors.')
    for name, value in anchors.items():
        points = np.asarray(value, dtype=float)
        if points.shape != (2, 2) or not np.isfinite(points).all() or np.linalg.norm(points[0] - points[1]) < 1e-8:
            raise ValueError(f'{name} groove anchors must be two distinct finite aligned XY points.')
    config['schema_version'] = 1
    return config


def load_case(source, config_path):
    """Return native vertices/faces and the reviewed frame/map after validation."""
    if source is None or config_path is None:
        raise ValueError('Both source STL and reviewed case config are required.')
    config = config_path if isinstance(config_path, dict) else json.loads(Path(config_path).read_text(encoding='utf-8-sig'))
    config = validate_config(config)
    vertices, faces, digest = read_binary_stl(source)
    if digest != config['sha256']:
        raise ValueError('Source STL SHA-256 differs from the reviewed configuration. Inspect and review this exact file.')
    groups = face_components(faces)
    for region in config['regions']:
        selected = faces[groups == region['group']]
        if not len(selected):
            raise ValueError(f"Component {region['group']} for {region['id']} is absent from this mesh.")
        report = topology(vertices, selected)
        if (not report['closed'] or not report['outward'] or report['inconsistent_winding_edges'] or report['degenerate_triangles']):
            raise ValueError(f"{region['id']} component {region['group']} is not a closed, outward, consistently oriented nondegenerate crown: {report}. This pipeline does not segment fused scans or repair meshes.")
    center = np.asarray(config['alignment']['center'], dtype=float)
    axes = np.asarray(config['alignment']['axes_columns'], dtype=float)
    config['source'] = str(Path(source).resolve())
    return dict(v=vertices, f=faces, q=(vertices - center) @ axes, center=center,
                axes=axes, face_groups=groups, config=config)


def tooth_region(mesh, config):
    face_ids = np.flatnonzero(mesh['face_groups'] == config['group'])
    return np.unique(mesh['f'][face_ids]), face_ids
