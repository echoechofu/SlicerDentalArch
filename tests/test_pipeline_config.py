"""Synthetic-only tests for source binding and reviewed crown-map requirements."""
from pathlib import Path
import copy
import hashlib
import struct
import sys
import tempfile
import unittest

import numpy as np

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / 'algorithms'))
import case_io


def synthetic_mesh():
    # Four outward faces of a tetrahedron, repeated as 14 disconnected objects.
    vertices = np.array([[0., 0., 0.], [1., 0., 0.], [0., 1., 0.], [0., 0., 1.]])
    faces = np.array([[0, 2, 1], [0, 1, 3], [0, 3, 2], [1, 2, 3]])
    return np.concatenate([vertices[faces] + [3 * i, 0, 0] for i in range(14)])


def binary_stl(triangles):
    records = [struct.pack('<12fH', *([0.] * 3 + tri.flatten().tolist()), 0) for tri in triangles]
    return b'Synthetic geometry only'.ljust(80, b' ') + struct.pack('<I', len(records)) + b''.join(records)


class PipelineConfigTests(unittest.TestCase):
    def setUp(self):
        self.tmp = tempfile.TemporaryDirectory()
        self.source = Path(self.tmp.name) / 'synthetic.stl'
        self.raw = binary_stl(synthetic_mesh())
        self.source.write_bytes(self.raw)
        self.config = dict(schema_version=1, sha256=hashlib.sha256(self.raw).hexdigest(),
                           alignment=dict(center=[3., 2., 1.], axes_columns=[[0., 1., 0.], [-1., 0., 0.], [0., 0., 1.]]),
                           orientation_reviewed=True, tooth_mapping_reviewed=True,
                           regions=[dict(id=name, kind=kind, group=i) for i, (name, kind) in enumerate(case_io.EXPECTED_REGIONS)])

    def tearDown(self):
        self.tmp.cleanup()

    def test_supplied_frame_and_original_source_are_preserved(self):
        mesh = case_io.load_case(self.source, self.config)
        np.testing.assert_array_equal(mesh['center'], self.config['alignment']['center'])
        np.testing.assert_array_equal(mesh['axes'], self.config['alignment']['axes_columns'])
        np.testing.assert_allclose(mesh['q'] @ mesh['axes'].T + mesh['center'], mesh['v'])
        self.assertEqual(self.source.read_bytes(), self.raw)
        self.assertNotIn('source', self.config)

    def test_region_order_is_canonical(self):
        config = copy.deepcopy(self.config)
        config['regions'].reverse()
        mesh = case_io.load_case(self.source, config)
        self.assertEqual([r['id'] for r in mesh['config']['regions']], [r[0] for r in case_io.EXPECTED_REGIONS])

    def test_unreviewed_template_cannot_run(self):
        for flag in ('orientation_reviewed', 'tooth_mapping_reviewed'):
            config = copy.deepcopy(self.config)
            config[flag] = False
            with self.subTest(flag=flag), self.assertRaisesRegex(ValueError, flag):
                case_io.load_case(self.source, config)

    def test_different_source_is_rejected(self):
        self.source.write_bytes(b'X' + self.raw[1:])
        with self.assertRaisesRegex(ValueError, 'SHA-256 differs'):
            case_io.load_case(self.source, self.config)

    def test_missing_duplicate_and_unassigned_regions_are_rejected(self):
        for mode in ('missing', 'duplicate', 'unassigned', 'absent', 'wrong_kind'):
            config = copy.deepcopy(self.config)
            if mode == 'missing': config['regions'].pop()
            if mode == 'duplicate': config['regions'][1]['group'] = 0
            if mode == 'unassigned': config['regions'][0]['group'] = None
            if mode == 'absent': config['regions'][0]['group'] = 99
            if mode == 'wrong_kind': config['regions'][0]['kind'] = 'molar'
            with self.subTest(mode=mode), self.assertRaises(ValueError):
                case_io.load_case(self.source, config)

    def test_scale_shear_and_reflection_are_rejected(self):
        for axes in (np.eye(3) * 2, [[1, 1, 0], [0, 1, 0], [0, 0, 1]], np.diag([-1, 1, 1])):
            config = copy.deepcopy(self.config)
            config['alignment']['axes_columns'] = np.asarray(axes).tolist()
            with self.subTest(axes=axes), self.assertRaisesRegex(ValueError, 'orthonormal and right-handed'):
                case_io.load_case(self.source, config)

    def test_open_inward_and_inconsistent_crowns_are_rejected(self):
        for mode in ('open', 'inward', 'inconsistent'):
            triangles = synthetic_mesh()
            if mode == 'open': triangles = triangles[1:]
            if mode == 'inward': triangles[:4] = triangles[:4, ::-1]
            if mode == 'inconsistent': triangles[0] = triangles[0, ::-1]
            raw = binary_stl(triangles)
            self.source.write_bytes(raw)
            config = copy.deepcopy(self.config)
            config['sha256'] = hashlib.sha256(raw).hexdigest()
            with self.subTest(mode=mode), self.assertRaisesRegex(ValueError, 'closed, outward'):
                case_io.load_case(self.source, config)

    def test_ascii_truncated_and_nonfinite_inputs_are_rejected(self):
        invalid = [b'solid synthetic\nendsolid', self.raw[:-1]]
        triangles = synthetic_mesh()
        triangles[0, 0, 0] = np.nan
        invalid.append(binary_stl(triangles))
        for raw in invalid:
            self.source.write_bytes(raw)
            with self.subTest(size=len(raw)), self.assertRaises(ValueError):
                case_io.read_binary_stl(self.source)

    def test_invalid_groove_anchors_are_rejected(self):
        for anchors in ({'VL_M1': [[0, 0], [1, 1]]}, {'VL_P1': [[0, 0], [0, 0]]}, {'VL_P1': [[0, 0, 1], [1, 1, 2]]}):
            config = copy.deepcopy(self.config)
            config['groove_anchors'] = anchors
            with self.subTest(anchors=anchors), self.assertRaises(ValueError):
                case_io.load_case(self.source, config)


if __name__ == '__main__':
    unittest.main()
