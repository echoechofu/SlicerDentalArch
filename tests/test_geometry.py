# SPDX-License-Identifier: BSD-3-Clause
"""Synthetic-only numerical tests; runnable without Slicer or patient files."""

from copy import deepcopy
from pathlib import Path
import sys
import unittest

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import BSpline

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "DentalArchEditor" / "DentalArchEditorLib"))
from geometry import (arc_length, basis, evaluate, fit_arch, footpoints,
                      refit_state, synchronize_endpoints)


def synthetic_state():
    teeth = {}
    included = ["Tooth" + str(i) for i in range(10)]
    for i, name in enumerate(included):
        positions = [[i + 0.15, 0.0, 2.0], [i + 0.85, 0.0, 2.0]]
        distal, mesial = positions if i < 5 else positions[::-1]
        teeth[name] = {"M": mesial, "D": distal, "kind": "synthetic"}
    teeth["FirstBoundary"] = {"M": [0, 0, 1], "D": [-2, 0, 1]}
    teeth["LastBoundary"] = {"M": [10, 0, 4], "D": [12, 0, 4]}
    return {
        "schema_version": 1,
        "frame": {"center": [0, 0, 0], "axes": np.eye(3).tolist(), "plane_z": 2.0},
        "arch": {"control_points_xy": [[0, 0], [10 / 3, 0], [20 / 3, 0], [10, 0]],
                 "knots": [0, 0, 0, 0, 1, 1, 1, 1], "degree": 3},
        "teeth": teeth, "included_teeth": included,
        "boundaries": [{"tooth": "FirstBoundary", "point": "M"},
                       {"tooth": "LastBoundary", "point": "M"}],
        "fit": {"smooth": 1e-4, "huber_fraction": 0.015},
        "units": "synthetic coordinate units",
    }


class GeometryTests(unittest.TestCase):
    def test_signed_space_and_scope(self):
        state = synthetic_state()
        result = evaluate(state)
        self.assertTrue(result["valid"], result["errors"])
        self.assertAlmostEqual(result["R"], 7.0, places=12)
        self.assertAlmostEqual(result["A"], 10.0, places=12)
        self.assertAlmostEqual(result["C"], -3.0, places=12)
        self.assertAlmostEqual(result["D"], 3.0, places=12)
        # Boundary teeth are measurable but their widths are outside the R scope.
        state["teeth"]["LastBoundary"]["D"] = [20, 0, 4]
        changed = evaluate(state)
        self.assertAlmostEqual(changed["widths"]["LastBoundary"], 10)
        self.assertAlmostEqual(changed["R"], result["R"])
        self.assertAlmostEqual(changed["A"], result["A"])

    def test_width_edit_keeps_manual_arch(self):
        state = synthetic_state()
        controls = deepcopy(state["arch"])
        state["teeth"]["Tooth4"]["M"][0] += 0.2
        before = deepcopy(state)
        result = evaluate(state)
        self.assertTrue(result["valid"])
        self.assertEqual(state, before)
        self.assertEqual(state["arch"], controls)
        self.assertAlmostEqual(result["R"], 7.2)
        self.assertAlmostEqual(result["A"], 10)

    def test_boundary_update_keeps_interior_controls_and_fixed_plane(self):
        state = synthetic_state()
        original_interior = deepcopy(state["arch"]["control_points_xy"][1:-1])
        state["teeth"]["LastBoundary"]["M"] = [11, 0, 9]
        original = deepcopy(state)
        result = evaluate(state)
        self.assertEqual(state, original)
        self.assertTrue(result["valid"], result["errors"])
        self.assertEqual(result["controls_xy"][1:-1], original_interior)
        self.assertEqual(result["controls_native"][-1], [11, 0, 2])
        self.assertAlmostEqual(result["A"], 11)
        self.assertAlmostEqual(result["R"], 7)
        synchronize_endpoints(state)
        self.assertEqual(state["arch"]["control_points_xy"][-1], [11, 0])
        self.assertEqual(state["arch"]["control_points_xy"][1:-1], original_interior)
        self.assertEqual(state["frame"], original["frame"])

    def test_rigid_frame_invariance_and_native_coordinate_output(self):
        state = synthetic_state()
        angle = 0.71
        rotation = np.array([[np.cos(angle), 0, np.sin(angle)],
                             [0, 1, 0], [-np.sin(angle), 0, np.cos(angle)]])
        shift = np.array([12.2, -8.7, 31.1])
        expected = evaluate(state)
        for tooth in state["teeth"].values():
            for role in ("M", "D"):
                tooth[role] = (rotation @ tooth[role] + shift).tolist()
        state["frame"]["center"] = shift.tolist()
        state["frame"]["axes"] = rotation.tolist()
        result = evaluate(state)
        self.assertTrue(result["valid"], result["errors"])
        for key in ("R", "A", "C", "D"):
            self.assertAlmostEqual(result[key], expected[key], places=10)
        np.testing.assert_allclose(result["curve_native"],
                                   np.array(expected["curve_native"]) @ rotation.T + shift,
                                   atol=1e-12)

    def test_height_difference_affects_width_but_not_projected_arch(self):
        state = synthetic_state()
        state["teeth"]["Tooth0"]["M"][2] += 1
        result = evaluate(state)
        self.assertTrue(result["valid"])
        self.assertAlmostEqual(result["R"], 6.3 + np.sqrt(1.49))
        self.assertAlmostEqual(result["A"], 10)

    def test_invalid_input_clears_all_totals(self):
        mutations = [
            lambda s: s["teeth"]["Tooth0"].pop("M"),
            lambda s: s["teeth"]["Tooth0"].update(M=[np.nan, 0, 0]),
            lambda s: s["teeth"]["Tooth0"].update(M=s["teeth"]["Tooth0"]["D"]),
            lambda s: s["included_teeth"].append("Tooth0"),
            lambda s: s["frame"].update(axes=(2 * np.eye(3)).tolist()),
            lambda s: s["arch"].update(knots=[0, 0, 0, 0, 1, 1, 1]),
            lambda s: s["arch"].update(control_points_xy=[[0, 0], [1, 1]]),
            lambda s: s["arch"].update(degree=2),
            lambda s: s["teeth"]["LastBoundary"].update(M=[0, 0, 99]),
        ]
        for mutation in mutations:
            with self.subTest(mutation=mutation):
                state = synthetic_state()
                mutation(state)
                result = evaluate(state)
                self.assertFalse(result["valid"])
                self.assertTrue(result["errors"])
                for key in ("R", "A", "C", "D"):
                    self.assertIsNone(result[key])
        for state in [None, [], "invalid", {"schema_version": 2}]:
            result = evaluate(state)
            self.assertFalse(result["valid"])
            self.assertIsNone(result["C"])

    def test_foldback_is_invalid_even_with_coarse_display_sampling(self):
        state = synthetic_state()
        state["arch"]["control_points_xy"] = [[0, 0], [40, 8], [-30, 8], [10, 0]]
        result = evaluate(state, samples=2)
        self.assertFalse(result["valid"])
        self.assertFalse(result["shape_diagnostics"]["no_self_intersection_certified"])
        self.assertIsNone(result["C"])
        self.assertEqual(len(result["curve_native"]), 2)

    def test_stationary_curve_is_invalid(self):
        state = synthetic_state()
        # x(t)=40(t-.5)^3+5 has a zero-speed point inside the span.
        state["arch"]["control_points_xy"] = [[0, 0], [10, 0], [0, 0], [10, 0]]
        result = evaluate(state)
        self.assertFalse(result["valid"])
        self.assertFalse(result["shape_diagnostics"]["nonzero_speed_certified"])

    def test_reverse_bending_is_warning_not_invalid(self):
        state = synthetic_state()
        state["arch"]["control_points_xy"] = [[0, 0], [10 / 3, 4], [20 / 3, -4], [10, 0]]
        result = evaluate(state)
        self.assertTrue(result["valid"], result["errors"])
        self.assertTrue(result["shape_diagnostics"]["reverse_bending"])
        self.assertTrue(result["warnings"])

    def test_analytic_length_is_independent_of_display_sampling(self):
        state = synthetic_state()
        # Cubic Bezier representation of x=10t, y=10t(1-t).
        state["arch"]["control_points_xy"] = [[0, 0], [10 / 3, 10 / 3],
                                                [20 / 3, 10 / 3], [10, 0]]
        independent = quad(lambda t: np.sqrt(100 + (10 - 20 * t)**2), 0, 1)[0]
        coarse, fine = evaluate(state, 2), evaluate(state, 1001)
        self.assertTrue(coarse["valid"])
        self.assertAlmostEqual(coarse["A"], independent, places=10)
        self.assertEqual(coarse["A"], fine["A"])
        self.assertGreater(coarse["A"], 10)

    def test_degenerate_cubic_nearest_points(self):
        knots = np.r_[np.zeros(4), np.ones(4)]
        line = BSpline(knots, np.c_[np.arange(4) / 3, np.zeros(4)], 3)
        t, d = footpoints(line, [[0.5, 1], [-2, 2], [2, 2]])
        np.testing.assert_allclose(t, [0.5, 0, 1], atol=1e-12)
        np.testing.assert_allclose(d, [1, np.sqrt(8), np.sqrt(5)], atol=1e-12)
        constant = BSpline(knots, np.zeros((4, 2)), 3)
        _, distance = footpoints(constant, [[1, 1]])
        self.assertAlmostEqual(distance[0], np.sqrt(2))

    def test_explicit_refit_retains_teeth_frame_and_input(self):
        state = synthetic_state()
        state["arch"]["control_points_xy"][1][1] = 2
        original = deepcopy(state)
        updated = refit_state(state)
        self.assertEqual(state, original)
        self.assertEqual(updated["teeth"], original["teeth"])
        self.assertEqual(updated["frame"], original["frame"])
        self.assertTrue(updated["last_refit"]["converged"])
        self.assertAlmostEqual(evaluate(updated)["A"], 10, places=7)
        self.assertLess(evaluate(updated)["A"], evaluate(state)["A"])

    def test_fit_scale_and_translation_equivariance(self):
        t = np.linspace(0.05, 0.95, 20)
        points = np.c_[10 * t, 4 * t * (1 - t)]
        endpoints = np.array([[0, 0], [10, 0]])
        first = fit_arch(points, endpoints)
        scale, shift = 2.3, np.array([-3.2, 7.8])
        second = fit_arch(points * scale + shift, endpoints * scale + shift)
        self.assertTrue(first["converged"])
        self.assertTrue(second["converged"])
        np.testing.assert_allclose(second["control"], first["control"] * scale + shift,
                                   atol=2e-6, rtol=0)
        self.assertAlmostEqual(second["length"], first["length"] * scale, places=6)
        np.testing.assert_allclose(first["control"][[0, -1]], endpoints, atol=1e-12)


if __name__ == "__main__":
    unittest.main()
