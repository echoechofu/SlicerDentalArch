# SPDX-License-Identifier: BSD-3-Clause
"""Case-independent planar arch and native-coordinate space measurements.

This module needs NumPy/SciPy, not Slicer. It does not segment teeth, infer tooth
identities, snap points to surfaces, convert physical units, or modify files.
``evaluate`` never mutates its input. Its signed result is C = R - A: positive
means a net shortage under the chosen reference-arch protocol.

The editing domain deliberately requires an arch that progresses monotonically
from its first boundary toward its second boundary in the transverse direction.
This conservative constraint certifies absence of self-intersections, including
small loops that could be missed by a sampled-polyline intersection test. It
also rejects some non-self-intersecting curves that turn back transversely.
"""

from copy import deepcopy

import numpy as np
from scipy.integrate import quad
from scipy.interpolate import BSpline
from scipy.optimize import minimize


TOOTH_ORDER = [
    "VL_P2", "VL_P1", "VL_C", "VL_I2", "VL_I1",
    "VR_I1", "VR_I2", "VR_C", "VR_P1", "VR_P2",
]
DEFAULT_BOUNDARIES = [
    {"tooth": "VL_M1", "point": "M"},
    {"tooth": "VR_M1", "point": "M"},
]


def _array(value, shape, name):
    result = np.asarray(value, dtype=float)
    if result.shape != shape or not np.isfinite(result).all():
        raise ValueError(f"{name} must contain finite coordinates with shape {shape}")
    return result


def _parts(state):
    """Validate the persisted schema and return independent numeric arrays."""
    if not isinstance(state, dict):
        raise ValueError("State must be a JSON object")
    if state.get("schema_version") != 1:
        raise ValueError("Unsupported or missing state schema_version")
    frame = state["frame"]
    center = _array(frame["center"], (3,), "frame.center")
    axes = _array(frame["axes"], (3, 3), "frame.axes")
    if not np.allclose(axes.T @ axes, np.eye(3), atol=1e-8, rtol=0):
        raise ValueError("Reference frame axes must be orthonormal; scaling is not supported")
    if not np.isclose(np.linalg.det(axes), 1.0, atol=1e-8, rtol=0):
        raise ValueError("Reference frame must be right-handed")
    z = float(frame["plane_z"])
    if not np.isfinite(z):
        raise ValueError("Reference plane height must be finite")

    teeth = {}
    if not isinstance(state["teeth"], dict):
        raise ValueError("teeth must map tooth identifiers to M/D points")
    for tooth_id, item in state["teeth"].items():
        if not isinstance(tooth_id, str) or not tooth_id:
            raise ValueError("Tooth identifiers must be nonempty strings")
        teeth[tooth_id] = {role: _array(item[role], (3,), f"{tooth_id}.{role}")
                           for role in ("M", "D")}
        if np.linalg.norm(teeth[tooth_id]["M"] - teeth[tooth_id]["D"]) <= 1e-12:
            raise ValueError(f"{tooth_id}: M and D must be distinct")
    included = list(state.get("included_teeth", TOOTH_ORDER))
    if not included or len(set(included)) != len(included):
        raise ValueError("included_teeth must be nonempty and contain no duplicates")
    if any(tooth_id not in teeth for tooth_id in included):
        raise ValueError("A tooth required for space analysis is missing")
    boundary_spec = state.get("boundaries", DEFAULT_BOUNDARIES)
    if len(boundary_spec) != 2:
        raise ValueError("Exactly two arch boundaries are required")
    endpoints_native = []
    for spec in boundary_spec:
        if spec["point"] not in ("M", "D") or spec["tooth"] not in teeth:
            raise ValueError("An arch boundary does not reference a defined tooth point")
        endpoints_native.append(teeth[spec["tooth"]][spec["point"]])
    endpoints_xy = ((np.asarray(endpoints_native) - center) @ axes)[:, :2]
    if np.linalg.norm(endpoints_xy[1] - endpoints_xy[0]) <= 1e-10:
        raise ValueError("The projected arch boundaries must be distinct")

    arch = state["arch"]
    controls = np.array(arch["control_points_xy"], dtype=float, copy=True)
    if controls.ndim != 2 or controls.shape[1] != 2 or len(controls) < 4:
        raise ValueError("The arch requires at least four planar control points")
    if not np.isfinite(controls).all():
        raise ValueError("Arch control points must be finite")
    degree = arch["degree"]
    if degree != 3:
        raise ValueError("Only cubic B-spline reference arches are supported")
    knots = _array(arch["knots"], (len(controls) + 4,), "arch.knots")
    if np.any(np.diff(knots) < 0) or knots[3] >= knots[-4]:
        raise ValueError("Arch knots must be ordered over a nonzero domain")
    if not (np.all(knots[:4] == knots[3]) and np.all(knots[-4:] == knots[-4])):
        raise ValueError("Arch must have four repeated knots at each endpoint")
    interior = knots[4:-4]
    if len(interior) and not (np.all(np.diff(interior) > 0)
                             and interior[0] > knots[3] and interior[-1] < knots[-4]):
        raise ValueError("Interior knots must be simple and strictly inside the domain")
    controls[[0, -1]] = endpoints_xy
    return center, axes, z, teeth, included, controls, knots


def synchronize_endpoints(state):
    """Update only the two bound control points, retaining all interior controls.

    Unlike ``evaluate``, this helper intentionally updates and returns ``state``.
    It never recalculates the reference frame or changes any tooth point.
    """
    *_, controls, knots = _parts(state)
    state["arch"]["control_points_xy"] = controls.tolist()
    return state


def _real_roots_in_unit_interval(coefficients):
    coefficients = np.asarray(coefficients, dtype=float)
    size = float(np.max(np.abs(coefficients))) if len(coefficients) else 0.0
    if size == 0:
        return []
    while len(coefficients) > 1 and abs(coefficients[-1]) <= 1e-12 * size:
        coefficients = coefficients[:-1]
    return [float(root.real) for root in np.polynomial.polynomial.polyroots(coefficients)
            if abs(root.imag) < 1e-8 and 0 < root.real < 1]


def _polynomial_range(coefficients):
    points = [0.0, 1.0] + _real_roots_in_unit_interval(
        np.polynomial.polynomial.polyder(coefficients))
    values = np.polynomial.polynomial.polyval(points, coefficients)
    return float(np.min(values)), float(np.max(values))


def _span_coefficients(curve, a, b):
    """Polynomial coefficients in s=(t-a)/(b-a), ascending powers."""
    span = b - a
    return np.array([curve(a), curve(a, nu=1) * span,
                     curve(a, nu=2) * span**2 / 2,
                     curve(a, nu=3) * span**3 / 6])


def shape_diagnostics(curve):
    """Check polynomial-span extrema, not only the display sampling grid."""
    breaks = np.unique(curve.t[3:-3])
    endpoints = curve([breaks[0], breaks[-1]])
    separation = float(np.linalg.norm(endpoints[1] - endpoints[0]))
    transverse = (endpoints[1] - endpoints[0]) / separation
    minimum_speed_squared = np.inf
    minimum_transverse_speed = np.inf
    minimum_progress = np.inf
    turn_minimum, turn_maximum = np.inf, -np.inf
    for a, b in zip(breaks[:-1], breaks[1:]):
        co = _span_coefficients(curve, a, b)
        velocity = np.array([co[1], 2 * co[2], 3 * co[3]]) / (b - a)
        acceleration = np.array([velocity[1], 2 * velocity[2]]) / (b - a)
        square = np.polynomial.polynomial.polyadd(
            np.polynomial.polynomial.polymul(velocity[:, 0], velocity[:, 0]),
            np.polynomial.polynomial.polymul(velocity[:, 1], velocity[:, 1]))
        minimum_speed_squared = min(minimum_speed_squared, _polynomial_range(square)[0])
        minimum_transverse_speed = min(minimum_transverse_speed,
                                       _polynomial_range(velocity @ transverse)[0])
        minimum_progress = min(minimum_progress,
                               float((curve(b) - curve(a)) @ transverse))
        cross = np.polynomial.polynomial.polysub(
            np.polynomial.polynomial.polymul(velocity[:, 0], acceleration[:, 1]),
            np.polynomial.polynomial.polymul(velocity[:, 1], acceleration[:, 0]))
        low, high = _polynomial_range(cross)
        turn_minimum, turn_maximum = min(turn_minimum, low), max(turn_maximum, high)
    duration = breaks[-1] - breaks[0]
    speed_scale = separation / duration
    tolerance = speed_scale * 1e-10
    monotone = (minimum_transverse_speed >= -tolerance
                and minimum_progress > separation * 1e-10)
    regular = minimum_speed_squared > (speed_scale * 1e-9)**2
    turn_tolerance = separation**2 / duration**3 * 1e-9
    return {
        "min_speed": float(np.sqrt(max(minimum_speed_squared, 0))),
        "min_transverse_derivative": minimum_transverse_speed,
        "min_span_transverse_progress": minimum_progress,
        "min_turn_cross_product": turn_minimum,
        "max_turn_cross_product": turn_maximum,
        "transverse_monotonicity_certified": bool(monotone),
        "nonzero_speed_certified": bool(regular),
        "no_self_intersection_certified": bool(monotone and regular),
        "reverse_bending": bool(turn_minimum < -turn_tolerance
                                 and turn_maximum > turn_tolerance),
        "validation_domain": "Transversely monotone planar cubic B-spline",
    }


def arc_length(curve):
    """Integrate the analytic curve, independent of display polyline density."""
    breaks = np.unique(curve.t[curve.k:len(curve.t) - curve.k])
    value, error = quad(lambda t: float(np.linalg.norm(curve(t, nu=1))),
                        float(breaks[0]), float(breaks[-1]),
                        points=breaks[1:-1], epsabs=1e-9, epsrel=1e-10, limit=100)
    return float(value), float(error)


def evaluate(state, samples=501):
    """Return current measurements, or ``valid=False`` with no stale totals.

    All coordinates and lengths use the input's native units. The boundary
    control points are projected into the saved plane, leaving five interior
    controls (for the default seven-control arch) untouched. Invalid curve
    geometry may still return a finite preview to let a user correct it.
    """
    result = {"valid": False, "errors": [], "warnings": [], "R": None,
              "A": None, "C": None, "D": None, "widths": {},
              "curve_native": [], "controls_native": [], "controls_xy": []}
    try:
        if not isinstance(samples, int) or samples < 2 or samples > 100001:
            raise ValueError("samples must be an integer between 2 and 100001")
        center, axes, z, teeth, included, controls, knots = _parts(state)
        curve = BSpline(knots, controls, 3)
        t = np.linspace(knots[3], knots[-4], samples)
        xy = curve(t)
        result["curve_native"] = (np.column_stack([xy, np.full(len(xy), z)])
                                   @ axes.T + center).tolist()
        result["controls_native"] = (np.column_stack([controls, np.full(len(controls), z)])
                                      @ axes.T + center).tolist()
        result["controls_xy"] = controls.tolist()
        result["shape_diagnostics"] = shape_diagnostics(curve)
        shape = result["shape_diagnostics"]
        if not shape["nonzero_speed_certified"]:
            result["errors"].append("The curve has a stationary or near-stationary point; adjust the handles")
        if not shape["transverse_monotonicity_certified"]:
            result["errors"].append("The arch folds back transversely; this shape is outside the supported editing domain")
        if shape["reverse_bending"]:
            result["warnings"].append("The arch changes bending direction; review the local shape")
        if result["errors"]:
            return result
        widths = {key: float(np.linalg.norm(p["M"] - p["D"])) for key, p in teeth.items()}
        required = float(sum(widths[key] for key in included))
        available, integration_error = arc_length(curve)
        if not np.isfinite([required, available]).all():
            raise ValueError("Non-finite measurement result")
        result.update(valid=True, widths=widths, R=required, A=available,
                      C=required - available, D=available - required,
                      integration_error_estimate=integration_error,
                      units=state.get("units", "native STL coordinate units"))
    except (ValueError, TypeError, KeyError, IndexError, np.linalg.LinAlgError) as error:
        result["errors"].append(str(error))
    return result


def basis(n=7):
    if not isinstance(n, (int, np.integer)) or n < 4:
        raise ValueError("At least four cubic controls are required")
    knots = np.r_[np.zeros(4), np.linspace(0, 1, n - 2)[1:-1], np.ones(4)]
    return BSpline(knots, np.eye(n), 3), knots


def bending_matrix(bs):
    matrix = np.zeros((bs.c.shape[0], bs.c.shape[0]))
    nodes, weights = np.polynomial.legendre.leggauss(3)
    breaks = np.unique(bs.t)
    for a, b in zip(breaks[:-1], breaks[1:]):
        t = (nodes + 1) * (b - a) / 2 + a
        derivative = bs(t, nu=2)
        matrix += derivative.T @ (derivative * (weights * (b - a) / 2)[:, None])
    return matrix


def footpoints(curve, points):
    """Nearest points over all cubic spans, including degenerate lower degrees."""
    points = np.asarray(points, dtype=float)
    if points.ndim != 2 or points.shape[1] != 2 or not np.isfinite(points).all():
        raise ValueError("Footpoint observations must be finite two-dimensional points")
    breaks = np.unique(curve.t[curve.k:len(curve.t) - curve.k])
    parameters, distances = [], []
    for point in points:
        candidates = list(breaks)
        for a, b in zip(breaks[:-1], breaks[1:]):
            coeff = _span_coefficients(curve, a, b)
            polynomial = np.zeros(6)
            for axis in range(2):
                offset = coeff[:, axis].copy()
                offset[0] -= point[axis]
                product = np.polynomial.polynomial.polymul(
                    offset, [coeff[1, axis], 2 * coeff[2, axis], 3 * coeff[3, axis]])
                polynomial[:len(product)] += product
            candidates.extend(a + root * (b - a)
                              for root in _real_roots_in_unit_interval(polynomial))
        candidates = np.asarray(candidates)
        square = np.sum((curve(candidates) - point)**2, axis=1)
        best = int(np.argmin(square))
        parameters.append(float(candidates[best]))
        distances.append(float(np.sqrt(square[best])))
    return np.array(parameters), np.array(distances)


def huber_values(r, delta):
    return np.where(r <= delta, 0.5 * r * r, delta * (r - 0.5 * delta))


def fit_arch(points, endpoints, weights=None, n_control=7, smooth=1e-4,
             huber_fraction=0.015, max_iter=300):
    """Endpoint-fixed robust fit; local convergence, not a global-optimum claim.

    Observation order follows the arch from its first boundary to its second.
    The smoothness penalty is parameter second derivative, not curvature squared.
    The former application-specific fitting mathematics is retained here without
    its source-file identity or tooth-segmentation assumptions.
    """
    points = np.asarray(points, dtype=float)
    ends = _array(endpoints, (2, 2), "endpoints")
    if points.ndim != 2 or points.shape[1] != 2 or len(points) < 2 or not np.isfinite(points).all():
        raise ValueError("At least two finite 2D observations are required")
    weights = np.full(len(points), 0.5) if weights is None else _array(weights, (len(points),), "weights")
    if np.any(weights <= 0) or not np.isfinite(smooth) or smooth <= 0:
        raise ValueError("Weights and smoothing must be positive")
    if not np.isfinite(huber_fraction) or huber_fraction <= 0:
        raise ValueError("Huber fraction must be positive")
    origin = ends.mean(0)
    scale = float(np.linalg.norm(ends[1] - ends[0]))
    if scale <= 1e-10:
        raise ValueError("Projected boundaries must be distinct")
    p, e = (points - origin) / scale, (ends - origin) / scale
    bs, knots = basis(n_control)
    penalty = bending_matrix(bs)
    sequence = np.vstack([e[0], p, e[1]])
    steps = np.linalg.norm(np.diff(sequence, axis=0), axis=1)
    t = np.cumsum(steps)[:-1] / sum(steps)
    robust = np.ones(len(points))
    control = np.zeros((n_control, 2))
    control[[0, -1]] = e
    free, fixed = np.arange(1, n_control - 1), np.array([0, n_control - 1])
    last = None
    for iteration in range(25):
        design, w = bs(t), weights * robust
        normal = design.T @ (w[:, None] * design) + smooth * penalty
        rhs = design.T @ (w[:, None] * p) - normal[:, fixed] @ e
        control[free] = np.linalg.solve(normal[np.ix_(free, free)], rhs[free])
        curve = BSpline(knots, control, 3)
        t, distance = footpoints(curve, p)
        robust = np.minimum(1, huber_fraction / np.maximum(distance, 1e-14))
        sample = curve(np.linspace(0, 1, 101))
        change = float(np.max(np.linalg.norm(sample - last, axis=1))) if last is not None else np.inf
        last = sample
        if change < 1e-9 and iteration > 3:
            break

    def objective_gradient(flat):
        ctrl = control.copy()
        ctrl[free] = flat.reshape(-1, 2)
        candidate = BSpline(knots, ctrl, 3)
        parameters, distances = footpoints(candidate, p)
        design = bs(parameters)
        residual = candidate(parameters) - p
        w = np.minimum(1, huber_fraction / np.maximum(distances, 1e-14)) * weights
        value = float(weights @ huber_values(distances, huber_fraction)
                      + 0.5 * smooth * np.sum(ctrl * (penalty @ ctrl)))
        gradient = design.T @ (w[:, None] * residual) + smooth * (penalty @ ctrl)
        return value, gradient[free].ravel()

    optimized = minimize(objective_gradient, control[free].ravel(), jac=True,
                         method="L-BFGS-B", options={"maxiter": max_iter, "ftol": 1e-14,
                                                    "gtol": 1e-9, "maxls": 40, "maxcor": 15})
    control[free] = optimized.x.reshape(-1, 2)
    objective, gradient = objective_gradient(optimized.x)
    converged = bool(optimized.success and np.linalg.norm(gradient, ord=np.inf) < 1e-7)
    native_control = control * scale + origin
    curve = BSpline(knots, native_control, 3)
    t, distances = footpoints(curve, points)
    length, error = arc_length(curve)
    return {"curve": curve, "control": native_control, "knots": knots, "t": t,
            "residuals": distances, "converged": converged, "length": length,
            "integration_error": error, "iterations": int(optimized.nit),
            "optimizer_message": str(optimized.message), "objective": objective,
            "gradient_infinity_norm": float(np.linalg.norm(gradient, ord=np.inf)),
            "shape_diagnostics": shape_diagnostics(curve)}


def refit_state(state):
    """Explicitly replace the manual arch in a deep copy of a valid input state.

    Tooth IDs are not parsed for side. ``included_teeth`` is ordered from one
    posterior boundary through the midline to the opposite boundary: D/M pairs
    on its first half, M/D pairs on its second half. An optional
    ``fit.point_order`` list of ``[tooth_id, role]`` pairs overrides this order.
    Failed optimization or invalid geometry raises ValueError; the input remains
    unchanged so the caller can retain the user's last manual curve.
    """
    center, axes, _, teeth, included, controls, _ = _parts(state)
    settings = state.get("fit", {})
    default_order = [(name, role) for i, name in enumerate(included)
                     for role in (("D", "M") if i < len(included) / 2 else ("M", "D"))]
    order = settings.get("point_order", default_order)
    expected = {(name, role) for name in included for role in ("M", "D")}
    if len(order) != len(expected) or {tuple(pair) for pair in order} != expected:
        raise ValueError("fit.point_order must list both points of every included tooth exactly once")
    points = np.array([(teeth[name][role] - center) @ axes for name, role in order])[:, :2]
    fitted = fit_arch(points, controls[[0, -1]], n_control=len(controls),
                      smooth=float(settings.get("smooth", 1e-4)),
                      huber_fraction=float(settings.get("huber_fraction", 0.015)))
    if not fitted["converged"]:
        raise ValueError("Reference arch refit did not converge: " + fitted["optimizer_message"])
    updated = deepcopy(state)
    updated["arch"] = {"control_points_xy": fitted["control"].tolist(),
                       "knots": fitted["knots"].tolist(), "degree": 3}
    checked = evaluate(updated)
    if not checked["valid"]:
        raise ValueError("Refitted arch is outside the valid editing domain: " + "; ".join(checked["errors"]))
    updated["last_refit"] = {"converged": True, "iterations": fitted["iterations"],
                             "rms_residual": float(np.sqrt(np.mean(fitted["residuals"]**2))),
                             "gradient_infinity_norm": fitted["gradient_infinity_norm"]}
    return updated
