"""Regression fixtures and independent floating-point pose checks for LW-159-R1."""

import copy
from fractions import Fraction as F
from itertools import product
import json
import math
from pathlib import Path
import unittest

from bounds import Box, S, circle_arcs
from check import exact_json, min_angle_bounds, verify

HERE = Path(__file__).resolve().parent


def reference_pose(g, a, b, c, t):
    """Start at D, rotate D->A clockwise by the law-of-cosines angle at D."""
    ax, ay = a * math.cos(t), a * math.sin(t)
    ray = math.atan2(ay, ax - g)
    d = math.hypot(ax - g, ay)
    phi = ray - math.acos((c * c + d * d - b * b) / (2 * c * d))
    ux, uy = math.cos(phi), math.sin(phi)
    return (ax, ay), (g + c * ux, c * uy), (g + (c + 18.3) * ux, (c + 18.3) * uy)


def point_to_segment(p, a, b):
    dx, dy = b[0] - a[0], b[1] - a[1]
    w = max(0, min(1, ((p[0] - a[0]) * dx + (p[1] - a[1]) * dy) / (dx * dx + dy * dy)))
    return math.hypot(p[0] - a[0] - w * dx, p[1] - a[1] - w * dy)


def distance(a, b, c, d):
    def cross(p, q, r):
        return (q[0] - p[0]) * (r[1] - p[1]) - (q[1] - p[1]) * (r[0] - p[0])

    if cross(a, b, c) * cross(a, b, d) < 0 and cross(c, d, a) * cross(c, d, b) < 0:
        return 0
    return min(
        point_to_segment(a, c, d),
        point_to_segment(b, c, d),
        point_to_segment(c, a, b),
        point_to_segment(d, a, b),
    )


class TestReference(unittest.TestCase):
    @classmethod
    def setUpClass(cls):
        cls.spec = exact_json(HERE.parent / "kit" / "design.json")
        cls.data = verify(cls.spec)

    def test_recomputed_receipt(self):
        self.assertEqual(self.data, json.loads((HERE / "receipt.json").read_text()))

    def test_rational_interval_primitives(self):
        for x1, x2, y1, y2 in [
            ("-1.333", "3.1", "1.12", "3.75"),
            ("-1", "-0.001", "11314", "22003"),
        ]:
            a, b = Box.enclose(x1, x2), Box.enclose(y1, y2)
            for u, v in product(map(F, (x1, x2)), map(F, (y1, y2))):
                for box, val in [
                    (a + b, u + v),
                    (a - b, u - v),
                    (a * b, u * v),
                    (a / b, u / v),
                ]:
                    self.assertLessEqual(box.lo, val * S)
                    self.assertGreaterEqual(box.hi, val * S)
        for a, b in circle_arcs(32):
            radius = a.square() + b.square()
            self.assertLessEqual(radius.lo, S)
            self.assertGreaterEqual(radius.hi, S)

    def test_fixed_length_realized_angles(self):
        rows = self.data["motion"]
        self.assertAlmostEqual(rows["nominal_swing_deg"], 59.9847637275, places=5)
        for idx, radius in enumerate((46.6 + 84.1, 84.1 - 46.6)):
            t = math.radians(rows["nominal_crank_reversal_deg_ccw"][idx])
            a, b, _ = reference_pose(93.8, 46.6, 84.1, 111.3, t)
            self.assertAlmostEqual(math.hypot(*b), radius, places=6)
            self.assertLess(
                abs(a[0] * b[1] - a[1] * b[0]) / (math.hypot(*a) * math.hypot(*b)), 3e-7
            )
            phi = math.degrees(math.atan2(b[1], b[0] - 93.8))
            self.assertAlmostEqual(
                phi, rows["nominal_rocker_extrema_deg"][idx], places=5
            )

    def test_shaft_length_corners_and_interior_poses(self):
        geometry = self.data["geometry"]
        dims = [
            [float(F(x)) for x in row]
            for row in geometry["effective_pin_axis_lengths_mm"]
        ]
        mm = self.data["geometry"]["nonincident_pair_clearance_lower_mm"][
            "crank / rocker"
        ]
        transmission = geometry["closure_and_transmission"][
            "acute_transmission_lower_deg"
        ]
        for g, a, b, c in product(*dims):
            for i in range(180):
                t = math.tau * i / 180
                A, B, T = reference_pose(g, a, b, c, t)
                O, D = (0, 0), (g, 0)
                tri = (
                    (A[0] - B[0]) * (D[0] - B[0]) + (A[1] - B[1]) * (D[1] - B[1])
                ) / (b * c)
                self.assertGreaterEqual(
                    math.degrees(math.acos(abs(tri))), transmission - 1e-7
                )
                self.assertGreaterEqual(distance(O, A, D, T) - 2 * 7.51, mm)
                self.assertAlmostEqual(math.dist(A, B), b, places=9)
                self.assertGreater(
                    (D[0] - A[0]) * (B[1] - A[1]) - (D[1] - A[1]) * (B[0] - A[0]), 0
                )

    def test_rounded_five_mm_toggle(self):
        dims = [F(s) for s in ("93.8", "46.6", "84.1", "111.3")]
        with self.assertRaisesRegex(ValueError, "strict closure"):
            min_angle_bounds([(v - 5, v + 5) for v in dims])

    def test_bad_fabrication_is_rejected(self):
        cases = [
            ("joint", "moving_bore", [F("2.96"), F("3.1")], "free clearance"),
            ("joint", "grip_blind_bore_diameter", [F("2.96"), F("3.1")], "slide"),
            ("joint", "washer_od_min", F("3.5"), "escape"),
            ("joint", "ground_spacer", F("1"), "moving_shafts_over_base"),
            ("base", "foot_height", F("5"), "ground_retention_over_table"),
            ("base", "ymax", F("124"), "top"),
            ("links", "width", F("29"), "crank.*rocker"),
            ("datum", "branch", "lower", "upper assembly"),
        ]
        for group, key, val, message in cases:
            with self.subTest(key=key):
                bad = copy.deepcopy(self.spec)
                bad[group][key] = val
                with self.assertRaisesRegex(ValueError, message):
                    verify(bad)


if __name__ == "__main__":
    unittest.main()
