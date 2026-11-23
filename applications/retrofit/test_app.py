import json
from pathlib import Path
import unittest
from app import solve
from field_replay import replay


class Decisions(unittest.TestCase):
    def setUp(self):
        self.d = json.loads(Path(__file__).with_name("example.json").read_text())

    def test_coheat_cannot_choose(self):
        r = solve(self.d)
        self.assertEqual(r["status"], "needs_measurement")
        self.assertEqual(r["next_test"], "Effective infiltration survey")
        self.assertAlmostEqual(
            r["test_comparison"][0]["worst_remaining_advantage_width"], 0.1152
        )
        self.assertEqual(r["test_comparison"][-1]["guaranteed_width_reduction"], 0)

    def test_same_total_opposite_decisions(self):
        for air, winner in (
            ([32, 48], "Insulate fabric"),
            ([170, 185], "Seal unintended leaks"),
        ):
            self.d["measured"] = [{"quantity": "infiltration", "bounds_W_K": air}]
            self.assertEqual(solve(self.d)["chosen_upgrade"], winner)

    def test_tie_point_and_inconsistent(self):
        self.d["actions"][1] = dict(self.d["actions"][0], id="Copy")
        self.d["total_W_K"] = [100, 100]
        self.d["measured"] = [{"quantity": "fabric", "bounds_W_K": [100, 100]}]
        self.assertEqual(solve(self.d)["status"], "equal_value")
        self.d["measured"][0]["bounds_W_K"] = [110, 120]
        self.assertEqual(solve(self.d)["status"], "model_mismatch")

    def test_dominance_only_touches_tie_boundary(self):
        self.d["actions"][1]["fraction_infiltration_saved"] = 0
        for _ in range(2):
            r = solve(self.d)
            self.assertEqual(r["status"], "one_upgrade_never_worse")
            self.assertEqual(r["chosen_upgrade"], "Insulate fabric")
            self.assertIsNone(r["next_test"])
            self.assertIn(0, r["first_minus_second_kWh_per_year_per_cost"])
            self.d["actions"].reverse()

    def test_unknown_axis(self):
        self.d["tests"][0]["quantity"] = "raw blower flow"
        with self.assertRaises(KeyError):
            solve(self.d)

    def test_published_field_cases(self):
        a, b = replay()
        self.assertEqual(a["independent_endpoint_change_W_K"], [-4, 88])
        self.assertEqual(b["independent_endpoint_change_W_K"], [-10, 30])
        self.assertEqual([a["measured_change_W_K"], b["measured_change_W_K"]], [42, 10])


if __name__ == "__main__":
    unittest.main()
