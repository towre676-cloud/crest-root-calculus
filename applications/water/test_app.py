import copy
import json
from pathlib import Path
import unittest
from app import solve
from field_replay import replay

ROOT = Path(__file__).parent


class Search(unittest.TestCase):
    def setUp(self):
        self.d = json.loads((ROOT / "example.json").read_text())

    def test_next_visit(self):
        r = solve(self.d)
        self.assertEqual(r["zones"], ["Mill", "Orchard"])
        self.assertEqual(r["next_logger"], "North")
        self.assertEqual(r["logger_comparison"][0]["guaranteed_reduction"], 120)
        self.d["observations"]["North"] = [2.95, 3.05]
        self.assertEqual(solve(self.d)["zones"], ["Orchard"])

    def test_joint_strength_not_independent_per_reading(self):
        self.d["zones"] = {"Orchard": 180}
        first = self.d["scenarios"][0]
        first["slope"]["North"] = -1
        first["offset"] = {}
        self.d["scenarios"] = [first]
        self.d["observations"]["North"] = [-3, -2.5]
        self.assertEqual(solve(self.d)["status"], "model_mismatch")
        self.d["observations"]["North"] = [-2.1, -2]
        self.assertEqual(
            solve(self.d)["surviving_scenarios"][0]["strength"], [2.0, 2.1]
        )

    def test_shared_endpoint_and_multiple_models_per_street(self):
        d = {
            "extent_unit": "sites",
            "sensors": {"h": {"cost": 1, "error": 0}},
            "zones": {"a": 1, "b": 1},
            "scenarios": [
                {"zone": "a", "slope": {"h": 1}, "strength": [0, 1]},
                {"zone": "b", "slope": {"h": 1}, "strength": [1, 2]},
            ],
        }
        d["scenarios"].append(copy.deepcopy(d["scenarios"][0]))
        self.assertIsNone(solve(d)["next_logger"])
        d["observations"] = {"h": [1, 1]}
        self.assertEqual(solve(d)["search_extent"], 2)
        d["observations"]["h"] = [3, 4]
        self.assertEqual(solve(d)["status"], "model_mismatch")
        with self.assertRaises(ValueError):
            d["zones"]["a"] = -2
            solve(d)

    def test_field(self):
        d = json.loads((ROOT / "field.json").read_text())
        r = replay(d)
        self.assertEqual(r["events_scored"], 8)
        self.assertEqual(len(r["excluded_events"]), 1)
        self.assertEqual(r["summary"]["selected"]["true_site_retained"], 3)
        self.assertEqual(
            r["later_day"]["after_calibration"]["selected"]["true_site_retained"], 4
        )
        self.assertEqual(
            r["later_day"]["after_calibration"]["fixed"]["total_sites_to_check"], 28
        )


if __name__ == "__main__":
    unittest.main()
