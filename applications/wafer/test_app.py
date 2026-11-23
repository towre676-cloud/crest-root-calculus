from pathlib import Path
import json
import unittest

from app import library, solve, suggestions

ROOT=Path(__file__).resolve().parent

class WaferTests(unittest.TestCase):
    def fixture(self,name):
        return json.loads((ROOT/name).read_text())
    def test_single_selects_informative_channel(self):
        r=solve(self.fixture('example.json'))
        self.assertEqual((r['within'],r['above'],r['next_reading']['channel']),(44,37,32))
    def test_two_reads_grade_without_identity(self):
        r=solve(self.fixture('two-readings.json'))
        self.assertEqual((r['status'],r['surviving_profiles'],r['bottom_width_range_nm']),('within_limit',3,[126.5,134.0]))
    def test_measured_inconsistent(self):
        r=solve(self.fixture('measured-die.json'))
        self.assertEqual(r['status'],'outside_model')
        self.assertAlmostEqual(r['smallest_maximum_intensity_discrepancy'],0.01401533621)
    def test_touching_intervals_remain_ambiguous(self):
        r=suggestions([0],[1],[[0],[2]],[1],{},1)
        self.assertEqual(r[0]['opposite_decision_pairs_left'],1)
    def test_above_or_prior(self):
        case=self.fixture('example.json')
        case.update(width_limit_nm=111,observations={})
        self.assertEqual(solve(case)['status'],'above_limit')
    def test_validation(self):
        for change in [{'error':-1},{'die':1},{'observations':{'85':0.2}},{'width_limit_nm':'NaN'}]:
            with self.subTest(change=change), self.assertRaises(ValueError):
                solve({**self.fixture('example.json'),**change})

if __name__=='__main__':
    unittest.main()
