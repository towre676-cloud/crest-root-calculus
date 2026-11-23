import json
import unittest
from app import ROOT, STAGES, evaluate, library
from benchmark import benchmark

def trace(ident, capacity, volts):
    return {'id':ident,'capacity_ah':capacity,'volts':{str(t):volts for t in STAGES}}

class GradeTests(unittest.TestCase):
    def setUp(self):
        self.q={'min_capacity_ah':4.8,'voltage_tolerance_v':0.01,'capacity_error_ah':0,
                'reference':[trace('low',4.7,3.31),trace('high',4.8,3.33)],'observations':{}}

    def test_closed_endpoints_and_grade_boundary(self):
        q=dict(self.q, observations={'60':3.32})
        self.assertEqual(evaluate(q)['candidate_count'],2)
        self.assertEqual(evaluate(dict(q,observations={'60':3.34}))['decision'],'meets_capacity_grade')
        self.assertEqual(evaluate(dict(q,observations={'60':3.30}))['decision'],'below_capacity_grade')

    def test_mismatch_and_exhausted_ambiguity(self):
        self.assertEqual(evaluate(dict(self.q,observations={'60':4.0}))['decision'],'reference_test')
        r=evaluate(dict(self.q,observations={'1800':3.32}))
        self.assertEqual((r['decision'],r['candidate_count']),('reference_test',2))

    def test_single_uncertain_capacity_needs_reference(self):
        r=evaluate(dict(self.q,capacity_error_ah=.1,observations={'60':3.33}))
        self.assertEqual(r['decision'],'reference_test')

    def test_budget_and_reference_examples(self):
        q=json.loads((ROOT/'grade.json').read_text())
        self.assertEqual(evaluate(q)['decision'],'meets_capacity_grade')
        self.assertEqual(evaluate(json.loads((ROOT/'deferred.json').read_text()))['decision'],'reference_test')
        self.assertEqual(evaluate(dict(q,observations={'300':3.4}))['decision'],'reference_test')
        with self.assertRaises(ValueError): evaluate(dict(q,observations={'600':3.4}))

    def test_bad_input(self):
        for update in [{'voltage_tolerance_v':-1},{'observations':{'500':3}}, {'reference':[]},
                       {'min_capacity_ah':'NaN'}, {'reference':[trace('same',4,3),trace('same',5,4)]}]:
            with self.assertRaises(ValueError): evaluate({**self.q,**update})

    def test_demo_and_heldout_receipts(self):
        self.assertEqual(len(library()),23)
        self.assertEqual(evaluate(json.loads((ROOT/'example.json').read_text())),json.loads((ROOT/'results.json').read_text()))
        self.assertEqual(benchmark(),json.loads((ROOT/'benchmark.json').read_text()))
if __name__=='__main__':unittest.main()
