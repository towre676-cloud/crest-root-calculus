from fractions import Fraction as F
import json
import unittest
from app import ROOT, between, cycles, evaluate, fit, prediction
from benchmark import bench

def candidate(epoch=(0,0),period=(1,1),dur=(2.4,2.4)):
    return {'name':'one','epoch_days':list(epoch),'period_days':list(period),'duration_hours':list(dur)}

class TimingTests(unittest.TestCase):
    def test_exact_full_and_partial(self):
        polys,duration=fit(candidate(),[])
        self.assertEqual(prediction(polys,duration,between([-.05,.05])),{1})
        self.assertEqual(prediction(polys,duration,between([-.01,.05])),{0})
        self.assertEqual(prediction(polys,duration,between([.1,.3])),{0})

    def test_ambiguous_duration_and_interval(self):
        polys,d=fit(candidate(epoch=(-.05,.05),dur=(1.2,2.4)),[])
        self.assertEqual(prediction(polys,d,between([-.03,.03])),{0,1})
        self.assertEqual(prediction(polys,d,between([-2,2])),{1})

    def test_joint_constraint(self):
        h=candidate(period=(1,1.2))
        self.assertTrue(fit(h,[between([1.1,1.1])])[0])
        self.assertTrue(fit(h,[between([2.4,2.4])])[0])
        self.assertFalse(fit(h,[between([1.1,1.1]),between([2.4,2.4])])[0])
        self.assertTrue(fit(h,[between([-2.2,-2.2]),between([1.1,1.1])])[0])

    def test_midpoint_and_non_detection(self):
        q=json.loads((ROOT/'example.json').read_text())
        rec=evaluate(q)
        choice=next(w for w in q['windows'] if w['name']==rec['recommended_window'])
        r=evaluate(dict(q,non_detections=[choice['days']]))
        self.assertEqual(r['surviving_hypotheses'],['daily alias'])
        self.assertEqual(evaluate(dict(q,windows=[]))['decision'],'no_guaranteed_separation')
        self.assertEqual(evaluate(dict(q,events=[[2455568.25,2455568.26]]))['decision'],'model_mismatch')

    def test_replay(self):
        r=bench()
        self.assertEqual(r,json.loads((ROOT/'benchmark.json').read_text()))
        self.assertEqual(r['measured_RV_timestamps'],320)
        self.assertEqual(evaluate(json.loads((ROOT/'confirmed.json').read_text())),json.loads((ROOT/'confirmation.json').read_text()))
    def test_input_validation(self):
        q=json.loads((ROOT/'example.json').read_text())
        for change in [{'hypotheses':[]},{'time_system':''},{'windows':[{'name':'x','available':True,'days':[2,1]}]}]:
            with self.assertRaises(ValueError):evaluate({**q,**change})
if __name__=='__main__':unittest.main()
