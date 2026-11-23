import json
from collections import defaultdict
from itertools import combinations
import unittest
from app import ROOT, canonical, distances, local_count, solve


class MaterialTests(unittest.TestCase):
    def test_bloom_and_local_measurement(self):
        case=json.loads((ROOT/'example.json').read_text())
        r=solve(case)
        self.assertEqual(len(r['candidates']),2)
        self.assertEqual(r['next_probe']['predicted_count_per_candidate'],[1,0])
        for count,p in zip([1,0],r['candidates']):
            a=solve({**case,'observations':[{'sites':3,'span':4,'count':count}]})
            self.assertEqual(a['candidates'],[p])
        self.assertIsNone(solve({**case,'probes':[]})['next_probe'])

    def test_all_small_exact_fibres_against_enumeration(self):
        for width in range(3,11):
            table=defaultdict(set)
            for pts in combinations(range(1,width),2):
                full=(0,*pts,width)
                table[tuple(distances(full))].add(canonical(full))
            for obs,actual in table.items():
                out=solve({'distances':list(obs)})
                self.assertTrue(out['complete'])
                self.assertEqual(set(map(tuple,out['candidates'])),actual)

    def test_uncertain_multiset_against_bruteforce(self):
        for obs in [[1,2,4],[0,2,4],[2,4,5],[1,2,3,5,6,8]]:
            n=3 if len(obs)==3 else 4
            oracle=set()
            for width in range(max(1,max(obs)-1),max(obs)+2):
                for interior in combinations(range(1,width),n-2):
                    z=(0,*interior,width)
                    if all(abs(a-b)<=1 for a,b in zip(distances(z),obs)):
                        oracle.add(canonical(z))
            r=solve({'distances':obs,'distance_error':1})
            self.assertEqual(set(map(tuple,r['candidates'])),oracle)

    def test_measured_powder_and_reference(self):
        data=json.loads((ROOT/'powder.json').read_text())
        measured=solve(data)
        self.assertEqual(measured['status'],'ambiguous')
        self.assertEqual({c:len(v) for c,v in measured['classes'].items()},{'P':3,'I':2,'F':1})
        filtered=solve(json.loads((ROOT/'powder_reference.json').read_text()))
        self.assertEqual(list(filtered['classes']),['P'])
        self.assertEqual(filtered['classes']['P'][0]['a_A_display_approx'],[4.1569065,4.1569259])
        data['cell_range_A']=[3,3.01]
        self.assertEqual(solve(data)['status'],'incompatible')
        with self.assertRaises(ValueError): solve({**data,'d_spacing_intervals_A':[[1,2],[2,3]]})

    def test_cubic_boundary_is_closed(self):
        data={'mode':'cubic_powder','cell_range_A':[1,1],'d_spacing_intervals_A':[[1,1]]}
        out=solve(data)
        self.assertEqual(out['classes']['P'][0]['a_squared_A2_exact'],['1','1'])
        self.assertNotIn('I',out['classes'])
        self.assertNotIn('F',out['classes'])

    def test_budgets_and_contradictions(self):
        current=solve({'distances':[1,2,3],'node_budget':1})
        self.assertEqual(current['status'],'search_incomplete')
        self.assertIsNone(current['next_probe'])
        self.assertEqual(solve({'distances':[1,1,8]})['status'],'incompatible')
        with self.assertRaises(ValueError): solve({'distances':[1,2]})

    def test_local_counts_by_combinations(self):
        sites=(0,1,2,4,7,9)
        for k in range(2,6):
            for w in range(10):
                self.assertEqual(local_count(sites,k,w),sum(v[-1]-v[0]<=w for v in combinations(sites,k)))

if __name__=='__main__': unittest.main()
