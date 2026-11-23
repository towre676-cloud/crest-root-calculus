from fractions import Fraction
import json
from itertools import product
import unittest
from algebra import characteristic, field, parse_form, point_count, simplicity_witness, smoothness
from app import ROOT, cyclic_map, solve


class AtlasTests(unittest.TestCase):
    def test_published_blind_receipt(self):
        request=json.loads((ROOT/'example.json').read_text())
        actual=solve(request)
        self.assertEqual(actual,json.loads((ROOT/'receipt.json').read_text()))
        self.assertEqual(actual['counts'],{'degree_three_elliptic_map_over_Q':1,'no_elliptic_map_over_Qbar':7,'unresolved':5})

    def test_klein_independent_target_invariant(self):
        q=cyclic_map(parse_form('x^3*y+y^3*z+z^3*x'))
        e,d,c,b,a=q['target_Y_squared_coefficients_ascending']
        I=12*a*e-3*b*d+c*c
        J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
        self.assertEqual(Fraction(6912*I**3,4*I**3-J*J),-3375)
        self.assertIsNotNone(q['rational_origin'])

    def test_direct_field_point_oracle(self):
        f=parse_form('x^4+3*x^2*y*z+x*y^3+y*z^3+z^4')
        for p,n in [(2,1),(2,3),(3,2),(5,1)]:
            a,m=field(p,n);q=p**n
            def test(x,y,z):
                accum=0
                for key,c in f.items():
                    t=c%p
                    for base,degree in zip((x,y,z),key):
                        for _ in range(degree): t=m[t][base]
                    accum=a[accum][t]
                return accum==0
            # Use a different affine/projective chart order.
            direct=sum(test(1,y,z) for y,z in product(range(q),repeat=2))
            direct+=sum(test(0,1,z) for z in range(q))+test(0,0,1)
            self.assertEqual(point_count(f,p,n),direct)

    def test_bad_reduction_and_bad_polynomial(self):
        for p in (2,3,5): self.assertIsNone(smoothness(parse_form('x^4'),p))
        self.assertIsNone(simplicity_witness([1,0,0,-2,0,0,1]))
        with self.assertRaises(ValueError): parse_form('x^3+z^4')
        with self.assertRaises(ValueError): parse_form('x^4++z^4')
        with self.assertRaises(ValueError): parse_form('__import__(x)')

if __name__=='__main__': unittest.main()
