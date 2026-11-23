"""Verify arithmetic inputs for the full rational spectrum of the a=3 quartic."""
from itertools import product
import json

import sympy as sp

from map_classes import count_classes


def main():
    r = sp.symbols("r")
    e = r*(r-1)*(r-3)
    g = (r-1)*(r-3)*(r+1)*(r-2)
    x = 2*(r-2)/(r-1)
    if sp.cancel(e.subs(r,x)*(r-1)**4+2*g) != 0:
        raise ArithmeticError("Quadratic-twist identification failed")
    inv = sp.Rational(256*(1-3+3**2)**3,(1-3)**2*3**2)
    if inv != sp.Rational(21952,9):
        raise ArithmeticError("Unexpected j-invariant")

    odd_primes = {3,5,7,11,13,17,19,37,43,67,163}
    excluded, records = set(), []
    for p in (5,7,11,17):
        if sp.discriminant(e,r)%p==0 or sp.discriminant(g,r)%p==0:
            raise ArithmeticError("Bad reduction in the sieve")
        def trace(poly, infinity):
            total = sum(int(sp.legendre_symbol(int(poly.subs(r,i)),p)) for i in range(p))
            return 1-infinity-total
        a,b = trace(e,1),trace(g,2)
        if b != int(sp.legendre_symbol(-2,p))*a or (p==5 and (a,b)!=(2,-2)):
            raise ArithmeticError("Wrong traces on E or its quadratic twist")
        disallowed = sorted(ell for ell in odd_primes if ell != p and
                            pow(a*a-4*p,(ell-1)//2,ell) == ell-1)
        excluded.update(disallowed)
        records.append({"p":p,"trace_E":a,"trace_twist":b,"excluded_primes":disallowed})
    if excluded != odd_primes:
        raise ArithmeticError("The odd prime-isogeny sieve has a gap")
    # Odd valuations of a rational function on P^1, including infinity.
    def odd_places(expr):
        num,den = sp.fraction(sp.cancel(expr))
        out = set()
        deg = []
        for factor in (num,den):
            pol = sp.Poly(factor,r)
            roots = pol.ground_roots()
            if sum(roots.values()) != pol.degree():
                raise ArithmeticError("Unexpected nonrational branch place")
            out ^= {str(val) for val, exponent in roots.items() if exponent%2}
            deg.append(pol.degree())
        if (deg[0]-deg[1])%2:
            out.add("infinity")
        return out
    b1 = odd_places(e)
    b2 = odd_places(e.subs(r,r+1))
    trivial_classes = (set(),b1,b2,b1^b2)
    rows = {}
    # Multiply the actual pullbacks before factoring. The pole at r=1 of
    # the third fraction is retained in the rational function as necessary.
    for bits in product((0,1),repeat=3):
        if not any(bits):
            continue
        aa,bb,cc = bits
        killed = []
        for root in (0,1,3):
            value = (r-root)**aa * (r+1-root)**bb * (x-root)**cc
            if odd_places(value) in trivial_classes:
                killed.append(root)
        expected = [3 if aa else 1] if bb and cc else []
        if killed != expected:
            raise ArithmeticError("Torsion kernel parity rule failed")
        rows["".join(map(str,bits))] = killed
    sample = [count_classes(d) for d in (2,3,4)]
    if [(entry["all_geometric_classes"],entry["geometrically_primitive_classes"],
         entry["primitive_classes_defined_over_Q"]) for entry in sample] != [(5,5,3),(4,4,0),(19,4,2)]:
        raise ArithmeticError("Unexpected geometric map-class counts")
    print(json.dumps({"family_parameter":3,"target":"y^2=x*(x-1)*(x-3)",
                      "third_double_quotient_twist":-2,"j":str(inv),
                      "nonzero_two_torsion_killed_by_direction_mod_2": rows,
                      "geometric_class_counts":sample,
                      "geometric_optimal_degree":"2*(a^2+b^2+c^2)/(2 if b and c odd else 1)",
                      "all_geometric_degrees":"positive even integers and odd multiples of positive d = 3 mod 8",
                      "rational_isogeny_sieve":records,
                      "all_rational_elliptic_degrees":"2*(m^2+n^2), (m,n) != (0,0)",
                      "primitive_rational_elliptic_degrees":"2*(m^2+n^2), gcd(m,n)=1"}, indent=2))


if __name__=="__main__":
    main()
