"""Exact K-points on two families of covers of C6 from a certified C6(K) list.

The baseline fixed_curve module supplies C6 points. A caller with an equivalent
field-matched result can pass it explicitly as the base argument. The completeness
status then refers to that result. All affine points are checked independently.
"""
import argparse
import json
from fractions import Fraction
from math import gcd
from cypari2 import Pari

pari = Pari()
w = pari("w")
t = pari("x")


def positive_integer(n):
    if type(n) is not int or n < 1:
        raise ValueError("The exponents must be positive integers.")


def lift_fibres(field, base, horizontal, vertical, infinity_order):
    if str(field) != base["field"] or base["status"] not in {"complete", "partial"}:
        raise ValueError("The input point list needs a matching field and completeness status.")
    nf = pari.nfinit(field)
    one = pari.Mod(1, field)
    affine = set()
    for u_string, v_string in base["verified_affine_points"]:
        u = pari(u_string) * one
        v = pari(v_string) * one
        if v**2 != u**6 - 1:
            raise ArithmeticError("The input point fails the C6 equation.")
        xs = pari.nfroots(nf, horizontal - u)
        ys = pari.nfroots(nf, vertical - v)
        for x in xs:
            for y in ys:
                U = horizontal.subst(t, x) * one
                V = vertical.subst(t, y) * one
                if (U, V) != (u, v) or U**6 - V**2 != 1:
                    raise ArithmeticError("An affine fibre does not lie above the claimed point.")
                affine.add((str(pari.lift(x)), str(pari.lift(y))))
    roots = pari.nfroots(nf, t**infinity_order - 1)
    infinity = {str(pari.lift(root)) for root in roots if (one * root)**infinity_order == 1}
    if len(infinity) != len(roots):
        raise ArithmeticError("An infinity label failed its root-of-unity check.")
    return {
        "field": str(field),
        "status": base["status"],
        "verified_affine_points": sorted(affine),
        "infinity_labels": sorted(infinity),
        "point_count": len(affine) + len(infinity),
    }


def field_and_base(coefficients, base):
    field = pari.Polrev(coefficients, w)
    if base is None:
        from solve_field import classify
        base = classify(coefficients)
    return field, base


def power_cover(coefficients, m, n, base=None):
    positive_integer(m)
    positive_integer(n)
    field, base = field_and_base(coefficients, base)
    a = 2 * m
    b = 6 * n
    h = gcd(a, b)
    result = lift_fibres(field, base, t**n, t**m, h)
    result.update(m=m, n=n, genus=1 + (a*b - a - b - h)//2,
                  infinity_coordinate=f"y^{a//h}/x^{b//h}",
                  geometrically_hyperelliptic=(m == 1))
    return result


def plane_cover(coefficients, n, a, b, base=None):
    positive_integer(n)
    if n < 2:
        raise ValueError("The two-parameter plane family uses n at least two.")
    a, b = Fraction(a), Fraction(b)
    if max(abs(a), abs(b)) > Fraction(1, 2):
        raise ValueError("The smoothness proof applies to rational parameters in [-1/2,1/2].")
    field, base = field_and_base(coefficients, base)
    horizontal = t**n + pari(str(a))*t
    vertical = t**(3*n) + pari(str(b))*t
    degree = 6*n
    result = lift_fibres(field, base, horizontal, vertical, degree)
    result.update(n=n, a=str(a), b=str(b), degree=degree, genus=(degree-1)*(degree-2)//2,
                  infinity_coordinate="Y/X", geometrically_hyperelliptic=False)
    return result


def verify():
    rational = [0, 1]
    fields = [rational, [1, 0, 1], [-7, 0, 1], [7, 0, 1], [-17, 0, 1], [1, 0, -1, 0, 1]]
    from solve_field import classify
    bases = {tuple(f): classify(f) for f in fields}
    tests = 0
    for f in fields:
        base = bases[tuple(f)]
        for m, n in ((1,1),(1,2),(2,1),(2,2),(3,1),(4,3),(5,4)):
            r = power_cover(f, m, n, base)
            if r["status"] != "complete":
                raise ArithmeticError("Unexpected partial input.")
            if f in ([-7,0,1], [7,0,1]):
                predicted = 8 if (m,n) == (1,1) else 4
            else:
                roots_order = {tuple(rational):2,(1,0,1):4,(-17,0,1):2,(1,0,-1,0,1):12}[tuple(f)]
                v = gcd(2*m, roots_order)
                predicted = (gcd(6*n,roots_order) + gcd(gcd(2*m,6*n),roots_order)
                             + (v if (roots_order//v)%2==0 else 0))
            if r["point_count"] != predicted:
                raise ArithmeticError((f, m, n, predicted, r))
            tests += 1
    expected = {tuple(rational):2, (-17,0,1):4, (-7,0,1):2}
    for f, count in expected.items():
        r = plane_cover(list(f), 2, "1/2", "1/2", bases[f])
        if (r["genus"], r["point_count"], r["status"]) != (55,count,"complete"):
            raise ArithmeticError(r)
        tests += 1
    partial = dict(bases[tuple(rational)], status="partial")
    if plane_cover(rational,2,"1/2","1/2",partial)["status"] != "partial":
        raise ArithmeticError("A partial point list was promoted to complete.")
    # The genus formula agrees independently with adjunction on smooth plane Fermat curves.
    for n in range(1, 51):
        m = 3*n
        a, b = 2*m, 6*n
        h = gcd(a,b)
        if 1+(a*b-a-b-h)//2 != (b-1)*(b-2)//2:
            raise ArithmeticError("Genus formula on the diagonal.")
    print(f"{tests} complete field/curve lists, a partial-list guard, and 50 genus comparisons passed.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    actions = parser.add_subparsers(dest="action", required=True)
    actions.add_parser("verify")
    powers = actions.add_parser("powers")
    powers.add_argument("field", help="JSON integer coefficients in ascending order")
    powers.add_argument("m", type=int)
    powers.add_argument("n", type=int)
    plane = actions.add_parser("plane")
    plane.add_argument("field")
    plane.add_argument("n", type=int)
    plane.add_argument("a")
    plane.add_argument("b")
    args = parser.parse_args()
    if args.action == "verify":
        verify()
    else:
        f = json.loads(args.field)
        answer = (power_cover(f,args.m,args.n) if args.action == "powers"
                  else plane_cover(f,args.n,args.a,args.b))
        print(json.dumps(answer, indent=2))
