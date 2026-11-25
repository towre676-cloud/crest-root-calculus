"""Replay the exact equation solver and independent finite-torsion checks."""

import argparse
from itertools import product
import json
from math import gcd
from pathlib import Path

import sympy as sp

from reconstruct import kernel_order, leading_constant, reconstruct, spectrum, x


HERE = Path(__file__).resolve().parent
FIXED = {"identity": 3, "transposition": 1, "3-cycle": 0}
MATRICES = {3: ((1, 0), (0, 1)), 1: ((0, 1), (1, 0)), 0: ((0, 1), (1, 1))}
EVEN_FIXTURES = (
    (2, 1, [-328483, -769986, -837355, -516180, -186645, -36386, -2717]),
    (2, 0, [-27491309912, -46078252758, -24337164815, -2901666120, 838147080, 112025216, -7483696]),
    (4, 1, [-478289, -909486, -675605, -238140, -33215, 1602, 693]),
    (4, 0, [90733720, 193787046, 182394583, 95708424, 29349432, 4969856, 363440]),
)


def reference_order(a, b, k, matrix):
    """Count points in the graph condition in Smith coordinates of the isogeny."""
    n = 2*gcd(a, k)
    count = 0
    for p, q in product(range(n), repeat=2):
        v = (2*a*p, 2*a*q)
        w = (2*b*p, 2*b*k*q)
        if any(value % n for value in v+w):
            continue
        image = tuple((i*v[0]//n+j*v[1]//n) % 2 for i, j in matrix)
        if image == tuple((value//n) % 2 for value in w):
            count += 1
    return count


def check_finite_law():
    for k in (1, 2, 3, 4, 5, 7, 8, 9, 11, 12, 13, 15, 18, 27):
        choices = MATRICES if k % 2 else {1: MATRICES[3], 0: MATRICES[1]}
        n = 2*k
        factor = sp.prod(1-sp.Rational(1, p*p) for p in sp.factorint(n))
        for fixed_count, matrix in choices.items():
            total = 0
            for a, b in product(range(n), repeat=2):
                if gcd(gcd(a, b), n) != 1:
                    continue
                order = reference_order(a, b, k, matrix)
                total += order
                # Add multiples of n to construct coprime representatives of
                # the same residue. With a=0 use n as its positive lift.
                aa, bb = a or n, b
                while gcd(aa, bb) != 1:
                    bb += n
                assert order == kernel_order(aa, bb, k, fixed_count)
            assert 3*total/(2*n*n*factor) == leading_constant(k, fixed_count)
    print("Verified the degree and leading-constant laws by finite enumeration for 14 isogeny degrees.")


def check_equations(gp):
    original = json.loads((HERE / "reconstructions.json").read_text(encoding="utf-8"))
    for case in original:
        report = reconstruct(case["f"], 400, gp)
        expected = case["expect"]
        assert len(report["quotients"]) == len(case["quotients"])
        for result, old in zip(report["quotients"], case["quotients"]):
            assert result["status"] == "certified"
            assert (result["k"], result["pi_type"]) == (expected["k"], expected["type"])
            assert result["fixed_nonzero_two_torsion"] == FIXED[expected["type"]]
            assert result["spectrum"] == {int(d): m for d, m in old["spectrum"].items()}
            assert result["spectrum"][2] == 2*report["geometric_reduced_involutions"]
        algebraic = report["algebraic_reduced_involutions"]
        f = sp.Poly.from_list(case["f"][::-1], x)
        T = sp.Symbol("T")
        known = report["rational_reduced_involutions"]
        for item in algebraic:
            field = sp.Poly(sp.sympify(item["minimal_polynomial"]), T, domain=sp.QQ)
            assert field.is_irreducible
            known += field.degree()
            A, B, C = map(sp.sympify, item["coefficients"])
            delta = A*A+B*C
            assert field.gcd(sp.Poly(delta, T)).degree() == 0
            identity = sp.Poly(sum(f.nth(i)*(A*x+B)**i*(C*x-A)**(6-i) for i in range(7))
                               - delta**3*f.as_expr(), x)
            assert all(sp.Poly(coefficient, T).rem(field).is_zero
                       for coefficient in identity.all_coeffs())
        assert known == report["geometric_reduced_involutions"]
        print(case["name"], known, "geometric reduced involutions")

    for k, fixed_count, coefficients in EVEN_FIXTURES:
        report = reconstruct(coefficients, 400, gp)
        assert report["rational_reduced_involutions"] == 1
        result, = report["quotients"]
        assert (result["k"], result["fixed_nonzero_two_torsion"]) == (k, fixed_count)
        assert result["spectrum"][2] == 2*report["geometric_reduced_involutions"]
        assert result["spectrum"] != spectrum(k, 1-fixed_count, 400)
        print("Even reconstruction", k, fixed_count, list(result["spectrum"].items())[:5])

    f = sp.Poly.from_list(original[0]["f"][::-1], x)
    quintic = sp.Poly(sp.cancel(x**6*f.as_expr().subs(x, 1+1/x)), x)
    assert quintic.degree() == 5
    report = reconstruct([int(c) for c in quintic.all_coeffs()[::-1]], 100, gp)
    assert len(report["quotients"]) == 2
    assert all(q["spectrum"] == spectrum(1, 1, 100) for q in report["quotients"])

    conjugate = reconstruct([8, -32, 46, -42, 23, -8, 1], 6, gp)
    assert conjugate["quotients"][0]["status"] == "geometric elliptic factors are nonisogenous"


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp", default="gp")
    args = parser.parse_args()
    check_finite_law()
    check_equations(args.gp)
    print("All exact bielliptic reconstruction checks passed.")
