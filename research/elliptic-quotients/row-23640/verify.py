"""Verify the elliptic degree spectra of a family through Sutherland row 23640.

Run with Python 3.11+ and SymPy. All comparisons use exact arithmetic.
"""
import json
from pathlib import Path

import sympy as sp

from discover import discover
from verify_involution import count_curves, require, verify_equations, verify_reductions


HERE = Path(__file__).resolve().parent


def verify_family():
    u, v, t, k, w = sp.symbols("u v t k w")
    a = u**2 + 7*u*v + 2*v**2
    b = u*v*(u**2 + 3*u*v + v**2)
    delta = a**2 - 12*b
    a_k = a + 138*k*u*v
    b_k = b + 23*k*u*v*a + 1587*k**2*u**2*v**2
    quartic = 3*t**4 + a_k*t**2 + b_k

    require(sp.expand((6*t**2 + a_k)**2 - delta - 12*quartic) == 0,
            "The family does not have the claimed fixed elliptic quotient")
    change = sp.Poly(quartic - (3*t**4 + a*t**2 + b), u, v, t, k)
    require(all(c % 23 == 0 for c in change.coeffs()), "Reduction varies in the family")
    require(quartic.subs({u: 1, v: 0, t: 0}) == 0, "Missing rational fixed point")

    # The homogeneous gradient ideal has empty projective zero locus mod 23.
    coords = (u, v, t)
    gradients = [quartic.subs(k, 0).diff(x) for x in coords]
    for x in coords:
        basis = sp.groebner([g.subs(x, 1) for g in gradients],
                            *(y for y in coords if y != x), modulus=23)
        require(list(basis) == [1], f"Smoothness check failed on chart {x}=1")

    # The elliptic target's quartic-to-cubic change, on v=1.
    q = u**2 + u + 8
    require(sp.expand(delta.subs(v, 1) - (q**2 - 60)) == 0, "Quartic model mismatch")
    h = 2*(w + q)
    v_new = h*(2*u + 1)
    diff = v_new**2 - h*(h-15)*(h-16)
    require(sp.rem(diff, w**2 - q**2 + 60, w) == 0, "Cubic model mismatch")
    return {"base_parameter": "integer k", "good_reduction": 23,
            "fixed_point_u_v_t": [1, 0, 0], "elliptic_model": "Y^2=X*(X-1)*(X+15)"}


def verify_absolute_simplicity():
    p = 23
    counts = [count_curves(p, d) for d in (1, 2, 3)]
    c_counts, e_counts = map(list, zip(*counts))
    require(c_counts == [30, 608, 12186], "Unexpected source counts")
    trace = p + 1 - e_counts[0]
    require(trace == 0 and e_counts == [24, 576, 12168], "Unexpected quotient counts")
    s1, s2, s3 = [p**i + 1 - count for i, count in enumerate(c_counts, 1)]
    e2 = sp.Rational(s1**2 - s2, 2)
    e3 = sp.Rational(s1**3 - 3*s1*s2 + 2*s3, 6)
    T, U = sp.symbols("T U")
    jacobian = T**6 - s1*T**5 + e2*T**4 - e3*T**3 + p*e2*T**2 - p**2*s1*T + p**3
    elliptic = T**2 - trace*T + p
    prym = T**4 + 6*T**3 + 34*T**2 + 138*T + 529
    require(sp.expand(jacobian - elliptic*prym) == 0, "Unexpected Prym polynomial")

    # The nontrivial roots of this resultant are ratios of distinct Frobenius roots.
    first = 23*U**4 + 12*U**3 + 14*U**2 + 12*U + 23
    second = 529*U**4 + 736*U**3 + 558*U**2 + 736*U + 529
    eliminant = sp.resultant(prym, prym.subs(T, U*T), T)
    require(sp.expand(eliminant - 23**4*(U-1)**4*first**2*second) == 0,
            "Root-ratio resultant mismatch")
    for expr, var in ((prym, T), (first, U), (second, U)):
        poly = sp.Poly(expr, var)
        reduced = sp.Poly(expr, var, modulus=11)
        require(poly.content() == 1 and reduced.degree() == 4 and reduced.is_irreducible,
                "Irreducibility witness modulo 11 failed")
    require(all(sp.LC(sp.Poly(f, U)) > 1 for f in (first, second)),
            "The root-ratio polynomials need to be primitive and nonmonic")
    return {"p": p, "C_counts": c_counts, "E_counts": e_counts,
            "elliptic_frobenius": str(elliptic), "prym_frobenius": str(prym),
            "irreducibility_prime": 11, "root_ratio_factors": [str(first), str(second)]}


def verify_rational_isogenies():
    witnesses = {}
    primes = (7, 17, 23)
    odd_mazur_primes = (3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163)
    for p in primes:
        require((15*16) % p != 0, "The cubic reduction is singular")
        squares = [0]*p
        for y in range(p):
            squares[y*y % p] += 1
        count = 1 + sum(squares[x*(x-1)*(x+15) % p] for x in range(p))
        trace = p + 1 - count
        delta = trace**2 - 4*p
        excluded = [ell for ell in odd_mazur_primes
                    if p != ell and pow(delta, (ell-1)//2, ell) == ell-1]
        witnesses[p] = {"trace": trace, "discriminant": delta, "excludes": excluded}
    covered = {ell for data in witnesses.values() for ell in data["excludes"]}
    require(covered == set(odd_mazur_primes), "Odd rational isogeny sieve incomplete")
    return witnesses


def main():
    data = json.loads((HERE / "quartic.json").read_text(encoding="utf-8"))
    x, y, z = sp.symbols("x y z")
    # Input text is compared with a rendering of the monomials, never evaluated.
    source = sum(c*x**i*y**j*z**k for i, j, k, c in data["terms"])
    printed = str(source).replace(" ", "").replace("**", "^")
    require(printed == data["quartic"], "Source equation and term list disagree")
    catalogue = HERE.parents[2] / "applications" / "atlas" / "example.json"
    if catalogue.is_file():
        cases = json.loads(catalogue.read_text(encoding="utf-8"))["curves"]
        row = next(case for case in cases if case.get("row") == 23640)
        require(all(data[key] == value for key, value in row.items()),
                "Saved row differs from the atlas input")
    found = discover(data, height=1)
    require(len(found) == 1 and found[0]["center"] == (1, 1, -1)
            and found[0]["axis"] == (0, 1, 0), "Height-one discovery failed")
    proof = {"row": data["row"], "height_one_discovery": found, "quotient": verify_equations(data["terms"]),
             "family": verify_family(), "prym": verify_absolute_simplicity(),
             "independent_row_check": verify_reductions(),
             "rational_isogeny_sieve": verify_rational_isogenies(),
             "geometric_map_degrees": "2*n (n>=1)", "geometric_degree_three_classes": 0,
             "rational_map_degrees": "2*n^2 or 4*n^2 (n>=1)"}
    print(json.dumps(proof, indent=2))


if __name__ == "__main__":
    main()
