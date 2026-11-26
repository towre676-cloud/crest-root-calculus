"""Replay the cyclic cubic model, quotient classes and multiplicities."""

import argparse
from pathlib import Path
import sys

import sympy as sp

from cubic_genus_two import compute, degree_counts
from reconstruct import spectrum


def verify_case(source):
    x, t = sp.symbols("x t")
    values = list(map(str, sp.Poly(source, x).all_coeffs()[::-1]))
    result = compute(values, 100)
    chart = result["chart"]

    def poly(values):
        return sum(sp.Rational(c)*x**i for i, c in enumerate(values))

    A, C = (poly(chart[key]) for key in ("A", "C"))
    den = poly(chart["xi_denominator"])
    assert sp.expand(source**chart["source_y_power"]*C-den**3*A) == 0
    P, Q, R = (sp.Rational(chart["sextic_coefficients"][i]) for i in (6, 3, 0))
    eta = sp.Rational(chart["eta_constant"]) + t*poly(chart["eta_xi_cubed_coefficient"])
    c2 = sp.Poly(C, x).nth(2)
    assert sp.expand(eta**2-(P*t*t+Q*t+R)-4*c2*t*(t*C-A)) == 0
    return result


def verify_rational_absence(cases):
    """Independent finite-field obstruction to a rational elliptic factor."""
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "exact-completeness"))
    from absence import good_model
    from jac import weil_poly

    for case, prime in zip(cases, (5, 7, 5), strict=True):
        coeffs = map(sp.Rational, case["chart"]["sextic_coefficients"])
        model = [int(c.p)*pow(int(c.q), -1, prime) % prime for c in coeffs]
        assert good_model(model, prime)
        polynomial = sp.Poly.from_list(weil_poly(model, prime), sp.Symbol("T"), domain=sp.QQ)
        assert polynomial.is_irreducible
        print(f"No rational factor for this model at p={prime}. Weil polynomial {polynomial.as_expr()}.")


def main(finite_fields=False):
    x = sp.Symbol("x")
    first = verify_case((1-2*x)*(x*x-1)**2)
    second = verify_case((1-2*x)**2*(x*x-1)*(x+3)**3)
    third = verify_case((1-2*x)*(x*x-2)**2)
    nonsquare = verify_case((1-2*x)*(x*x-8)**2)
    cyclotomic = verify_case((1-2*x)*(x*x+sp.Rational(1,27))**2)
    cm = verify_case((-2*x)*(x*x-1)**2)
    assert set(first["quotient_j"]) == set(second["quotient_j"]) == {"6912", "790272"}
    assert [first["chart"]["source_y_power"], second["chart"]["source_y_power"]] == [1, 2]
    for result in (first, second, third, nonsquare, cyclotomic):
        assert result["status"] == "certified"
        assert result["minimal_isogeny_degree"] == result["fixed_nonzero_two_torsion"] == 3
        assert result["geometric_optimal_counts"] == degree_counts(100)
        assert result["equivariant_elliptic_maps"] is False
    assert third["non_cm_witness"]["test"] == "nonintegral j"
    assert [case["rational_optimal_counts"] for case in (first, second, third, nonsquare, cyclotomic)] == [
        {2: 2}, {2: 2}, {}, {}, {}]
    R = sp.Matrix([[0, -1], [1, -1]])
    I = sp.Matrix([[0, 1], [1, 0]])
    W = R-R**2
    assert R**2+R+sp.eye(2) == sp.zeros(2) and I*R*I == R**2
    assert W*I == -I*W and W**2 == -3*sp.eye(2) and (W*I)**2 == 3*sp.eye(2)
    assert cm["status"] == "non-CM not certified" and cm["quotient_j"] == ["54000"]*2
    assert degree_counts(10_000) == spectrum(3, 3, 10_000)
    if finite_fields:
        verify_rational_absence((third, nonsquare, cyclotomic))
    print("Verified cyclic cubic normal forms, non-CM spectra and the Eisenstein formula through degree 10,000.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--finite-fields", action="store_true",
                        help="Also check the rational absence examples using Python-FLINT.")
    main(parser.parse_args().finite_fields)
