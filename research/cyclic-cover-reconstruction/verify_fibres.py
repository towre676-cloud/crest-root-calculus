"""Check bounded fibre searches, open conditions and original map identities."""

from copy import deepcopy

import sympy as sp

from cyclic import X, coefficients, polynomial
from fibres import certified_stratum, compute, plan, verify
from rational import complementary


def poly(expression):
    return sp.Poly(expression, X, domain=sp.QQ)


def verify_map(source, chart, row, record, exponent):
    numerator = polynomial(record["numerator"])
    denominator = polynomial(record["denominator"])
    target = polynomial(record["outer"])
    multiplier = polynomial(record["vertical"])
    vertical_denominator = (polynomial(record["vertical_pole_denominator"])
                            * polynomial(chart["vertical_denominator"]))
    k = chart["target_exponent"]
    source_power = chart["vertical_input_power"]*k//exponent
    source_side = source**source_power * multiplier**k * denominator**target.degree()
    target_side = sum(target.nth(i)*numerator**i*denominator**(target.degree()-i)
                      for i in range(target.degree()+1)) * vertical_denominator**k
    if source_side != target_side or numerator.gcd(denominator).degree() != 0:
        raise ArithmeticError("A replayed quotient fails its original-coordinate identity.")
    assert numerator.degree() == row["numerator_degree"]
    assert row["degree"] == (exponent//k)*max(numerator.degree(), denominator.degree())
    if row["infinity_image"] == "root":
        assert target.nth(0) == 0 and numerator.degree() < denominator.degree()
    else:
        assert numerator.degree() > denominator.degree()


def replay(source, exponent, bound, **kwargs):
    vals = coefficients(source)
    report = compute(vals, exponent, bound, **kwargs)
    verify(vals, exponent, report)
    outputs = []
    for chart in report["charts"]:
        for row in chart["strata"]:
            for result in row["maps"]:
                verify_map(source, chart, row, result, exponent)
                outputs.append((row, result))
    return report, outputs


def one_map(outputs, numerator, denominator, outer, vertical):
    wanted = tuple(coefficients(poly(expr)) for expr in (numerator, denominator, outer, vertical))
    names = ("numerator", "denominator", "outer", "vertical")
    (row, result), = [(row, result) for row, result in outputs
                      if tuple(result[key] for key in names) == wanted]
    return row, result


def targeted(source, exponent, bound, predicate, **kwargs):
    _, charts = plan(coefficients(source), exponent, bound, exponent, 1, "both", **kwargs)
    chart, branch, strata = charts[0]
    s, = filter(predicate, strata)
    row = certified_stratum(branch, exponent, 1, s)
    assert row == certified_stratum(branch, exponent, 1, s, supplied=row["certificate"], check=True)
    for result in row["maps"]:
        verify_map(source, chart, row, result, exponent)
    return [(row, result) for result in row["maps"]]


def independent_unbranched_identity():
    h, a0, a1, a2, g0, g1 = sp.symbols("h a0 a1 a2 g0 g1")
    A, Q = X**4+a2*X*X+a1*X+a0, (X+h)**3
    source = (X**4+X+2)**2-X**6
    equation = sp.Poly(A*A+g1*A*Q+g0*Q**2-source, X)
    basis = sp.groebner(equation.all_coeffs(), g1, g0, a2, a1, a0, h,
                        domain=sp.QQ, order="grevlex")
    assert list(basis) == [g1, g0+1, a2, a1-1, a0-2, h]


def complementary_control():
    # Doubling on the elliptic target produces an order-six unbranched pole.
    t = X+1
    q = poly(t**6)
    f = poly(t**4-1)
    a = poly(t**8+18*t**4-27)-28*q
    v = poly(t**4-9)
    answer = complementary(f, q, poly(t**4), 3, 8, 2, 4,
                           [sp.QQ(1), sp.QQ(6), sp.QQ(4), sp.QQ(-8)], sp.QQ)
    assert answer == (a, v, poly(0))
    assert (a+28*q)**2-64*q**2 == f*v**3


def non_equivariant_control():
    x, xi, eta = sp.symbols("x xi eta")
    t = (1-2*x)/(x*x-1)
    assert sp.cancel((1+x*t)**2-t*t-t-1) == 0
    horizontal = xi+1/xi
    multiplier = (xi+1)/xi**2
    assert sp.cancel((xi**6+xi**3+1)*multiplier**2
                     -(horizontal+2)*(horizontal**3-3*horizontal+1)) == 0
    w = sp.Symbol("w")
    branch = (w+2)*(w**3-3*w+1)
    cubic = sp.cancel(w**4*branch.subs(w, -1/w-2))
    assert cubic == w**3+9*w**2+6*w+1
    j = sp.Rational(1728*4*(-21)**3, 4*(-21)**3+27*37**2)
    assert j == 790272
    assert sp.cancel(multiplier.subs(xi, 1/xi)/xi**3-multiplier) == 0


def main():
    a, q = X**4+X+2, X**3
    source = poly(a*a-q*q)
    report, rows = replay(source, 3, 4, target=3)
    assert len(rows) == 3 and len(report["charts"]) == 2
    row, result = one_map(rows, a, q, X**2-1, 1)
    assert row["infinity_image"] == "pole" and result["unbranched_pole_factor"] == ["0", "1"]

    missing = deepcopy(report)
    missing["charts"][0]["strata"].pop()
    try:
        verify(coefficients(source), 3, missing)
    except (ArithmeticError, ValueError):
        pass
    else:
        raise AssertionError("A report with a missing degree stratum was accepted.")

    _, rows = replay(source+poly(X), 3, 4, target=3, character=1)
    assert rows == []
    _, rows = replay(source*poly((X+7)**6*(X+8)**3), 6, 8, target=3, character=1)
    row, _ = one_map(rows, a, q, X**2-1, 1)
    assert row["degree"] == 8

    q = X*(X-1)**3
    _, rows = replay(poly(X*(1+q)), 3, 4, target=3, character=1)
    row, result = one_map(rows, 1, q, X**2+X, 1)
    assert row["infinity_image"] == "root"
    assert row["branch_pole_denominator"] == ["0", "1"]
    assert result["unbranched_pole_factor"] == ["-1", "1"]

    _, higher = replay(poly(X*(X**2+1)*(X**4+X**2+1)), 2, None, higher_genus_only=True)
    assert len(higher) == 6 and all(row["degree"] == 2 and row["target_degree"] == 5 for row, _ in higher)

    mixed = X*(X-1)*(X-2)**2
    _, rows = replay(poly(mixed.subs(X, X**2)), 3, 2, character=1, target=3, include_mixed=True)
    row, _ = one_map(rows, X**2, 1, mixed, 1)
    assert row["target_profile"] == [[1, 2], [2, 1]]

    q = (X-1)**3
    # A double target root and a double horizontal zero yield a nonconstant V.
    rows = targeted(poly(X*(X**4-q*q)), 3, 3,
                    lambda s: (s.infinity, s.first, s.unbranched_degree,
                               s.target_degree, s.numerator_degree) == ("root", 2, 1, 4, 2)
                    and s.profile == ((1, 2), (2, 1)) and s.poles.horizontal.degree() == 0,
                    include_mixed=True)
    one_map(rows, X**2, q, X**2*(X**2-1), X)

    q = (X-1)**6
    rows = targeted(poly(X*(X+q)), 3, 6,
                    lambda s: (s.infinity, s.target_degree, s.unbranched_degree,
                               s.numerator_degree, s.degree) == ("root", 2, 2, 1, 6))
    one_map(rows, X, q, X**2+X, 1)

    independent_unbranched_identity()
    complementary_control()
    non_equivariant_control()
    print("Verified arbitrary fibres, mixed target branching, subfields, no-map and non-equivariant controls.")


if __name__ == "__main__":
    main()
