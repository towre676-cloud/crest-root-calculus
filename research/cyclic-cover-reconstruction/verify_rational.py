"""Independent identities and completeness checks for branch-supported poles."""

from copy import deepcopy

import sympy as sp

from rational import X, coefficients, compute, polynomial, verify


def poly(expr):
    return sp.Poly(expr, X, domain=sp.QQ)


def replay(branch, k, maximum, denominator=None):
    values = coefficients(branch)
    pole_values = coefficients(denominator) if denominator is not None else None
    report = compute(values, k, k, None, maximum, pole_values)
    verify(values, k, report)
    rows = []
    for chart in report["charts"]:
        power = chart["vertical_input_power"]
        source_denominator = polynomial(chart["vertical_denominator"])
        for row in chart["strata"]:
            q = polynomial(row["denominator"])
            for record in row["maps"]:
                a = polynomial(record["numerator"])
                g = polynomial(record["outer"])
                v = polynomial(record["vertical"])
                d = polynomial(record["vertical_pole_denominator"])*source_denominator
                left = branch**power*v**k*q**g.degree()
                right = sum(g.nth(i)*a**i*q**(g.degree()-i) for i in range(g.degree()+1))*d**k
                assert left == right
            rows.append(row)
    return report, rows


def independent_cubic(branch):
    """Eliminate the original homogeneous coefficient identity, without derivatives."""
    b, c, a4, a3, a2, a0, g1, g0 = sp.symbols("b c a4 a3 a2 a0 g1 g0")
    a = X**5+a4*X**4+a3*X**3+a2*X**2+a0
    v = X**2+b*X+c
    equation = sp.Poly(a**2+g1*a*X+g0*X**2-branch.exquo(poly(X)).as_expr()*v**3, X)
    assignments = {}
    for variable, power in zip((a4, a3, a2, g1, a0, g0), (9, 8, 7, 6, 5, 2), strict=True):
        coefficient = sp.expand(equation.nth(power).subs(assignments))
        divisor = coefficient.coeff(variable)
        assert divisor in (1, 2)
        assignments[variable] = -coefficient.subs(variable, 0)/divisor
    residuals = [sp.expand(term.subs(assignments)) for term in equation.all_coeffs()]
    return list(sp.groebner(residuals, b, c, domain=sp.QQ))


def main():
    q = poly(X)
    hyperelliptic = poly(X*(X**2+1)*(X**4+X**2+1))
    _, rows = replay(hyperelliptic, 2, 2)
    assert [(row["degrees"], row["geometric_points"]) for row in rows] == [
        ([2, 3, 0], 1), ([2, 5, 2], 1)]
    assert all(row["maps"][0]["numerator"] == ["1", "0", "1"] for row in rows)

    cubic = polynomial(["0", "32/9", "-34/3", "1060/81", "-55/9", "1"])
    report, rows = replay(cubic, 3, 5, q)
    assert [(row["degrees"], row["geometric_points"]) for row in rows] == [
        ([2, 2, 0], 0), ([5, 2, 2], 1)]
    assert rows[1]["factor_parameters"] == 1 and rows[1]["uses_complement"]
    assert rows[1]["maps"][0]["vertical"] == ["2", "-3", "1"]
    assert independent_cubic(cubic) == [sp.Symbol("b")+3, sp.Symbol("c")-2]
    negative = cubic+q
    _, negative_rows = replay(negative, 3, 5, q)
    assert len(negative_rows) == 2 and all(row["geometric_points"] == 0 for row in negative_rows)
    assert independent_cubic(negative) == [1]

    for denominator, numerator, removed, expected_vertical in (
        (X*(X-1)**2, X**4+2, X*(X-1)**2, X*(X-1)**2),
        (X**4, X**5+2, X, X**3),
    ):
        a, den, r = map(poly, (numerator, denominator, removed))
        source = r*(a**2-den**2)
        _, rows = replay(source, 3, a.degree(), den)
        matches = [record for row in rows for record in row["maps"]]
        assert sum(row["geometric_points"] for row in rows) == 1 and len(matches) == 1
        record = matches[0]
        assert polynomial(record["numerator"]) == a
        assert polynomial(record["vertical_pole_denominator"]) == poly(expected_vertical)

    removed_stratum = deepcopy(report)
    removed_stratum["charts"][0]["strata"].pop()
    try:
        verify(coefficients(cubic), 3, removed_stratum)
    except (ArithmeticError, ValueError):
        pass
    else:
        raise AssertionError("A report with a missing stratum was accepted.")

    print("Verified hyperelliptic, cubic, repeated-pole and negative cases.")


if __name__ == "__main__":
    main()
