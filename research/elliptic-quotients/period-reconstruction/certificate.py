"""Exact checks and vertical coordinates for the period reconstruction receipts."""

import sympy as sp


x = sp.Symbol("x")


def certify_record(coefficients, record):
    """Return the map's degree, coefficient field, differential and Y/y."""
    variables = {"x": x, "w": sp.I}
    gaussian = sp.QQ_I

    def poly(expression):
        return sp.Poly(sp.sympify(expression, locals=variables), x, domain=gaussian)

    numerator, denominator = (poly(record[key]) for key in ("N", "D"))
    g2, g3 = (poly(record[key]).nth(0) for key in ("g2", "g3"))
    linear = record["differential"].removeprefix("(").removesuffix(") x dx/y")
    a, b = (poly(value).nth(0) for value in linear.split(") dx/y + ("))
    differential = poly(a+b*x)
    f = sp.Poly.from_list(coefficients[::-1], x, domain=gaussian)

    if denominator.is_zero or numerator.gcd(denominator).degree() != 0:
        raise ValueError("Use a reduced horizontal rational function.")
    degree = max(numerator.degree(), denominator.degree())
    if degree < 1 or differential.is_zero or g2**3-27*g3**2 == 0:
        raise ValueError("The map must be nonconstant with nonsingular elliptic target.")
    derivative = numerator.diff()*denominator-numerator*denominator.diff()
    cubic = 4*numerator**3-g2*numerator*denominator**2-g3*denominator**3
    if f*derivative**2 != differential**2*denominator*cubic:
        raise ValueError("The vertical coordinate fails exact substitution.")
    vertical = sp.cancel(derivative.as_expr()/(differential*denominator**2).as_expr())
    scalars = [a, b, g2, g3, *numerator.all_coeffs(), *denominator.all_coeffs()]
    field = "Q" if all(value.is_Rational for value in scalars) else "Q(i)"
    slope = (sp.Integer(1), sp.cancel(b/a)) if a else (sp.Integer(0), sp.Integer(1))
    return degree, field, slope, vertical
