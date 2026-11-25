"""Find every rational degree-two quotient of a hyperelliptic equation."""

import argparse
import json
from pathlib import Path

import sympy as sp

from discover import X, coefficients, genus, normalization, polynomial
from finite import rational_points, verify as verify_finite


T = sp.Symbol("t")
CENTER, SCALE = sp.symbols("a k")


def source_model(values):
    square, source = normalization(polynomial(values))
    if genus(source) < 2:
        raise ValueError("The normalized source must have genus at least two.")
    return square, source


def affine_center(source):
    n = source.degree()
    if n % 2:
        return None
    center = -source.nth(n - 1) / (n * source.LC())
    return center if source.compose(sp.Poly(2 * center - X, X, domain=sp.QQ)) == source else None


def nonaffine_equations(source):
    half_degree = genus(source) + 1
    shifted = sp.Poly(source.as_expr().subs(X, X + CENTER), X)
    return [shifted.nth(2 * half_degree - j) * SCALE**(half_degree - j) - shifted.nth(j)
            for j in range(half_degree)]


def deck_transformations(source, points):
    center = affine_center(source)
    if center is not None:
        for sign in (1, -1):
            yield 2 * center - X, sp.Integer(sign)
    h = genus(source) + 1
    for a, k in points:
        a, k = sp.Rational(a), sp.Rational(k)
        if not k:
            raise ArithmeticError("A singular projective transformation survived.")
        square = k**h
        if square < 0:
            continue
        numerator, integral_num = sp.integer_nthroot(square.p, 2)
        denominator, integral_den = sp.integer_nthroot(square.q, 2)
        if integral_num and integral_den:
            root = sp.Rational(numerator, denominator)
            for sign in (1, -1):
                yield a + k / (X - a), sign * root / (X - a)**h


def rational_pair(expression, variable):
    numerator, denominator = sp.cancel(expression).as_numer_denom()
    p = sp.Poly(numerator, variable, domain=sp.QQ)
    q = sp.Poly(denominator, variable, domain=sp.QQ)
    return p.quo_ground(q.LC()), q.monic()


def function_record(expression):
    top, bottom = rational_pair(expression, X)
    return {"numerator": coefficients(top), "denominator": coefficients(bottom)}


def descend(horizontal, value):
    """Reduce a fixed rational function in the quadratic algebra over Q(t)."""
    top, bottom = rational_pair(horizontal, X)
    field = sp.QQ.frac_field(T)
    relation = sp.Poly(top.as_expr() - T * bottom.as_expr(), X, domain=field)
    if relation.degree() != 2:
        raise ArithmeticError("The horizontal invariant lacks degree two.")
    p, q = rational_pair(value, X)
    residual = (p.set_domain(field) * q.set_domain(field).invert(relation)).rem(relation)
    if residual.degree() > 0:
        raise ArithmeticError("The vertical square is not in the invariant horizontal field.")
    return rational_pair(residual.nth(0), T)


def quotient(source, square, m, c):
    # The invariant projection of y is zero for c = -1. In that case
    # project xy, which is nonzero because m is nonidentity.
    u = sp.cancel(X + m)
    if not u.has(X):
        u = sp.cancel(X * m)
    factor = sp.cancel(1 + c)
    if not factor:
        factor = sp.cancel(X + c * m)
    if not factor or (sp.cancel(m.subs(X, m) - X) != 0
                      or sp.cancel(c * c.subs(X, m) - 1) != 0
                      or sp.cancel(source.as_expr().subs(X, m) - c**2 * source.as_expr()) != 0
                      or sp.cancel(u.subs(X, m) - u) != 0
                      or sp.cancel(c * factor.subs(X, m) - factor) != 0):
        raise ArithmeticError("The involution or its fixed generators failed an identity.")
    numerator, denominator = descend(u, source.as_expr() * factor**2)
    radicand = sp.Poly((numerator * denominator).as_expr().subs(T, X), X, domain=sp.QQ)
    removed, target = normalization(radicand)
    vertical = sp.cancel(factor * denominator.as_expr().subs(T, u) / removed.as_expr().subs(X, u)
                         / square.as_expr())
    if not 1 <= genus(target) < genus(source):
        raise ArithmeticError("The fixed curve has an inadmissible genus.")
    if sp.cancel((source * square**2).as_expr() * vertical**2 - target.as_expr().subs(X, u)) != 0:
        raise ArithmeticError("The quotient map failed the original equation.")
    return {"involution": {"x": function_record(m), "y_multiplier": function_record(c)},
            "target": coefficients(target), "target_genus": genus(target), "degree": 2,
            "map": {"u": function_record(u), "v_over_input_y": function_record(vertical)}}


def quotient_list(source, square, points):
    return [quotient(source, square, m, c) for m, c in deck_transformations(source, points)]


def recover(values):
    square, source = source_model(values)
    certificate = rational_points(nonaffine_equations(source), [CENTER, SCALE])
    return {"source": coefficients(source), "square": coefficients(square), "source_genus": genus(source),
            "nonaffine_certificate": certificate,
            "quotients": quotient_list(source, square, certificate["rational_points"])}


def verify(values, result):
    square, source = source_model(values)
    if (result["source"], result["square"], result["source_genus"]) != (
            coefficients(source), coefficients(square), genus(source)):
        raise ArithmeticError("The normalized model differs from the source input.")
    certificate = result["nonaffine_certificate"]
    if (certificate["variables"] != list(map(str, (CENTER, SCALE)))
            or certificate["equations"] != list(map(str, nonaffine_equations(source)))):
        raise ArithmeticError("The nonaffine certificate omits or changes chart equations.")
    verify_finite(certificate)
    if result["quotients"] != quotient_list(source, square, certificate["rational_points"]):
        raise ArithmeticError("The quotient list differs from the certified involution list.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON array of coefficients, constant term first")
    parser.add_argument("--verify", type=Path, help="verify a saved certificate for the given input")
    args = parser.parse_args()
    values = json.loads(args.input.read_text(encoding="utf-8"))
    if args.verify:
        verify(values, json.loads(args.verify.read_text(encoding="utf-8")))
        print("Verified all rational degree-two quotients.")
    else:
        print(json.dumps(recover(values), indent=2))


if __name__ == "__main__":
    main()
