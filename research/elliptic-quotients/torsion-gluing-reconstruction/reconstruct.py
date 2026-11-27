"""Exact rational degree-three quotients and the graph of a complementary pair."""
from __future__ import annotations

from dataclasses import dataclass
from itertools import combinations

import sympy as sp
from cypari2 import Pari


x, t = sp.symbols("x t")
pari = Pari()


@dataclass(frozen=True)
class Cover:
    """X = numerator / poles, Y = lift * y / poles**2, with Y**2 = cubic(X)."""

    numerator: sp.Poly
    poles: sp.Poly
    lift: sp.Poly
    cubic: sp.Poly


def integral_primitive(poly):
    return poly.clear_denoms(convert=True)[1].primitive()[1].set_domain(sp.QQ)


def rational_zeros(equations, variables):
    """Extract rational points in a finite affine chart by lex elimination."""
    if not variables:
        return [()] if not any(equations) else []
    basis = sp.groebner(equations, *variables, domain=sp.QQ)
    if list(basis) == [1]:
        return []
    if not basis.is_zero_dimensional:
        raise ValueError("Positive-dimensional fibre chart requires a separate solver.")
    last = variables[-1]
    univariate = sp.Poly(next(p for p in reversed(list(basis)) if p.free_symbols <= {last}), last)
    return [
        partial + (root,)
        for root in univariate.ground_roots()
        for partial in rational_zeros([p.subs(last, root) for p in basis], variables[:-1])
    ]


def cover_identity(source, cover):
    p, d, h, g = cover.numerator, cover.poles, cover.lift, cover.cubic
    cleared = sum(p**k * d**(3-k) * g.nth(k) for k in range(4))
    return (
        p.gcd(d).degree() == 0 and max(p.degree(), d.degree()) == 3
        and g.degree() == 3 and g.gcd(g.diff()).degree() == 0
        and d * cleared == source * h**2
    )


def short_cover(source, p, d):
    """Lift the candidate map, choosing an exact Q-isomorphic short target."""
    complement = source.exquo(d)
    if p.gcd(d).degree() != 0:
        return None
    g = sp.Poly(sp.resultant(complement.as_expr(), t*d.as_expr()-p.as_expr(), x), t).monic()
    if g.degree() != 3 or g.gcd(g.diff()).degree() != 0:
        return None
    cleared = sum(p**k * d**(3-k) * g.nth(k) for k in range(4))
    square, remainder = cleared.div(complement)
    if not remainder.is_zero or square.is_zero:
        return None
    scale, factors = square.factor_list()
    if any(exponent % 2 for _, exponent in factors):
        return None
    h = sp.prod(
        (factor**(exponent//2) for factor, exponent in factors),
        start=sp.Poly(1, x, domain=sp.QQ),
    )
    # Complete the cube and choose short coordinates obtained from a minimal model.
    a = g.nth(2) / 3
    short = sp.Poly(scale**3 * g.as_expr().subs(t, t/scale - a), t)
    _, change = pari.ellinit([pari(str(short.nth(k))) for k in (1, 0)]).ellminimalmodel()
    u = sp.Rational(str(change[0])) / 6
    cover = Cover(
        (p + d*a).mul_ground(scale / u**2),
        d,
        h.mul_ground(scale**2 / u**3),
        sp.Poly(u**(-6) * short.as_expr().subs(t, u**2*t), t),
    )
    if not cover_identity(source, cover):
        raise ArithmeticError("Recovered map failed its polynomial identity.")
    return cover


def discover(source):
    """Enumerate rational degree-three maps modulo target isomorphism and translation.

    The model is squarefree of degree six; all source Weierstrass points are finite.
    Each of the four target branch fibres contains an odd number of Weierstrass
    points. The unique triple makes its cubical denominator rational.
    """
    if source.degree() != 6 or source.gcd(source.diff()).degree() != 0:
        raise ValueError("Expected a squarefree rational sextic.")
    factors = [factor for factor, _ in source.factor_list()[1]]
    a, b, c = coefficients = sp.symbols("a b c")
    numerator = a + b*x + c*x*x
    maps = []
    for size in range(1, len(factors)):
        for subset in combinations(factors, size):
            poles = integral_primitive(sp.prod(subset))
            if poles.degree() != 3:
                continue
            complement = source.exquo(poles).as_expr().subs(x, t)
            d = poles.as_expr()
            fibre = d.subs(x, t)*numerator - numerator.subs(x, t)*d
            quadratic, residue = sp.div(fibre, x-t, x)
            assert residue == 0
            lead, middle, constant = sp.Poly(quadratic, x).all_coeffs()
            equations = sp.Poly(sp.rem(middle**2-4*lead*constant, complement, t), t).all_coeffs()
            for i, coefficient in enumerate(coefficients):
                fixed = {v: 0 for v in coefficients[:i]} | {coefficient: 1}
                free = coefficients[i+1:]
                for root in rational_zeros([e.subs(fixed) for e in equations], free):
                    p = integral_primitive(sp.Poly(numerator.subs(fixed | dict(zip(free, root))), x))
                    cover = short_cover(source, p, poles)
                    if cover is not None:
                        maps.append(cover)
    return maps


def division3(cubic):
    """Monic x-coordinate polynomial of the eight nonzero 3-torsion points."""
    a, b = cubic.nth(1), cubic.nth(0)
    return sp.Poly(t**4 + 2*a*t*t + 4*b*t - a*a/3, t, domain=sp.QQ)


def scale_by(poly, number):
    return sp.Poly.new(poly.rep.mul_ground(number), poly.gen)


def evaluate(poly, number, field):
    return poly.set_domain(field).rep.eval(number)


def odd_chart(poly, branch, weight, field):
    """The polynomial u**weight * poly(branch + 1/u), written in x = u."""
    u = sp.Poly(x, x, domain=field)
    linear = scale_by(u, branch) + 1
    return sum(
        scale_by(linear**i * u**(weight-i), field.convert(poly.nth(i)))
        for i in range(poly.degree()+1)
    )


def reduced_pullback(odd_source, cover, branch, torsion_x, field):
    """Mumford pair (U, w V) for f*(P-O), with w**2 = cubic(torsion_x)."""
    p = odd_chart(cover.numerator, branch, 3, field)
    d = odd_chart(cover.poles, branch, 3, field)
    h = odd_chart(cover.lift, branch, 3, field)
    fibre = (p - scale_by(d, torsion_x)).monic()
    assert fibre.degree() == 3 and fibre.gcd(d).degree() == 0
    assert fibre.gcd(h).degree() == 0
    # The pole divisor is supported at Weierstrass points, which are 2-torsion
    # relative to the Weierstrass infinity. Its negative has the same v=0.
    v = d * (d * h.invert(fibre)).rem(fibre)
    u = (fibre * d).monic()
    g_value = evaluate(cover.cubic, torsion_x, field)
    while u.degree() > 2:
        u = (odd_source - scale_by(v*v, g_value)).exquo(u).monic()
        v = (-v).rem(u)
    assert (odd_source - scale_by(v*v, g_value)).rem(u).is_zero
    return u, v


def pari_polynomial(poly):
    return pari.Polrev([pari(str(c)) for c in reversed(poly.all_coeffs())])


def linear_roots(poly, field):
    for factor, _ in poly.set_domain(field).factor_list()[1]:
        if factor.degree() == 1:
            lead, constant = factor.rep.to_list()
            yield -constant / lead


def with_branch(source, factor, field, alpha):
    """Choose a Weierstrass point after at most one compositum extension."""
    branch = next(linear_roots(source, field), None)
    if branch is not None:
        return field, alpha, branch
    branch_poly = source.factor_list()[1][0][0].monic()
    defining, alpha_image, branch_image, _ = pari.polcompositum(
        pari_polynomial(factor), pari_polynomial(branch_poly), 1
    )[0]

    def as_poly(expression):
        return sp.Poly.from_list(
            [sp.Rational(str(c)) for c in expression.lift().Vec()], t, domain=sp.QQ
        )

    extension = sp.QQ.alg_field_from_poly(as_poly(defining), alias="z")
    alpha = evaluate(as_poly(alpha_image), extension.unit, extension)
    branch = evaluate(as_poly(branch_image), extension.unit, extension)
    assert evaluate(factor, alpha, extension) == extension.zero
    assert evaluate(source, branch, extension) == extension.zero
    return extension, alpha, branch


def descend(value, alpha, degree, field):
    """Express a temporary field element in the original x-coordinate basis."""
    def column(element):
        coefficients = element.to_list()
        return [sp.QQ.to_sympy(c) for c in reversed(coefficients)] + [0] * (
            field.mod.degree()-len(coefficients)
        )

    matrix = sp.Matrix([column(alpha**k) for k in range(degree)]).T
    coefficients, _ = matrix.gauss_jordan_solve(sp.Matrix(column(value)))
    return sp.Poly.from_list(list(reversed(coefficients)), t, domain=sp.QQ)


def recover_graph(source, first, second):
    """Return (factor, A, B) on each component of E1[3]/±.

    A and B describe (X,Y) -> (A(X),B(X)*Y) in the graph of ker(i1+i2).
    """
    if not (cover_identity(source, first) and cover_identity(source, second)):
        raise ValueError("Both covers must satisfy their defining polynomial identities.")
    factors = division3(first.cubic).factor_list()[1]
    second_division = division3(second.cubic)
    result = []
    for factor, _ in factors:
        factor = factor.monic()
        field = sp.QQ.alg_field_from_poly(factor, alias="alpha")
        alpha = field.unit if factor.degree() > 1 else field.convert(-factor.nth(0))
        extension, alpha, branch = with_branch(source, factor, field, alpha)
        odd = odd_chart(source, branch, 6, extension)
        assert odd.degree() == 5
        u1, v1 = reduced_pullback(odd, first, branch, alpha, extension)
        matches = []
        for beta in linear_roots(second_division, field):
            beta_poly = sp.Poly.from_list(beta.to_list(), t, domain=sp.QQ)
            beta = evaluate(beta_poly, alpha, extension)
            u2, v2 = reduced_pullback(odd, second, branch, beta, extension)
            if u1 != u2 or v1.is_zero or v2.is_zero:
                continue
            ratio = -v1.rep.LC() / v2.rep.LC()
            if v1 != scale_by(v2, -ratio):
                continue
            if evaluate(first.cubic, alpha, extension)*ratio**2 == evaluate(second.cubic, beta, extension):
                matches.append((beta_poly, descend(ratio, alpha, factor.degree(), extension)))
        if len(matches) != 1:
            raise ValueError("No unique graph component for the chosen pair of covers.")
        result.append((factor, *matches[0]))
    return result


def coefficients(poly):
    """Rational coefficients in ascending order."""
    return [str(c) for c in reversed(poly.all_coeffs())]


def serialise(covers, components):
    return {
        "maps": [
            {k: coefficients(getattr(cover, k)) for k in ("numerator", "poles", "lift", "cubic")}
            for cover in covers
        ],
        "kernel": [
            dict(zip(("factor", "x", "y_multiplier"), map(coefficients, component)))
            for component in components
        ],
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sextic", type=json.loads, help="JSON coefficients in ascending order")
    parser.add_argument("--graph", action="store_true", help="Recover the graph for two complementary maps")
    args = parser.parse_args()
    source = sp.Poly.from_list([sp.Rational(c) for c in reversed(args.sextic)], x, domain=sp.QQ)
    covers = discover(source)
    if args.graph and len(covers) != 2:
        parser.error("--graph expects exactly two rational degree-three covers")
    components = recover_graph(source, *covers) if args.graph else []
    print(json.dumps(serialise(covers, components), indent=2))
