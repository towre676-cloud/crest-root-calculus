"""Exact degree-four reconstruction in the paired branch chart of a genus-two sextic."""
from __future__ import annotations

from dataclasses import dataclass
from functools import reduce
from itertools import product
from pathlib import Path
import sys

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "torsion-gluing-reconstruction"))
from reconstruct import (
    coefficients, descend, evaluate, integral_primitive, linear_roots, odd_chart,
    pari, scale_by, t, with_branch, x,
)


@dataclass(frozen=True)
class Cover:
    """A map (N/S², H*y/S³) to the short Weierstrass cubic."""

    numerator: sp.Poly
    half_poles: sp.Poly
    lift: sp.Poly
    cubic: sp.Poly


def quadratic_polys(source):
    if source.degree() != 6 or source.gcd(source.diff()).degree() != 0:
        raise ValueError("Expected a squarefree sextic.")
    factors = [integral_primitive(q) for q, _ in source.factor_list()[1]]
    linear = [p for p in factors if p.degree() == 1]
    quadratics = [p for p in factors if p.degree() == 2]
    if len(linear) == 2 and len(quadratics) == 2:
        return [linear[0]*linear[1], *quadratics]
    if len(quadratics) == 3:
        return quadratics
    raise ValueError("Expected three quadratic branch factors with at most one split pair.")


def branch_values(p, d, quadratics):
    roots = []
    for q in quadratics:
        r = (p * d.invert(q)).rem(q)
        if r.degree() > 0:
            return None
        roots.append(r.nth(0))
    return roots if len(set(roots)) == 3 else None


def certify_identity(source, cover):
    p, d, h, g = cover.numerator, cover.half_poles**2, cover.lift, cover.cubic
    assert p.gcd(d).degree() == 0 and max(p.degree(), d.degree()) == 4
    assert g.nth(3) == 1 and g.nth(2) == 0 and g.gcd(g.diff()).degree() == 0
    assert sum(p**i * d**(3-i) * g.nth(i) for i in range(4)) == source*h*h


def build_map(source, p, s, quadratics):
    d = s*s
    if max(p.degree(), d.degree()) != 4 or s.gcd(source).degree() or p.gcd(d).degree():
        return None
    roots = branch_values(p, d, quadratics)
    if roots is None:
        return None
    g = sp.Poly(sp.prod(t-r for r in roots), t, domain=sp.QQ)
    return build_from_cubic(source, p, s, g)


def build_from_cubic(source, p, s, g):
    d = s*s
    if max(p.degree(), d.degree()) != 4 or s.gcd(source).degree() or p.gcd(d).degree():
        return None
    if g.degree() != 3 or g.nth(3) != 1 or not g.is_sqf:
        return None
    quotient, remainder = sum(p**i*d**(3-i)*g.nth(i) for i in range(4)).div(source)
    if not remainder.is_zero or quotient.is_zero:
        return None
    scalar, factors = quotient.factor_list()
    if any(exponent % 2 for _, exponent in factors):
        return None
    h = sp.prod((q**(m//2) for q, m in factors), start=sp.Poly(1, x, domain=sp.QQ))
    shift = g.nth(2) / 3
    short = sp.Poly(scalar**3*g.as_expr().subs(t, t/scalar-shift), t)
    _, transform = pari.ellinit([pari(str(short.nth(k))) for k in (1, 0)]).ellminimalmodel()
    u = sp.Rational(str(transform[0]))/6
    cover = Cover(
        (p + d*shift).mul_ground(scalar/u**2),
        s,
        h.mul_ground(scalar**2/u**3),
        sp.Poly(short.as_expr().subs(t, u**2*t)/u**6, t, domain=sp.QQ),
    )
    certify_identity(source, cover)
    return cover


def incidence(quadratics, degree):
    """Eliminate the numerator in a normalized rank-three chart."""
    a, b = sp.symbols("a b")
    variables, s = [((), sp.Integer(1)), ((a,), x+a), ((a, b), x*x+a*x+b)][degree]
    d = s*s
    exponents = [i for i in range(5) if i != 2*degree]
    qs = [q.as_expr() for q in quadratics]
    residues = [sp.Poly(sp.rem(d, q, x), x) for q in qs]
    remainders = [[sp.Poly(sp.rem(x**k, q, x), x) for k in exponents] for q in qs]
    matrix = sp.Matrix([
        [sp.expand(r.nth(0)*e.nth(1)-r.nth(1)*e.nth(0)) for r in row]
        for row, e in zip(remainders, residues)
    ])
    minors = [sp.expand((-1)**k*matrix[:, [j for j in range(4) if j != k]].det()) for k in range(4)]
    p = sum(coefficient*x**k for k, coefficient in zip(exponents, minors))
    discriminants = []
    for q, e in zip(qs, residues):
        r = sp.Poly(sp.rem(p, q, x), x)
        divided, remainder = sp.div(e.nth(0)*p-r.nth(0)*d, q, x)
        assert remainder == 0
        h = sp.Poly(divided, x)
        discriminants.append(sp.factor(h.nth(1)**2-4*h.nth(0)*h.nth(2)))
    return variables, s, p, discriminants


def solve_two(first, second, variables):
    """Rational intersections of factors with nonzero resultant."""
    a, b = variables
    for l, r in product(sp.factor_list(first, a, b)[1], sp.factor_list(second, a, b)[1]):
        f, g = l[0], r[0]
        resultant = sp.resultant(f, g, b)
        if resultant == 0:
            raise ValueError("A common component requires a separate incidence chart.")
        for u in sp.Poly(resultant, a).ground_roots():
            common = sp.gcd(sp.Poly(f.subs(a, u), b), sp.Poly(g.subs(a, u), b))
            for v in common.ground_roots():
                yield {a: u, b: v}


def pole_candidates(discriminants, variables):
    if not variables:
        if not any(discriminants):
            yield {}
    elif len(variables) == 1:
        a, = variables
        common = reduce(sp.Poly.gcd, (sp.Poly(d, a) for d in discriminants))
        if common.is_zero:
            raise ValueError("The one-dimensional chart degenerates identically.")
        yield from ({a: root} for root in common.ground_roots())
    else:
        if any(disc == 0 for disc in discriminants):
            raise ValueError("The quadratic chart has an identically zero discriminant.")
        for values in solve_two(*discriminants[:2], variables):
            if discriminants[2].subs(values) == 0:
                yield values


def discover(source):
    quadratics = quadratic_polys(source)
    seen = set()
    covers = []
    for degree in (0, 1, 2):
        variables, s, p, discriminants = incidence(quadratics, degree)
        for values in pole_candidates(discriminants, variables):
            poles = sp.Poly(s.subs(values), x, domain=sp.QQ)
            if poles in seen:
                continue
            numerator = integral_primitive(sp.Poly(p.subs(values), x))
            if numerator.is_zero:
                continue
            cover = build_map(source, numerator, poles, quadratics)
            if cover is not None:
                seen.add(poles)
                covers.append(cover)
    return sorted(covers, key=lambda cover: (cover.half_poles.degree(), cover.half_poles.nth(0)), reverse=True)


def sextic_chart(source):
    """Use the original equation when all six branch points are finite."""
    if source.degree() == 6:
        return source, None
    if source.degree() != 5 or not source.is_sqf:
        raise ValueError("Expected a squarefree quintic or sextic.")
    a = next(i for i in range(6) if source.eval(i) != 0)
    return odd_chart(source, sp.QQ.convert(a), 6, sp.QQ), a


def exact_four_x(cover):
    """x-coordinates of points of exact order four on the short model."""
    g = cover.cubic
    a, b = g.nth(1), g.nth(0)
    return sp.Poly(t**6+5*a*t**4+20*b*t**3-5*a*a*t*t-4*a*b*t-8*b*b-a**3, t,
                   domain=sp.QQ)


def reduced_pullback(odd, cover, branch, coordinate, field):
    p = odd_chart(cover.numerator, branch, 4, field)
    d = odd_chart(cover.half_poles, branch, 2, field)
    h = odd_chart(cover.lift, branch, 3, field)
    u = (p-scale_by(d*d, coordinate)).monic()
    assert u.degree() == 4 and u.gcd(d*h).degree() == 0
    v = (d**3 * h.invert(u)).rem(u)
    value = evaluate(cover.cubic, coordinate, field)
    # The origin's empty branch fibre is twice a hyperelliptic fibre.
    while u.degree() > 2:
        u = (odd-scale_by(v*v, value)).exquo(u).monic()
        v = (-v).rem(u)
    assert (odd-scale_by(v*v, value)).rem(u).is_zero
    return u, v


def recover_kernel(source, first, second):
    """Return the 2-torsion labels and signed points of exact order four."""
    quadratics = quadratic_polys(source)
    for cover in (first, second):
        certify_identity(source, cover)
    labels = [branch_values(c.numerator, c.half_poles**2, quadratics) for c in (first, second)]
    return list(zip(*labels, strict=True)), recover_order_four(source, first, second)


def recover_order_four(source, first, second):
    """Match the signed pullback classes above the exact order-four abscissas."""
    second_poly = exact_four_x(second)
    components = []
    for factor, _ in exact_four_x(first).factor_list()[1]:
        factor = factor.monic()
        field = sp.QQ.alg_field_from_poly(factor, alias="alpha")
        alpha = field.unit if factor.degree() > 1 else field.convert(-factor.nth(0))
        extension, alpha, branch = with_branch(source, factor, field, alpha)
        odd = odd_chart(source, branch, 6, extension)
        first_u, first_v = reduced_pullback(odd, first, branch, alpha, extension)
        matches = []
        for beta in linear_roots(second_poly, field):
            beta_poly = sp.Poly.from_list(beta.to_list(), t, domain=sp.QQ)
            beta = evaluate(beta_poly, alpha, extension)
            second_u, second_v = reduced_pullback(odd, second, branch, beta, extension)
            if first_u != second_u or first_v.is_zero or second_v.is_zero:
                continue
            multiplier = -first_v.rep.LC()/second_v.rep.LC()
            if first_v != scale_by(second_v, -multiplier):
                continue
            if evaluate(first.cubic, alpha, extension)*multiplier**2 != evaluate(second.cubic, beta, extension):
                continue
            matches.append((beta_poly, descend(multiplier, alpha, factor.degree(), extension)))
        if len(matches) != 1:
            raise ValueError("No unique signed graph above an order-four abscissa.")
        components.append((factor, *matches[0]))
    assert sum(factor.degree() for factor, _, _ in components) == 6
    return components


def solve(source):
    source, center = sextic_chart(source)
    covers = discover(source)
    if len(covers) != 2:
        raise ValueError("The incidence must give two complementary covers for kernel reconstruction.")
    two, four = recover_kernel(source, *covers)
    return {
        "reciprocal_chart_center": center,
        "maps": [
            {name: coefficients(getattr(cover, name))
             for name in ("numerator", "half_poles", "lift", "cubic")}
            for cover in covers
        ],
        "order_two_x_pairs": [[str(a), str(b)] for a, b in two],
        "order_four": [
            dict(zip(("factor", "x", "y_multiplier"), map(coefficients, component)))
            for component in four
        ],
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sextic", type=json.loads, help="JSON list of ascending quintic or sextic coefficients")
    args = parser.parse_args()
    source = sp.Poly.from_list([sp.Rational(c) for c in reversed(args.sextic)], x, domain=sp.QQ)
    print(json.dumps(solve(source), indent=2))
