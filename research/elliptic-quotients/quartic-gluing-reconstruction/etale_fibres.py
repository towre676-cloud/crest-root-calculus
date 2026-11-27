"""Rational degree-four quotients with arbitrary Galois action on their Weierstrass pairs."""
from __future__ import annotations

from dataclasses import dataclass
from functools import reduce

import sympy as sp

from reconstruct_four import (
    build_from_cubic, pari, recover_order_four, sextic_chart, t, x,
)
from reconstruct import coefficients, pari_polynomial, rational_zeros

u = sp.Symbol("u")


@dataclass(frozen=True)
class Pairing:
    modulus: sp.Poly
    generator: sp.Poly
    trace: sp.Poly
    norm: sp.Poly


def pair_indices(indices):
    if not indices:
        yield ()
        return
    for j in indices[1:]:
        for rest in pair_indices(tuple(i for i in indices[1:] if i != j)):
            yield ((indices[0], j), *rest)


def lagrange_basis(source):
    """Roots and Lagrange idempotents over the exact splitting field."""
    f = pari_polynomial(source.monic())
    field = pari.nfsplitting(f).subst("x", "y")
    if int(field.poldegree()) == 1:
        roots = [pari(str(r)) for r in source.ground_roots()]
    else:
        roots = list(pari.nfroots(pari.nfinit(field), f))
    assert len(roots) == 6
    idempotents = []
    for r in roots:
        q, remainder = pari.divrem(f, pari("x")-r, "x")
        assert remainder == 0
        idempotents.append(q / f.deriv().subst("x", r))
    return roots, idempotents


def rational_poly(poly):
    """Descend a polynomial in x from the splitting field, or reject it."""
    values = [pari.lift(pari.polcoef(poly, k, "x")) for k in range(6)]
    if any(value.poldegree("y") > 0 for value in values):
        return None
    return sp.Poly.from_list([sp.Rational(str(c)) for c in reversed(values)], x, domain=sp.QQ)


def vector(poly, length):
    return sp.Matrix([poly.nth(i) for i in range(length)])


def paired_algebras(source):
    """Enumerate all Galois-stable partitions into pairs and their cubic algebras."""
    # A Galois-stable pairing embeds the permutation group in a group of order 48.
    if source.is_irreducible:
        order = int(pari.polgalois(pari_polynomial(source.monic()))[0])
        if 48 % order:
            return
    roots, idempotents = lagrange_basis(source)
    for pairs in pair_indices(tuple(range(6))):
        trace = rational_poly(sum((roots[i]+roots[j])*(idempotents[i]+idempotents[j]) for i, j in pairs))
        if trace is None:
            continue
        norm = rational_poly(sum(roots[i]*roots[j]*(idempotents[i]+idempotents[j]) for i, j in pairs))
        if norm is None:
            continue
        for k in range(4):
            generator = trace+norm*k
            powers = [(generator**i).rem(source) for i in range(4)]
            basis = sp.Matrix.hstack(*(vector(p, 6) for p in powers[:3]))
            if basis.rank() == 3:
                break
        else:
            raise ArithmeticError("No generator separates the three pairs.")
        def in_basis(poly):
            c, _ = basis.gauss_jordan_solve(vector(poly, 6))
            return sp.Poly.from_list(list(reversed(c)), u, domain=sp.QQ)
        modulus = sp.Poly(u**3, u, domain=sp.QQ)-in_basis(powers[3])
        yield Pairing(modulus, generator, in_basis(trace), in_basis(norm))


def normalize(poly, variable):
    return sp.Poly(poly, variable, domain=sp.QQ)


def cubic_value(pairing, value):
    return normalize(sp.resultant(pairing.modulus.as_expr(), t-value, u), t)


def finite_rational_points(equations, variables, nonzero):
    """Saturate by the map degeneracies before extracting rational points."""
    inverse = sp.Dummy("inverse")
    polys = [*equations, inverse*nonzero-1]
    basis = sp.groebner(polys, inverse, *variables, domain=sp.QQ, order="grevlex")
    if list(basis) == [1]:
        return []
    if not basis.is_zero_dimensional:
        raise ArithmeticError("The saturated degree-four chart failed the finiteness check.")
    lex = basis.fglm(order="lex")
    relations = [f for f in lex if not f.has(inverse)]
    return rational_zeros(relations, variables)


def bounded_intersection(equations, variables, nonzero):
    """Finite saturated rational solutions in at most two variables."""
    if any(e != 0 and not e.free_symbols for e in equations):
        return []
    if not variables:
        return [()] if nonzero != 0 else []
    eqs = [sp.Poly(e, *variables, domain=sp.QQ) for e in equations if e != 0]
    if not eqs:
        return finite_rational_points(equations, variables, nonzero)
    if len(variables) == 1:
        common = reduce(sp.Poly.gcd, eqs)
        return [(r,) for r in common.ground_roots() if nonzero.subs(variables[0], r)]
    a,b = variables
    unit = sp.Poly(nonzero,*variables,domain=sp.QQ)
    if unit.is_zero:
        return []
    eqs.sort(key=sp.Poly.total_degree)
    found = set()
    for factor,_ in eqs[0].factor_list()[1]:
        if unit.rem(factor).is_zero:
            continue
        other = next((eq for eq in eqs[1:] if not eq.rem(factor).is_zero),None)
        if other is None:
            return finite_rational_points(equations,variables,nonzero)
        for g,_ in other.factor_list()[1]:
            if unit.rem(g).is_zero:
                continue
            res = sp.Poly(sp.resultant(factor.as_expr(),g.as_expr(),b),a)
            if res.is_zero:
                return finite_rational_points(equations,variables,nonzero)
            for first in res.ground_roots():
                remaining = sp.gcd(sp.Poly(factor.as_expr().subs(a,first),b),
                                   sp.Poly(g.as_expr().subs(a,first),b))
                for last in remaining.ground_roots():
                    vals = {a:first,b:last}
                    if nonzero.subs(vals) and not any(eq.as_expr().subs(vals) for eq in eqs):
                        found.add((first,last))
    return sorted(found)


def map_conditions(source, pairing, s, n, value, pole_degree):
    leading = sp.Poly(n,x).nth(5)
    n = sp.expand(n-leading*x**5)
    q = x*x-pairing.trace.as_expr()*x+pairing.norm.as_expr()
    p = sp.Poly(sp.div(n-value*s*s,q,x)[0],x)
    disc = p.nth(1)**2-4*p.nth(0)*p.nth(2)
    coefficients_in_u = sp.Poly(sp.rem(disc,pairing.modulus.as_expr(),u),u).all_coeffs()
    # Three distinct values of b are equivalent to three distinct roots of its characteristic cubic.
    values_disc = sp.discriminant(sp.resultant(pairing.modulus.as_expr(),t-value,u),t)
    conditions = (leading,*coefficients_in_u)
    unit = sp.resultant(s,source.as_expr(),x)*sp.resultant(s,n,x)*values_disc
    if pole_degree < 2:
        unit *= sp.Poly(n,x).nth(4)
    return n,conditions,unit


def pairing_candidates(source, pairing, degree):
    a, b = sp.symbols("a b")
    variables, s = [((), sp.Integer(1)), ((a,), x+a), ((a, b), x*x+a*x+b)][degree]
    d = s*s
    generator = pairing.generator.as_expr()
    n1, n2 = (sp.rem(d*generator**i, source.as_expr(), x) for i in (1, 2))
    c1, c2 = (sp.Poly(n, x).nth(5) for n in (n1, n2))
    if c1 or c2:
        # Rank loss gives c1=c2=0 and is excluded by the branch-value discriminant.
        value = c2*u-c1*u*u
        n, equations, nonzero = map_conditions(source, pairing, s, c2*n1-c1*n2, value, degree)
        for root in bounded_intersection(equations, variables, nonzero):
            yield s, n, value, dict(zip(variables, root, strict=True))


def recover(source):
    """All rational optimal degree-four classes, with their branch pairings."""
    source, center = sextic_chart(source)
    output = []
    for pairing in paired_algebras(source):
        for degree in range(3):
            for s, n, b, subs in pairing_candidates(source, pairing, degree):
                half_poles,numerator = (normalize(expr.subs(subs), x) for expr in (s, n))
                g = cubic_value(pairing,b.subs(subs))
                cover = build_from_cubic(source,numerator,half_poles,g)
                if cover is None:
                    raise ArithmeticError("A saturated point failed the exact cover identity.")
                output.append((cover, pairing))
    return source,center,output


def class_value(source, cover, pairing):
    powers = [(pairing.generator**i).rem(source) for i in range(3)]
    abscissa = (cover.numerator*(cover.half_poles**2).invert(source)).rem(source)
    basis = sp.Matrix.hstack(*(vector(p,6) for p in powers))
    c,_ = basis.gauss_jordan_solve(vector(abscissa,6))
    return sp.Poly.from_list(list(reversed(c)),u,domain=sp.QQ)


def recover_kernel(source, first, second, pairing):
    v,w = (class_value(source, cover, pairing) for cover in (first, second))
    basis = sp.Matrix.hstack(*(vector((v**i).rem(pairing.modulus), 3) for i in range(3)))
    c, _ = basis.gauss_jordan_solve(vector(w, 3))
    at_two = sp.Poly.from_list(list(reversed(c)),t,domain=sp.QQ)
    assert cubic_value(pairing,v.as_expr()) == first.cubic
    assert cubic_value(pairing,w.as_expr()) == second.cubic
    return at_two, recover_order_four(source, first, second)


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("sextic", type=json.loads, help="Ascending coefficients of a rational quintic or sextic.")
    parser.add_argument("--graph", action="store_true", help="Also verify the graph for a complementary pair.")
    args = parser.parse_args()
    source, center, maps = recover(sp.Poly.from_list([sp.Rational(c) for c in args.sextic[::-1]], x, domain=sp.QQ))
    result = {
        "reciprocal_chart_center": center,
        "maps": [
            {name: coefficients(getattr(cover, name)) for name in (
                "numerator", "half_poles", "lift", "cubic"
            )} | {"pairing_generator": coefficients(pairing.generator),
                   "pairing_modulus": coefficients(pairing.modulus)}
            for cover, pairing in maps
        ],
    }
    if args.graph:
        if len(maps) != 2 or maps[0][1] != maps[1][1]:
            parser.error("--graph requires two covers with the same branch pairing.")
        first, second = (cover for cover, _ in maps)
        two, four = recover_kernel(source, first, second, maps[0][1])
        result["order_two_x_map"] = coefficients(two)
        result["order_four"] = [
            dict(zip(("factor", "x_map", "y_multiplier"), map(coefficients, component)))
            for component in four
        ]
    print(json.dumps(result, indent=2))
