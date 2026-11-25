"""Independent coefficient identities for cyclic-cover cross-checks."""

import sympy as sp

from cyclic import X, basic_equations, discriminant_cut, exported_map
from finite import rational_points


def original(source, k, d, m, e):
    horizontal, vertical, outer = (sp.symbols(f"{prefix}0:{length}")
                                    for prefix, length in (("a", d-1), ("b", e), ("c", m)))
    variables = horizontal + vertical + outer
    domain = sp.QQ.poly_ring(*variables)
    U = sp.Poly.from_list([1, *horizontal, 0], X, domain=domain)
    V = sp.Poly.from_list([1, *vertical], X, domain=domain)
    G = sp.Poly.from_list([source.LC(), *outer], X, domain=domain)
    variables, equations = discriminant_cut(basic_equations((G.compose(U)-source.set_domain(domain)*V**k).coeffs()),
                                             variables, G)
    return variables, equations, U, V


def raw_solve(source, k, d, m, e, basis_solver):
    variables, equations, U, V = original(source, k, d, m, e)
    certificate = basis_solver(equations, variables)
    maps = [exported_map(source, k, m, U, V, dict(zip(variables, map(sp.Rational, point), strict=True)))
            for point in certificate["rational_points"]]
    return {"geometric_points": certificate["geometric_points"],
            "maps": sorted(maps, key=lambda row: (row["inner"], row["outer"], row["vertical"]))}


def triangular(source, k, d, m, e):
    """Eliminate U and G from the original identity by constant linear pivots."""
    a = sp.symbols(f"a1:{d}")
    b = sp.symbols(f"b1:{e+1}")
    c = sp.symbols(f"c0:{m}")
    U = X**d + sum(v*X**(d-i) for i, v in enumerate(a, 1))
    V = X**e + sum(v*X**(e-i) for i, v in enumerate(b, 1))
    Y = sp.Symbol("Y")
    G = source.LC()*Y**m + sum(v*Y**i for i, v in enumerate(c))
    equation = sp.Poly(G.subs(Y, U) - source.as_expr()*V**k, X)
    pivots = [(d*m-i, var) for i, var in enumerate(a, 1)] + [(d*i, c[i]) for i in range(m-1, -1, -1)]
    values = {}
    for power, var in pivots:
        occurrence = sp.Poly(equation.nth(power), var)
        if occurrence.degree() != 1 or not occurrence.LC().is_Rational:
            raise ArithmeticError("A leading coefficient lost its constant linear pivot.")
        value = -occurrence.nth(0)/occurrence.LC()
        values[var] = value
        equation = sp.Poly(equation.as_expr().subs(var, value), X)
    outer = sp.Poly(G.subs(values), Y)
    variables, equations = discriminant_cut(basic_equations(equation.coeffs()), b, outer)
    if variables:
        solution = rational_points(equations, variables)
        points, count = solution["rational_points"], solution["geometric_points"]
    else:
        points = [] if equations else [[]]
        count = len(points)
    maps = []
    for point in points:
        assignment = dict(zip(variables, map(sp.Rational, point), strict=True))
        maps.append(exported_map(source, k, m, sp.Poly(U.subs(values), X), sp.Poly(V, X), assignment))
    return {"geometric_points": count, "maps": sorted(maps, key=lambda row: (row["inner"], row["outer"], row["vertical"]))}
