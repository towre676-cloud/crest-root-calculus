"""Complete constant-vertical-multiplier maps by exact polynomial decomposition."""
import sympy as sp

from finite import coefficients, rational_points, verify
from maps import X, curve


def inner_and_outer(source, outer_degree):
    n = source.degree()
    if n % outer_degree:
        raise ValueError("The source degree is not a multiple of the target degree.")
    d = n // outer_degree
    monic = source.monic()
    inner = X**d
    for i in range(d - 1, 0, -1):
        power = (outer_degree - 1) * d + i
        missing = monic.nth(power) - sp.Poly(inner**outer_degree, X).nth(power)
        inner += missing * X**i / outer_degree
    remainder = monic.as_expr()
    outer = 0
    for i in range(outer_degree, -1, -1):
        value = sp.Poly(remainder, X).nth(i * d)
        outer += value * X**i
        remainder = sp.expand(remainder - value * inner**i)
    return sp.Poly(inner, X), sp.Poly(outer, X), sp.Poly(remainder, X)


def reduced_incidence(source, target):
    inner, outer, remainder = inner_and_outer(source, target.degree())
    if remainder:
        return inner, outer, remainder, [], [], None
    a, s, inverse = sp.symbols("a s inverse")
    m = target.degree()
    target_monic = target.monic()
    b = (a * outer.nth(m - 1) - target_monic.nth(m - 1)) / m
    difference = sp.Poly(target_monic.as_expr().subs(X, a * X + b)
                         - a**m * outer.as_expr(), X)
    constraints = [sp.Poly(q, a, domain=sp.QQ) for q in difference.all_coeffs() if q != 0]
    eliminant = constraints[0]
    for constraint in constraints[1:]:
        eliminant = sp.gcd(eliminant, constraint)
    # a is nonzero, so factors of a describe excluded constant maps.
    while eliminant.nth(0) == 0:
        eliminant = eliminant.exquo(sp.Poly(a, a))
    equations = [eliminant.monic().as_expr(),
                 source.LC() * s**2 - target.LC() * a**m,
                 a * inverse - 1]
    return inner, outer, remainder, equations, [a, s, inverse], b


def recover_constant(source_values, target_values, prime=None):
    source, target = curve(source_values), curve(target_values)
    if source.degree() % target.degree():
        raise ValueError("Constant-multiplier maps require divisible branch degrees.")
    d = source.degree() // target.degree()
    inner, outer, remainder, equations, variables, shift = reduced_incidence(source, target)
    result = {"source": coefficients(source), "target": coefficients(target),
              "degree": d, "inner": coefficients(inner), "outer": coefficients(outer),
              "remainder": coefficients(remainder), "geometric_maps": 0,
              "rational_maps": []}
    if remainder:
        return result
    certificate = rational_points(equations, variables, prime)
    maps = []
    for point in certificate["rational_points"]:
        assignment = dict(zip(variables, map(sp.Rational, point)))
        U = sp.Poly(assignment[variables[0]] * inner.as_expr() + shift.subs(assignment), X)
        V = assignment[variables[1]]
        if sp.expand(target.as_expr().subs(X, U.as_expr()) - source.as_expr() * V**2) != 0:
            raise ArithmeticError("A recovered cover failed the source and target identity.")
        maps.append({"u": coefficients(U), "v_over_y": [str(V)], "degree": d})
    return {**result, "certificate": certificate, "geometric_maps": certificate["geometric_points"],
            "rational_maps": maps}


def verify_constant(result):
    source, target = curve(result["source"]), curve(result["target"])
    if source.degree() % target.degree() or result["degree"] != source.degree() // target.degree():
        raise ArithmeticError("Incorrect constant-multiplier degree.")
    inner, outer, remainder, equations, variables, shift = reduced_incidence(source, target)
    if (coefficients(inner), coefficients(outer), coefficients(remainder)) != (result["inner"], result["outer"], result["remainder"]):
        raise ArithmeticError("Incorrect inner decomposition.")
    if remainder:
        if result["geometric_maps"] != 0 or result["rational_maps"] != []:
            raise ArithmeticError("An indecomposable source acquired a cover.")
        return
    certificate = result["certificate"]
    if certificate["equations"] != list(map(str, equations)) or certificate["variables"] != list(map(str, variables)):
        raise ArithmeticError("The certificate uses a different coefficient ideal.")
    verify(certificate)
    expected = []
    for point in certificate["rational_points"]:
        assignment = dict(zip(variables, map(sp.Rational, point)))
        U = sp.Poly(assignment[variables[0]] * inner.as_expr() + shift.subs(assignment), X)
        expected.append({"u": coefficients(U), "v_over_y": [str(assignment[variables[1]])],
                         "degree": result["degree"]})
    if expected != result["rational_maps"] or result["geometric_maps"] != certificate["geometric_points"]:
        raise ArithmeticError("Incorrect exported cover list.")
