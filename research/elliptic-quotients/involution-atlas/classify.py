"""Compute all geometric bielliptic involutions of a smooth rational plane quartic."""
import argparse
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.rings import ring

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "row-23640"))
from discover import parse_quartic


def derivatives(form):
    result = []
    for k in range(3):
        derivative = {}
        for exponents, c in form.items():
            if exponents[k]:
                key = list(exponents)
                key[k] -= 1
                derivative[tuple(key)] = exponents[k] * c
        result.append(derivative)
    return result


def evaluate(form, vector, one):
    value = 0 * one
    for exponents, coefficient in form.items():
        term = coefficient * one
        for x, degree in zip(vector, exponents):
            if degree:
                term *= x**degree
        value += term
    return value


def simplex(degree):
    return [(i, j, degree-i-j) for i in range(degree+1) for j in range(degree-i+1)]


def center_equations(point, gradient, one):
    """Ten evaluations determine a homogeneous cubic in the test vector."""
    polar = [evaluate(g, point, one) for g in gradient]
    equations = []
    for h in simplex(3):
        on_axis = [polar[(j+1)%3]*h[(j+2)%3] - polar[(j+2)%3]*h[(j+1)%3]
                   for j in range(3)]
        equations.append(sum(a*evaluate(g, on_axis, one) for a, g in zip(point, gradient)))
    return equations


def univariate_relation(basis, var):
    polys = [sp.Poly(f, var, domain=sp.QQ) for f in basis if f.free_symbols <= {var}]
    if not polys:
        raise ArithmeticError(f"Missing finite-coordinate relation for {var}")
    return min(polys, key=lambda f: f.degree()).monic()


def plane_orbits(equations, a, b, parameter):
    """Rational univariate representations of the distinct affine solutions."""
    basis = sp.groebner(equations, a, b, order="grevlex", domain=sp.QQ)
    if list(basis) == [1]:
        return []
    if not basis.is_zero_dimensional:
        raise ArithmeticError("Unexpected positive-dimensional center incidence")
    min_b = univariate_relation(basis.fglm("lex"), b)
    min_a = univariate_relation(sp.groebner(list(basis), b, a, order="lex", domain=sp.QQ), a)
    # Adding both squarefree relations makes the algebra a quotient of a finite
    # etale Q-algebra, hence reduced, while preserving every geometric point.
    rad = sp.groebner([*basis, min_a.sqf_part().as_expr(), min_b.sqf_part().as_expr()],
                      a, b, order="grevlex", domain=sp.QQ)
    bound = min_a.sqf_part().degree() * min_b.sqf_part().degree()
    for k in range(bound*(bound-1)//2 + 1):
        sheared = [f.subs(a, parameter - k*b) for f in rad]
        lex = sp.groebner(sheared, b, parameter, order="lex", domain=sp.QQ)
        linear = next((f for f in lex if sp.degree(f, b) == 1
                       and sp.diff(f, b) == 1), None)
        if linear is None:
            continue
        min_t = univariate_relation(lex, parameter)
        if not min_t.is_sqf:
            raise ArithmeticError("A reduced coordinate relation has a repeated root")
        orbits = []
        for factor, _ in min_t.factor_list()[1]:
            bval = sp.Poly(-linear.subs(b, 0), parameter, domain=sp.QQ).rem(factor)
            aval = sp.Poly(parameter-k*bval.as_expr(), parameter, domain=sp.QQ).rem(factor)
            orbits.append((factor, [aval.as_expr(), bval.as_expr(), sp.Integer(1)]))
        return orbits
    raise ArithmeticError("No separating linear form within the pair-count bound")


def quotient_model(form, modulus, point, polar, f, reduce):
    """Give the double-cover target over the center field in a linear frame."""
    T = modulus.gens[0]
    one = sp.Integer(1)
    pivot = next(j for j, value in enumerate(polar) if value != 0)
    inverse = sp.invert(sp.Poly(polar[pivot], T, domain=sp.QQ), modulus).as_expr()
    base = []
    for j in range(3):
        if j != pivot:
            column = [sp.Integer(0)] * 3
            column[j] = one
            column[pivot] = reduce(-polar[j]*inverse)
            base.append(column)
    def fiber(u, v, t):
        return [reduce(u*a + v*b + t*c) for a, b, c in zip(*base, point)]

    # In these coordinates the quartic is f*t^4 + A*t^2 + B.
    # Interpolate binary forms on v=1, then check the homogeneous identity.
    def interpolate(values):
        n = len(values)
        matrix = sp.Matrix([[i**j for j in range(n)] for i in range(n)])
        return [reduce(value) for value in matrix.inv()*sp.Matrix(values)]

    Bvalues = [reduce(evaluate(form, fiber(u, 1, 0), one)) for u in range(5)]
    Avalues = [reduce(evaluate(form, fiber(u, 1, 1), one) - f - Bvalues[u])
               for u in range(3)]
    A, B = interpolate(Avalues), interpolate(Bvalues)

    def binary(coeffs, u, v):
        return sum(c*u**i*v**(len(coeffs)-1-i) for i, c in enumerate(coeffs))

    for u, v, t in simplex(4):
        original = evaluate(form, fiber(u, v, t), one)
        if reduce(original-f*t**4-binary(A, u, v)*t**2-binary(B, u, v)) != 0:
            raise ArithmeticError("Even quartic model failed verification")
    D = [reduce(sum(A[i]*A[k-i] for i in range(3) if 0 <= k-i < 3) - 4*f*B[k])
         for k in range(5)]
    e, d, c, b, a = D
    I = reduce(12*a*e-3*b*d+c*c)
    J = reduce(72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3)
    discr = reduce(4*I**3-J**2)
    if discr == 0:
        raise ArithmeticError("Singular genus-one quotient")
    denom_inv = sp.invert(sp.Poly(discr, T, domain=sp.QQ), modulus).as_expr()
    return {"frame_columns": [[str(v) for v in col] for col in [*base, point]],
            "t4_coefficient": str(f), "A_ascending": list(map(str, A)),
            "B_ascending": list(map(str, B)), "D_ascending": list(map(str, D)),
            "map": "[u:v:t] -> [u:v:2*f*t^2+A(u,v)]", "j": str(reduce(6912*I**3*denom_inv))}


def verify_orbit(form, gradient, modulus, point):
    var = modulus.gens[0]

    def reduce(value):
        return sp.Poly(value, var, domain=sp.QQ).rem(modulus).as_expr()

    f = reduce(evaluate(form, point, sp.Integer(1)))
    if f == 0:
        return None
    polar = [reduce(evaluate(g, point, sp.Integer(1))) for g in gradient]
    matrix = [[reduce(2*f*int(i == j) - point[i]*polar[j]) for j in range(3)]
              for i in range(3)]
    # A homogeneous degree-four identity is determined by these 15 points.
    scale = reduce((2*f)**4)
    for test in simplex(4):
        transformed = [reduce(sum(m*v for m, v in zip(row, test))) for row in matrix]
        diff = evaluate(form, transformed, sp.Integer(1)) - scale*evaluate(form, test, 1)
        if reduce(diff) != 0:
            raise ArithmeticError("Projective automorphism identity failed")
    return {"field_polynomial": str(modulus.as_expr()), "orbit_size": modulus.degree(),
            "center": [str(v) for v in point], "axis": [str(v) for v in polar],
            "quotient": quotient_model(form, modulus, point, polar, f, reduce)}


def classify(data):
    form = parse_quartic(data)
    gradient = derivatives(form)
    a, b, T = sp.symbols("a b T")
    x, y, z = sp.symbols("x y z")
    coords = (x, y, z)
    partials = [evaluate(g, coords, sp.Integer(1)) for g in gradient]
    for k in range(3):
        chart = {coords[k]: 1}
        others = [v for v in coords if v != coords[k]]
        if list(sp.groebner([g.subs(chart) for g in partials], *others, domain=sp.QQ)) != [1]:
            raise ValueError("The plane quartic is singular over the algebraic closure")
    raw = []
    for dimension, gens in [(2, (a, b)), (1, (a,))]:
        R, *indeterminates = ring(gens, sp.QQ)
        point = [*indeterminates, R.one, *([R.zero]*(2-dimension))]
        eqs = [e.as_expr() for e in center_equations(point, gradient, R.one) if e]
        if dimension == 2:
            raw.extend(plane_orbits(eqs, a, b, T))
        else:
            if not eqs:
                raise ArithmeticError("Unexpected line of solutions at infinity")
            gcd = sp.Poly(sp.gcd_list(eqs).subs(a, T), T, domain=sp.QQ)
            for factor, _ in gcd.factor_list()[1]:
                coordinate = sp.Poly(T, T, domain=sp.QQ).rem(factor).as_expr()
                raw.append((factor, [coordinate, sp.Integer(1), sp.Integer(0)]))
    fixed = [sp.Integer(1), sp.Integer(0), sp.Integer(0)]
    if all(q == 0 for q in center_equations(fixed, gradient, sp.Integer(1))):
        raw.append((sp.Poly(T, T, domain=sp.QQ), fixed))
    orbits = []
    hyperflexes = 0
    for modulus, point in raw:
        record = verify_orbit(form, gradient, modulus, point)
        if record is None:
            hyperflexes += modulus.degree()
        else:
            orbits.append(record)
    return {"geometric_involutions": sum(r["orbit_size"] for r in orbits),
            "rational_involutions": sum(r["orbit_size"] == 1 for r in orbits),
            "hyperflexes_removed": hyperflexes, "orbits": orbits}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="Quartic monomials as in the row-23640 example")
    args = parser.parse_args()
    try:
        print(json.dumps(classify(json.loads(args.input.read_text(encoding="utf-8"))), indent=2))
    except (ValueError, ZeroDivisionError, OSError) as error:
        parser.exit(2, f"Input error: {error}\n")


if __name__ == "__main__":
    main()
