"""Crest–Root row 36518. Run python verify_points.py [integer k].

Requires SymPy. The only external arithmetic theorem used in the argument
that E(Q) consists of O and (0, +/-4) is Fermat's theorem for exponent 3.
See proof.md for the rationality and completeness deductions.
"""

import argparse
from math import isqrt

import sympy as sp


def require(condition, description):
    if not condition:
        raise ArithmeticError(description)


def check_identities():
    u, v, t, k = sp.symbols("u v t k")
    a = 4*u**2 - 8*u*v + 2*v**2
    b = -2*u**4 + 2*u**3*v - 4*u**2*v**2 + 8*u*v**3 - 3*v**4
    r = t**2 + 5*k*u*v
    f = r**2 + a*r + b
    w = 2*r + a
    d = (2*u)**3*(u-v) + 16*(u-v)**4
    require(sp.expand(w**2 - d - 4*f) == 0, "quotient identity")
    require(f.subs({u: 0, v: 0, t: 1}) == 1, "projective chart coverage")
    require(sp.expand(f.subs({u: 0, v: 1}) - (t**2-1)*(t**2+3)) == 0,
            "two finite elliptic fibres")
    require(sp.expand(f.subs({u: 1, v: 1}) - (t**2+5*k-1)**2) == 0,
            "elliptic fibre at infinity")

    # Smooth reduction modulo 5 implies smoothness in characteristic zero.
    f0 = f.subs(k, 0)
    diff = sp.Poly(f-f0, u, v, t, k)
    require(all(c % 5 == 0 for c in diff.coeffs()), "common reduction")
    coordinates = (u, v, t)
    for q in coordinates:
        derivatives = [f0.diff(z).subs(q, 1) for z in coordinates]
        variables = [z for z in coordinates if z != q]
        basis = sp.groebner(derivatives, *variables, modulus=5)
        require(list(basis) == [1], "smoothness in chart " + str(q))

    # Rational maps used in the reduction to the Fermat cubic.
    x, X, Y = sp.symbols("x X Y")
    identity = (x**3+16)*(x**3-128)**2 - (x**3+64)**3 + 432*x**6
    require(sp.expand(identity) == 0, "map to Y^2 = X^3 - 432")
    require(sp.expand((36+Y)**3+(36-Y)**3-(6*X)**3
                      - 216*(Y**2-X**3+432)) == 0, "map to Fermat cubic")
    require(all((n**3 - 12*n**2 + 64) % 5 != 0 for n in range(5)),
            "no rational preimage of X = 12 with x nonzero")
    return f, coordinates, k


def rational_points(parameter):
    points = [(0, 1, -1), (0, 1, 1)]
    radicand = 1 - 5*parameter
    if radicand > 0:
        root = isqrt(radicand)
        if root*root == radicand:
            points.extend([(1, 1, -root), (1, 1, root)])
    return points


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("k", type=int, nargs="?", default=0)
    args = parser.parse_args()
    equation, variables, parameter_symbol = check_identities()
    result = rational_points(args.k)
    for point in result:
        substitutions = dict(zip(variables, point))
        substitutions[parameter_symbol] = args.k
        require(equation.subs(substitutions) == 0, "point fails substitution")
    print("Polynomial identities and smooth reduction verified.")
    print(f"k = {args.k}. Rational points in [u:v:t] coordinates")
    for u0, v0, t0 in result:
        print(f"[{u0}:{v0}:{t0}]")
