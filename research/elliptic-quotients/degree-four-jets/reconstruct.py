"""Recover the Weierstrass x-coordinate from an unramified chart and four parameters."""
import argparse
from fractions import Fraction as Q
import json
from pathlib import Path

import sympy as sp


T, S = sp.symbols("t s")
PRECISION = 57
CUBICS = [(i, j) for j in range(4) for i in range(4-j)]


def coefficient(a, b, n):
    return sum(a[k]*b[n-k] for k in range(max(0, n+1-len(b)), min(n+1, len(a))))


def product(a, b):
    return [coefficient(a, b, n) for n in range(PRECISION)]


def branch_powers(poly):
    """Taylor expansions of 1,s,...,s^4 on the specified smooth local branch."""
    terms = {(i,j): Q(c) for (i,j),c in sp.Poly(poly,T,S).terms()}
    slope = terms.get((0,1), 0)
    if terms.get((0,0), 0) or not slope or any(i+j > 4 for i,j in terms):
        raise ValueError("Expected a plane quartic chart through the origin with F_s(0) nonzero")
    powers = [[0]*PRECISION for _ in range(5)]
    powers[0][0] = 1
    for n in range(1, PRECISION):
        for j in range(2, 5):
            powers[j][n] = coefficient(powers[j-1], powers[1], n)
        remainder = sum(c*powers[j][n-i] for (i,j),c in terms.items() if i<=n)
        powers[1][n] = -remainder/slope
    return terms, powers


def differential(terms, powers, u, v):
    slope = terms[0,1]
    denominator = [sum(j*c*powers[j-1][n-i] for (i,j),c in terms.items() if j and i<=n)
                   for n in range(PRECISION)]
    values = []
    for n in range(PRECISION):
        numerator = slope*(n==0) + u*(n==1) + v*powers[1][n]
        values.append((numerator-coefficient(denominator, values, n))/slope)
    return values


def weierstrass_jet(h, a, b):
    h2 = product(h, h)
    r, square, cubic_rhs, deriv = [Q(1)], [Q(1)], [Q(1)], [Q(-2)]
    for n in range(1, PRECISION):
        square.append(coefficient(r, r, n))
        value = coefficient(square, r, n)
        if n >= 4:
            value += a*r[n-4]
        if n == 6:
            value += b
        cubic_rhs.append(value)
        error = coefficient(deriv, deriv, n) - 4*coefficient(h2, cubic_rhs, n)
        new = error/(4*(n+1))
        r.append(new)
        square[n] += 2*new
        cubic_rhs[n] += 3*new
        deriv.append((n-2)*new)
    return r


def check_source(chart):
    z = sp.Symbol("z")
    form = sp.Poly(chart, T, S).homogenize(z).as_expr()
    if sp.Poly(form, T, S, z).total_degree() != 4:
        raise ValueError("Expected an equation of total degree four")
    coords = (T, S, z)
    for coord in coords:
        others = [x for x in coords if x != coord]
        partials = [sp.diff(form, x).subs(coord, 1) for x in coords]
        if list(sp.groebner(partials, *others, domain=sp.QQ)) != [1]:
            raise ValueError("Singular projective source")


def reconstruct(chart, u, v, a, b):
    check_source(chart)
    u, v, a, b = map(Q, (u, v, a, b))
    if not 4*a**3+27*b**2:
        raise ValueError("Singular elliptic target")
    terms, powers = branch_powers(chart)
    r = weierstrass_jet(differential(terms, powers, u, v), a, b)
    basis = [[0]*i + powers[j][:PRECISION-i] for i,j in CUBICS]
    columns = [[0,0]+col[:-2] for col in basis] + [[-v for v in product(r,col)] for col in basis]
    matrix = sp.Matrix(list(zip(*columns))).to_DM().to_field()
    kernel = matrix.nullspace()
    nullity = kernel.shape[0]
    result = {"rank": 20-nullity, "nullity": nullity, "degree_four": nullity in (2,3)}
    if not nullity:
        return result
    coeffs = list(kernel.to_Matrix().row(0))
    monomials = [T**i*S**j for i,j in CUBICS]
    n, d = [sum(k*m for k,m in zip(block,monomials)) for block in (coeffs[:10],coeffs[10:])]
    x = sp.cancel(n/d)
    n, d = sp.fraction(x)
    # Tangent derivation is F_s*d/dt - F_t*d/ds. Verify the curve equation
    # independently of the order test on the reconstructed rational function.
    tangent = lambda g: sp.diff(chart,S)*sp.diff(g,T)-sp.diff(chart,T)*sp.diff(g,S)
    wronskian = d*tangent(n)-n*tangent(d)
    ell = terms[0,1]+u*T+v*S
    residual = wronskian**2 - 4*ell**2*d*(n**3+a*n*d**2+b*d**3)
    if sp.rem(residual, chart, S).expand() != 0:
        raise ArithmeticError("Reconstructed functions failed the global Weierstrass identity")
    result["x"] = str(x)
    result["weierstrass_identity"] = True
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON chart and normalized coefficients")
    data = json.loads(parser.parse_args().input.read_text())
    chart = sum(Q(c)*T**i*S**j for i,j,c in data["chart"])
    print(json.dumps(reconstruct(chart, *data["parameters"]), indent=2))


if __name__ == "__main__":
    main()
