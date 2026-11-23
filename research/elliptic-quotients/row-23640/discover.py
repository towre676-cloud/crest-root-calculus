"""Find rational involutions of a smooth plane quartic by enumerating centers.

Input is [[i,j,k,coefficient], ...] for F(x,y,z), or {"terms": [...]}.
Rational coefficients can be integers or strings such as "3/5". Run
    python discover.py quartic.json --height 3
The height is the largest absolute coordinate in a primitive integer center.
The program proves each reported symmetry by polynomial substitution. An empty
list exhausts only the given height in the supplied coordinate system.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import product
import json
from math import gcd, lcm
from pathlib import Path
import sys


def parse_quartic(data):
    terms = data.get("terms") if isinstance(data, dict) else data
    if not isinstance(terms, list):
        raise ValueError("Expected a list of quartic monomials.")
    polynomial = defaultdict(Fraction)
    for term in terms:
        if not isinstance(term, list) or len(term) != 4:
            raise ValueError("Each monomial needs three exponents and a coefficient.")
        exponent = tuple(term[:3])
        if any(type(a) is not int or a < 0 for a in exponent) or sum(exponent) != 4:
            raise ValueError("Exponents must be nonnegative integers with sum four.")
        if type(term[3]) not in (int, str):
            raise ValueError("Use an integer or a rational string for each coefficient.")
        polynomial[exponent] += Fraction(term[3])
    denominator = lcm(*(c.denominator for c in polynomial.values()))
    coefficients = {m: int(c * denominator) for m, c in polynomial.items() if c}
    if not coefficients:
        raise ValueError("The quartic polynomial is zero.")
    return coefficients


def value(polynomial, point):
    x, y, z = point
    return sum(c * x**i * y**j * z**k for (i, j, k), c in polynomial.items())


def derivatives(polynomial):
    result = [{}, {}, {}]
    for monomial, c in polynomial.items():
        for j, exponent in enumerate(monomial):
            if exponent:
                index = list(monomial)
                index[j] -= 1
                result[j][tuple(index)] = exponent * c
    return result


def primitive(vector):
    divisor = gcd(*vector)
    if next(v for v in vector if v) < 0:
        divisor = -divisor
    return tuple(v // divisor for v in vector)


def vanishes_on_axis(gradient, center, axis):
    # The cubic D_P F must vanish on L=0. These four evaluations
    # determine its restriction to the binary line exactly.
    j = next(k for k, a in enumerate(axis) if a)
    basis = []
    for k in range(3):
        if k != j:
            vec = [0, 0, 0]
            vec[k], vec[j] = axis[j], -axis[k]
            basis.append(vec)
    u, v = basis
    for t in range(4):
        point = [a + t * b for a, b in zip(u, v)]
        if sum(c * value(g, point) for c, g in zip(center, gradient)):
            return False
    return True


def substitute(polynomial, matrix):
    """Expand F(M*(x,y,z)) over the integers, returning nonzero coefficients."""
    result = defaultdict(int)
    for exponents, c in polynomial.items():
        expanded = {(0, 0, 0): c}
        for row, power in zip(matrix, exponents):
            for _ in range(power):
                nxt = defaultdict(int)
                for monomial, coefficient in expanded.items():
                    for k, entry in enumerate(row):
                        if entry:
                            index = list(monomial)
                            index[k] += 1
                            nxt[tuple(index)] += coefficient * entry
                expanded = nxt
        for monomial, coefficient in expanded.items():
            result[monomial] += coefficient
    return {m: c for m, c in result.items() if c}


def discover(data, height):
    if type(height) is not int or height < 1:
        raise ValueError("The center height must be a positive integer.")
    quartic = parse_quartic(data)
    gradient = derivatives(quartic)
    symmetries = []
    for center in product(range(-height, height + 1), repeat=3):
        if gcd(*center) != 1 or next(c for c in center if c) < 0:
            continue
        f = value(quartic, center)
        if f == 0:
            continue
        # Euler gives dot(grad F(P), P)=4F(P). In coordinates with P as
        # center and a fixed plane as axis, the cubic-in-center coefficient
        # must vanish on that plane, forcing L to be the polar at P.
        polar = [value(g, center) for g in gradient]
        axis = primitive(polar)
        if not vanishes_on_axis(gradient, center, axis):
            continue
        entries = primitive([2*f*(i == j) - center[i]*polar[j]
                             for i in range(3) for j in range(3)])
        matrix = [entries[k:k + 3] for k in (0, 3, 6)]
        # Eigenvalues are a,a,-a, so a equals the trace in this normalization.
        trace = sum(matrix[i][i] for i in range(3))
        if substitute(quartic, matrix) == {m: c * trace**4 for m, c in quartic.items()}:
            symmetries.append({"center": center, "axis": axis, "matrix": matrix})
    return symmetries


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", help="Quartic monomial JSON file, or - for stdin.")
    parser.add_argument("--height", type=int, default=1)
    args = parser.parse_args()
    try:
        raw = sys.stdin.read() if args.input == "-" else Path(args.input).read_text(encoding="utf-8")
        result = discover(json.loads(raw), args.height)
    except (ValueError, ZeroDivisionError, OSError) as error:
        parser.exit(2, f"Input error: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
