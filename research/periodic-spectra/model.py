"""Exact cyclic Jacobi spectra; all polynomial coefficients are ascending integers.

The transfer/Chebyshev identity is classical.  This module validates the input,
discovers its actual coefficient period, and makes the reduction checkable.
"""
import argparse
import json
from pathlib import Path


class StructureError(ValueError):
    """Input is outside the supported nonzero-bond cyclic Jacobi class."""


def _integer(value):
    return type(value) is int


def expand_case(case):
    a, b, repeats = case["diagonal"], case["bonds"], case["repeats"]
    if (not isinstance(a, list) or not isinstance(b, list) or not a
            or len(a) != len(b) or not all(_integer(v) for v in a + b)
            or any(v == 0 for v in b) or not _integer(repeats) or repeats < 1):
        raise StructureError("Use equal nonempty integer arrays, nonzero bonds and positive repeats.")
    a, b = a * repeats, b * repeats
    if len(a) < 3:
        raise StructureError("A simple cyclic chain needs at least three sites.")
    if "defect" in case:
        defect = case["defect"]
        i, delta = defect["index"], defect["diagonal_delta"]
        if not _integer(i) or not 0 <= i < len(a) or not _integer(delta):
            raise StructureError("A diagonal defect needs an in-range integer index and delta.")
        a[i] += delta
    return a, b


def coefficient_period(a, b):
    """Least joint period in O(L) exact comparisons (prefix-function algorithm)."""
    if not a or len(a) != len(b):
        raise StructureError("Coefficient arrays must have the same positive length.")
    prefix = [0] * len(a)
    for i in range(1, len(a)):
        k = prefix[i - 1]
        while k and (a[i], b[i]) != (a[k], b[k]):
            k = prefix[k - 1]
        if (a[i], b[i]) == (a[k], b[k]):
            k += 1
        prefix[i] = k
    m = len(a) - prefix[-1]
    return m if len(a) % m == 0 else len(a)


def matrix_from_case(case):
    a, b = expand_case(case)
    size = len(a)
    matrix = [[0] * size for _ in range(size)]
    for i in range(size):
        matrix[i][i] = a[i]
        j = (i + 1) % size
        matrix[i][j] = matrix[j][i] = b[i]
    return matrix


def detect_period(matrix):
    """Validate the ENTIRE supplied matrix before inferring a coefficient period."""
    n = len(matrix)
    if n < 3 or any(len(row) != n for row in matrix):
        raise StructureError("Expected a square matrix of order at least three.")
    for i, row in enumerate(matrix):
        for j, value in enumerate(row):
            if not _integer(value) or value != matrix[j][i]:
                raise StructureError("Expected a symmetric integer matrix.")
            if j not in {i, (i - 1) % n, (i + 1) % n} and value != 0:
                raise StructureError("Off-cycle couplings are outside the supported class.")
    a = [matrix[i][i] for i in range(n)]
    b = [matrix[i][(i + 1) % n] for i in range(n)]
    if any(v == 0 for v in b):
        raise StructureError("Every cyclic bond must be nonzero.")
    m = coefficient_period(a, b)
    return a[:m], b[:m], n // m


def trim(p):
    while len(p) > 1 and p[-1] == 0:
        p.pop()
    return p


def add(p, q, scale=1):
    out = p[:] + [0] * max(0, len(q) - len(p))
    for i, value in enumerate(q):
        out[i] += scale * value
    return trim(out)


def multiply(p, q):
    out = [0] * (len(p) + len(q) - 1)
    for i, u in enumerate(p):
        if u:
            for j, v in enumerate(q):
                if v:
                    out[i + j] += u * v
    return trim(out)


def scale(p, c):
    return trim([v * c for v in p])


def discriminant(a, b):
    p, q, r, s = [1], [0], [0], [1]
    product = 1
    for i, diagonal in enumerate(a):
        oldp, oldq = p, q
        p = add(multiply([-diagonal, 1], p), r, -b[i - 1])
        q = add(multiply([-diagonal, 1], q), s, -b[i - 1])
        r, s = scale(oldp, b[i]), scale(oldq, b[i])
        product *= b[i]
    return add(p, s), product


def trace_power(d, product, repeats):
    previous, current = [2], d[:]
    for _ in range(2, repeats + 1):
        previous, current = current, add(multiply(d, current), previous, -product * product)
    return current


def certificate(case):
    a, b = expand_case(case)
    m = coefficient_period(a, b)
    repeats = len(a) // m
    d, product = discriminant(a[:m], b[:m])
    polynomial = trace_power(d, product, repeats)
    polynomial[0] -= 2 * product ** repeats
    return {"sites": len(a), "period": m, "repeats": repeats,
            "discriminant": d, "bond_product": product,
            "coefficients": polynomial}


def solve(case):
    return certificate(case)["coefficients"]


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path)
    args = parser.parse_args()
    print(json.dumps(certificate(json.loads(args.case.read_text())), indent=2))
