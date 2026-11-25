"""Enumerate primitive isotropic classes by exact positive definite degree slices."""

from collections import Counter
from fractions import Fraction
from math import gcd, isqrt

import sympy as sp


def degree_basis(row):
    """Unimodular columns taking the integer degree row to (gcd, 0, ...)."""
    rank = len(row)
    columns = [[int(i == j) for i in range(rank)] for j in range(rank)]
    pivot = int(row[0])
    for j in range(1, rank):
        entry = int(row[j])
        if not entry:
            continue
        u, v, g = map(int, sp.gcdex(pivot, entry))
        first, last = columns[0], columns[j]
        columns[0] = [u*a+v*b for a, b in zip(first, last)]
        columns[j] = [-(entry//g)*a+(pivot//g)*b for a, b in zip(first, last)]
        pivot = g
    if pivot < 0:
        columns[0] = [-a for a in columns[0]]
        pivot = -pivot
    if pivot == 0:
        raise ValueError("The polarization pairing is zero.")
    return pivot, sp.Matrix.hstack(*(sp.Matrix(c) for c in columns))


def spectrum_from_lattice(gram, theta, bound):
    """Return the degree counts and all oriented primitive classes through bound."""
    G = sp.Matrix(gram)
    polarization = sp.Matrix(theta)
    if G.rows < 2 or G != G.T or (polarization.T*G*polarization)[0] != 2 or bound < 0:
        raise ValueError("Invalid lattice, principal polarization, or bound.")
    g, U = degree_basis(list(polarization.T*G))
    A = U.T*G*U
    H = -A[1:, 1:]
    if not H.is_positive_definite:
        raise ValueError("The polarization perpendicular is not negative definite.")
    L, D = H.LDLdecomposition()
    b = A[1:, 0]
    shift = H.inv()*b
    norm = Fraction(A[0, 0]+(b.T*shift)[0])
    k = H.rows
    diag = [Fraction(D[i, i]) for i in range(k)]
    lower = [[Fraction(L[i, j]) for j in range(k)] for i in range(k)]
    shift = list(map(Fraction, shift))
    transform = U.tolist()
    counts, classes = Counter(), []
    for degree in range(g, bound+1, g):
        h = degree//g
        center = [h*x for x in shift]
        vectors = [0]*k

        def visit(i, remaining):
            offset = center[i] - sum(lower[j][i]*(vectors[j]-center[j]) for j in range(i+1, k))
            a, b = offset.numerator, offset.denominator
            limit = remaining*b*b/diag[i]
            radius = isqrt(limit.numerator//limit.denominator)
            if i == 0:
                if radius*radius != limit:
                    return
                numerators = (a-radius, a+radius) if radius else (a,)
                choices = [numerator//b for numerator in numerators if numerator % b == 0]
            else:
                choices = range(-((radius-a)//b), (a+radius)//b+1)
            for value in choices:
                delta = Fraction(b*value-a, b)
                used = diag[i]*delta*delta
                if used > remaining:
                    continue
                vectors[i] = value
                if i:
                    visit(i-1, remaining-used)
                elif used == remaining and gcd(h, *vectors) == 1:
                    coordinates = [row[0]*h+sum(c*v for c, v in zip(row[1:], vectors))
                                   for row in transform]
                    counts[degree] += 1
                    classes.append((degree, list(map(int, coordinates))))

        visit(k-1, h*h*norm)
    return dict(sorted(counts.items())), classes
