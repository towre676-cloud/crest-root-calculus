"""Exact enumeration of nonzero integer vectors in a positive definite ellipsoid."""
from fractions import Fraction
from math import isqrt


def ldl(form):
    n = len(form)
    lower = [[Fraction(i == j) for j in range(n)] for i in range(n)]
    diagonal = []
    for i in range(n):
        for j in range(i):
            cross = sum(lower[i][k] * lower[j][k] * diagonal[k] for k in range(j))
            lower[i][j] = (Fraction(form[i][j]) - cross) / diagonal[j]
        norm = Fraction(form[i][i]) - sum(lower[i][j] ** 2 * diagonal[j] for j in range(i))
        if norm <= 0:
            raise ArithmeticError("The exact form is not positive definite.")
        diagonal.append(norm)
    return lower, diagonal


def short_vectors(form, bound):
    """One representative per sign, including vectors exactly on the boundary."""
    lower, diagonal = ldl(form)
    vector = [0] * len(form)

    def visit(i, remaining, oriented):
        if i < 0:
            if oriented:
                yield tuple(vector)
            return
        offset = sum((lower[j][i] * vector[j] for j in range(i + 1, len(form))), Fraction())
        a, b = offset.numerator, offset.denominator
        limit = remaining * b * b / diagonal[i]
        radius = isqrt(limit.numerator // limit.denominator)
        least, greatest = -((radius + a) // b), (radius - a) // b
        if not oriented:
            least = max(0, least)
        for value in range(least, greatest + 1):
            vector[i] = value
            cost = diagonal[i] * Fraction(b * value + a, b) ** 2
            yield from visit(i - 1, remaining - cost, oriented or value != 0)

    yield from visit(len(form) - 1, Fraction(bound), False)
