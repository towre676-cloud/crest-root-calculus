"""Check the first two good fibres directly over finite extensions, using Newton identities."""
from itertools import product
from math import prod
import sympy as sp
from rank import CUBICS


def traces_at_degree(prime, degree, supports):
    variable = sp.Symbol("t")
    elements = list(product(range(prime), repeat=degree))
    zero = (0,) * degree
    modulus = next(c for c in elements if sp.Poly(variable**degree +
                   sum(v * variable**j for j, v in enumerate(c)), variable, modulus=prime).is_irreducible)

    def multiply(x, y):
        coefficients = [0] * (2 * degree - 1)
        for j, u in enumerate(x):
            for k, v in enumerate(y):
                coefficients[j+k] += u*v
        for i in range(len(coefficients)-1, degree-1, -1):
            leading = coefficients[i] % prime
            for j, u in enumerate(modulus):
                coefficients[i-degree+j] -= leading*u
        return tuple(c % prime for c in coefficients[:degree])

    squares = {multiply(x, x) for x in elements}
    sums = dict.fromkeys(supports, 0)
    for x in elements:
        x2 = multiply(x, x)
        x3 = multiply(x, x2)
        characters = []
        for A, B, C in CUBICS:
            value = tuple((x3[j] + A*x2[j] + B*x[j] + (C if j == 0 else 0)) % prime
                          for j in range(degree))
            characters.append(0 if value == zero else (1 if value in squares else -1))
        for support in supports:
            sums[support] += prod(characters[int(i)] for i in support)
    return {s: -total - (len(s) % 2 == 0) for s, total in sums.items()}


def verify(small_fibres):
    for prime in (3, 5):
        expected = small_fibres[str(prime)]
        traces = [traces_at_degree(prime, r, expected) for r in range(1, 6)]
        for support, order in expected.items():
            genus = (3*len(support)-1)//2
            coefficients = [1]
            for k in range(1, genus+1):
                value = -sum(coefficients[k-i]*traces[i-1][support] for i in range(1, k+1))
                q, rem = divmod(value, k)
                if rem:
                    raise ArithmeticError("Nonintegral Weil coefficient")
                coefficients.append(q)
            actual = sum(coefficients) + sum(prime**(genus-i)*coefficients[i] for i in range(genus))
            if actual != order:
                raise ArithmeticError("Direct point counts and the archived finite order differ")
