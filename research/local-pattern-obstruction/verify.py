"""Verify the Klein-four obstruction to reconstruction from point counts."""
import json
from fractions import Fraction
from itertools import product

import sympy as sp


X = sp.Symbol("x")
S_COMPONENTS = [X**2 - 13, X**2 - 17, X**2 - 221]
T_COMPONENTS = [X, X, X**4 - 60*X**2 + 16]


def fixed_points(component, prime, extension):
    """Count distinct roots over F_(prime**extension) by Frobenius gcd."""
    polynomial = sp.Poly(component, X, modulus=prime)
    power = sp.Poly(1, X, modulus=prime)
    base = sp.Poly(X, X, modulus=prime)
    exponent = prime**extension
    while exponent:
        if exponent & 1:
            power = (power*base).rem(polynomial)
        base = (base*base).rem(polynomial)
        exponent >>= 1
    return int(polynomial.gcd(power - sp.Poly(X, X, modulus=prime)).degree())


def rational_points(components):
    return sum(int(f.degree() == 1) for component in components
               for f, multiplicity in sp.Poly(component, X, domain=sp.QQ).factor_list()[1])


def verify():
    # Traces of every group element and all powers follow from these four rows.
    characters = []
    for first, second in product((-1, 1), repeat=2):
        left = 3 + first + second + first*second
        right = 2 + (1 + first)*(1 + second)
        if left != right:
            raise ArithmeticError("The two permutation characters differ")
        characters.append({"chi_13": first, "chi_17": second,
                           "S_fixed": left, "T_fixed": right})

    # The minimal polynomial of sqrt(13)+sqrt(17), checked by elimination.
    A = sp.Symbol("a")
    quartic = sp.resultant(A**2 - 13, (X - A)**2 - 17, A)
    if sp.expand(quartic) != T_COMPONENTS[2]:
        raise ArithmeticError("Wrong biquadratic generator polynomial")
    for polynomial in [*S_COMPONENTS, T_COMPONENTS[2]]:
        if not sp.Poly(polynomial, X, domain=sp.QQ).is_irreducible:
            raise ArithmeticError("A claimed field polynomial is reducible")
    if (rational_points(S_COMPONENTS), rational_points(T_COMPONENTS)) != (0, 2):
        raise ArithmeticError("Wrong rational-point count")

    # Each component is a separate scheme, even when coordinate values overlap.
    rows = []
    for prime in sp.primerange(3, 101):
        if prime in (13, 17):
            continue
        first = int(sp.legendre_symbol(13, prime))
        second = int(sp.legendre_symbol(17, prime))
        counts = []
        for extension in range(1, 7):
            expected = 3 + first**extension + second**extension + (first*second)**extension
            left = sum(fixed_points(f, prime, extension) for f in S_COMPONENTS)
            right = sum(fixed_points(f, prime, extension) for f in T_COMPONENTS)
            if (left, right) != (expected, expected):
                raise ArithmeticError(f"Frobenius check failed at {prime}^{extension}")
            counts.append(left)
        rows.append({"prime": prime, "extension_degrees": [1, 2, 3, 4, 5, 6],
                     "both_counts": counts})

    # At every completion, one of 13,17,221 is a square. The general odd-prime
    # argument is in the proof. These are the ramified-prime checks it requires.
    if 17 % 8 != 1 or pow(17 % 13, 6, 13) != 1 or pow(13, 8, 17) != 1:
        raise ArithmeticError("The exceptional local square tests failed")

    mean = sum(Fraction(row["S_fixed"], 4) for row in characters)
    return {
        "rational_points": {"S": 0, "T": 2},
        "geometric_points": {"S": 6, "T": 6},
        "rational_component_degrees": {"S": [2, 2, 2], "T": [1, 1, 4]},
        "common_character": characters,
        "mean_fixed_point_count": str(mean),
        "good_prime_model": "all primes outside 2, 13, 17",
        "biquadratic_minimal_polynomial": str(quartic),
        "local_existence_at_every_completion": {
            "real": "13 > 0",
            "2": "17 = 1 modulo 8",
            "13": "17 = 4 modulo 13",
            "17": "13 = 8^2 modulo 17",
            "other_odd_primes": "at least one of a, b, a*b equals +1",
        },
        "independent_frobenius_checks": rows,
    }


if __name__ == "__main__":
    print(json.dumps(verify(), indent=2))
