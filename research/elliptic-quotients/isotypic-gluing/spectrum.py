"""Exact optimal elliptic degree spectra for two genus-two families."""

import argparse
from collections import Counter
from fractions import Fraction
from functools import reduce
from itertools import product
from math import gcd, isqrt


def kernel_size(vector, modulus, kernel):
    """Count E[modulus] mapping into the supplied kernel in E^r."""
    return sum(
        tuple(((v*x) % modulus, (v*y) % modulus) for v in vector) in kernel
        for x, y in product(range(modulus), repeat=2)
    )


def residue_weights(rank, modulus, kernel):
    """Compute the finite factor in the degree distribution."""
    return {
        vector: kernel_size(vector, modulus, kernel)
        for vector in product(range(modulus), repeat=rank)
        if reduce(gcd, vector, modulus) == 1
    }


def diagonal_kernel():
    """Graph of the transposition of the two coordinate vectors on E[2]."""
    return {
        ((x, y), (y, x))
        for x, y in product(range(2), repeat=2)
    }


def cyclic_kernel():
    """Graph of a three-cycle on the nonzero points of E[2]."""
    return {
        ((x, y), (y, (x+y) % 2))
        for x, y in product(range(2), repeat=2)
    }


def family_degree(a, b):
    if gcd(a, b) != 1:
        raise ValueError("Use a primitive integer direction.")
    return (a*a+b*b) if a % 2 and b % 2 else 2*(a*a+b*b)


def companion_degree(a, b):
    return 2*(a*a+b*b)


def map_recipe(a, b):
    """Express the optimal map as a sum of centered degree-two maps."""
    degree = family_degree(a, b)
    if a % 2 and b % 2:
        return degree, ("g_plus", "g_minus"), ((a+b)//2, (a-b)//2)
    return degree, ("f_plus", "f_minus"), (a, b)


def enumerate_spectrum(bound, degree_of=family_degree):
    """Count primitive directions modulo overall sign."""
    radius = isqrt(bound)
    counts = Counter()
    for a in range(1, radius+1):
        max_b = isqrt(bound-a*a)
        for b in range(-max_b, max_b+1):
            if gcd(a, b) == 1:
                degree = degree_of(a, b)
                if degree <= bound:
                    counts[degree] += 1
    if bound >= 2:
        counts[2] += 1  # (0,1)
    return counts


def factor_formula(degree, prime_divisors):
    """Closed class count given the distinct prime divisors of degree/2."""
    if degree <= 0 or degree % 2:
        return 0
    if any(prime % 4 != 1 for prime in prime_divisors):
        return 0
    return 2**(len(prime_divisors)+2)


def describe_parameter(t):
    """The normalized cubic quotient has rational Legendre parameter 1+t^2."""
    t = Fraction(t)
    if t <= 1:
        raise ValueError("Choose a rational t greater than one.")
    lam = 1+t*t
    j = 256*(1-lam+lam*lam)**3/(lam*lam*(1-lam)**2)
    return t, lam, j


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter", default="2")
    parser.add_argument("--bound", type=int, default=100)
    parser.add_argument("--companion", action="store_true", help="Use the three-cycle gluing")
    parser.add_argument("--direction", type=int, nargs=2, metavar=("A", "B"))
    opts = parser.parse_args()
    if opts.bound < 0:
        parser.error("Choose a nonnegative degree bound.")
    try:
        t, lam, j = describe_parameter(opts.parameter)
    except (ValueError, ZeroDivisionError) as error:
        parser.error(str(error))
    if opts.direction is not None:
        try:
            if opts.companion:
                if gcd(*opts.direction) != 1:
                    raise ValueError("Use a primitive integer direction.")
                degree = companion_degree(*opts.direction)
                names = ("first companion map", "second companion map")
                weights = opts.direction
            else:
                degree, names, weights = map_recipe(*opts.direction)
        except ValueError as error:
            parser.error(str(error))
        print(f"Optimal degree {degree}: [{weights[0]}]{names[0]} + [{weights[1]}]{names[1]}")
        return
    print(f"t = {t}, Legendre parameter = {lam}, j = {j}")
    counts = enumerate_spectrum(opts.bound, companion_degree if opts.companion else family_degree)
    print("degree  geometric classes")
    for degree, count in sorted(counts.items()):
        print(f"{degree:6}  {count:17}")
    print(f"{sum(counts.values())} geometric classes through degree {opts.bound}")


if __name__ == "__main__":
    main()
