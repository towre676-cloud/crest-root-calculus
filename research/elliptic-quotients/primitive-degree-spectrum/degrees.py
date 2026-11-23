"""Exact primitive degree counts in a congruence family of smooth plane quartics."""
import argparse
import json
from math import gcd, isqrt, prod

import sympy as sp


SIEVE_PRIMES = (7, 11, 13, 17, 19)
RATIONAL_MODULUS = 15*prod(SIEVE_PRIMES)


def primitive_count(norm):
    factors = sp.factorint(norm)
    has_inert_prime = any(p % 4 == 3 for p in factors)
    directions = 0 if has_inert_prime or norm % 4 == 0 else 2**(1+sum(p % 4 == 1 for p in factors))
    return directions + int(norm == 1)


def spectrum(a, degree, list_directions=False):
    if type(a) is not int or a % 15 != 4:
        raise ValueError("This classification requires an integer a congruent to 4 modulo 15")
    if type(degree) is not int or degree < 1:
        raise ValueError("Degree must be a positive integer")
    result = {"a": a, "degree": degree}
    count = 0
    allowed = degree % 2 == 0
    n = degree//2
    norm_possible = False
    if allowed:
        factors = sp.factorint(n)
        inert = [e for p, e in factors.items() if p % 4 == 3]
        norm_possible = all(e % 2 == 0 for e in inert)
        count = primitive_count(n)
    geometric_count = sum(primitive_count(int(k)) * int(sp.divisor_sigma(n//k))
                          for k in sp.divisors(n)) if allowed else 0
    result.update(geometric_map=allowed, all_geometric_classes=geometric_count,
                  primitive_classes_over_Q_and_Qbar=count)
    if a % RATIONAL_MODULUS == 4:
        result["rational_map"] = allowed and norm_possible
    elif not allowed or norm_possible:
        result["rational_map"] = allowed
    else:
        result["rational_map"] = None
    if list_directions:
        # Quotient by simultaneous sign. The remaining labels are distinct.
        directions = []
        if count:
            for m in range(isqrt(n)+1):
                delta = n-m*m
                j = isqrt(delta)
                if j*j == delta and gcd(m, j) == 1:
                    directions.append([m, j])
                    if m and j:
                        directions.append([m, -j])
        result["directions_m_n"] = directions
        result["third_quotient"] = degree == 2
        if len(directions) + int(degree == 2) != count:
            raise ArithmeticError("Primitive lattice enumeration disagrees with the count")
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("a", type=int)
    parser.add_argument("degree", type=int)
    parser.add_argument("--directions", action="store_true")
    args = parser.parse_args()
    try:
        print(json.dumps(spectrum(args.a, args.degree, args.directions), indent=2))
    except ValueError as error:
        parser.exit(2, f"Input error: {error}\n")


if __name__ == "__main__":
    main()
