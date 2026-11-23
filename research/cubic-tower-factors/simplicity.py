"""Exact sufficient tests for absolute simplicity and geometric nonisogeny."""
from functools import cache
from math import lcm

import numpy as np
import sympy as sp


@cache
def power_bound(degree):
    # phi(m) >= sqrt(m/2). A ratio between two roots of polynomials of
    # degree <= degree has at most degree^2 algebraic conjugates.
    bound = degree*degree
    limit = 2*bound*bound
    phi = list(range(limit+1))
    for prime in range(2, limit+1):
        if phi[prime] == prime:
            for n in range(prime, limit+1, prime):
                phi[n] -= phi[n]//prime
    return lcm(*(n for n in range(1, limit+1) if phi[n] <= bound))


@cache
def powered_reduction(coefficients, ell=1_000_003):
    """charpoly(F^M) mod ell, where coefficients list det(1-TF) from low to high."""
    degree = len(coefficients)-1
    if degree > 8 or not sp.isprime(ell) or not degree < ell <= 1_000_003:
        raise ValueError("Certificate arithmetic supports degrees up to 8 at primes below 1000004")
    modulus = power_bound(8)
    companion = np.zeros((degree, degree), dtype=np.int64)
    companion[1:, :-1] = np.eye(degree-1, dtype=np.int64)
    companion[:, -1] = [(-int(a)) % ell for a in reversed(coefficients[1:])]
    power = np.eye(degree, dtype=np.int64)
    while modulus:
        if modulus & 1:
            power = power @ companion % ell
        companion = companion @ companion % ell
        modulus >>= 1
    iterate = np.eye(degree, dtype=np.int64)
    traces, polynomial = [], [1]
    for j in range(1, degree+1):
        iterate = iterate @ power % ell
        traces.append(int(iterate.trace()) % ell)
        polynomial.append((-sum(a*b for a,b in zip(polynomial, reversed(traces)))
                           *pow(j, -1, ell)) % ell)
    return sp.Poly.from_list(polynomial, sp.Symbol("T"), modulus=ell)


def absolutely_simple(coefficients):
    ordinary = sp.Poly.from_list(coefficients, sp.Symbol("T"))
    if not ordinary.is_irreducible:
        return False
    return powered_reduction(tuple(coefficients)).is_sqf


def disjoint(a, b):
    f, g = powered_reduction(tuple(a)), powered_reduction(tuple(b))
    return f.gcd(g).degree() == 0
