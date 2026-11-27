"""Resolve rational-abscissa components of z² = d(x⁶ − 1) by Eisenstein descent.

solve_residual_pair(d, a, b) returns a finite set of all nonbranch rational
abscissas in that component, or None if this test leaves a possibility open.
The ordered squareclasses a,b belong to x²+x+1 and x²−x+1, respectively.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import isqrt, prod

from cypari2 import Pari
from sympy import sqrt_mod

pari = Pari()
MODULI = ((2, 16), (3, 9), (7, 7), (13, 13), (19, 19), (97, 97))
OMEGA_UNITS = ((1, 0), (-1, 0), (0, 1), (0, -1), (1, 1), (-1, -1))

@lru_cache(maxsize=256)
def factor(n):
    primes, powers = pari.factorint(n)
    if not all(pari.isprime(q) for q in primes):
        raise ArithmeticError("Unproved rational prime factor")
    return tuple((int(q), int(v)) for q, v in zip(primes, powers))

def squarefree(n):
    return prod(q for q, v in factor(n) if v % 2)

def mul(x, y):
    a, b = x
    c, d = y
    return (a*c-b*d, a*d+b*c-b*d)

def conjugate(x):
    return (x[0]-x[1], -x[1])

def norm(x):
    a, b = x
    return a*a-a*b+b*b

def reduce_mod(x, m):
    return (x[0] % m, x[1] % m)

@lru_cache(maxsize=256)
def split_prime(q):
    """Choose the lexicographically first nonnegative Eisenstein generator of norm q."""
    for root in sqrt_mod(-3, q, all_roots=True):
        r, t = q, int(root)
        while t*t > q:
            r, t = t, r % t
        u = isqrt((q-t*t)//3)
        if t*t+3*u*u != q:
            continue
        w = (t+u, 2*u)
        associates = [mul(unit, v) for unit in OMEGA_UNITS for v in (w, conjugate(w))]
        return min(v for v in associates if min(v) >= 0 and norm(v) == q)
    raise ArithmeticError("Cornacchia failed for split prime")

def representatives(n, k, generators=None):
    """Superset of squareclasses of an integral Eisenstein number of norm n times a square."""
    reps = [(1, 0), (-1, 0)]
    for q, _ in factor(6*k*n):
        if q == 3:
            choices = [(1, 2)] if n % 3 == 0 else [(1, 0)]
        elif q % 3 == 2:
            if n % q == 0:
                return []
            choices = [(1, 0), (q, 0)] if 6*k % q == 0 else [(1, 0)]
        elif n % q == 0:
            prime = generators[q] if generators and q in generators else split_prime(q)
            if norm(prime) != q:
                raise ArithmeticError("Invalid split generator")
            choices = [prime, conjugate(prime)]
        else:
            # Both split valuations have the same parity by the norm condition.
            choices = [(1, 0), (q, 0)] if 6*k % q == 0 else [(1, 0)]
        reps = [mul(old, item) for old in reps for item in choices]
    return reps

@lru_cache(maxsize=None)
def square_residues(m):
    return {reduce_mod(mul(w, w), m) for w in product(range(m), repeat=2)}

def local_filter(pairs, a, b, k, prime, modulus):
    square_int = {i*i % modulus for i in range(modulus)}
    square_eis = square_residues(modulus)
    left = list(pairs)
    survivors = set()
    charts = [(j, 1) for j in range(modulus)]
    charts.extend((1, j) for j in range(0, modulus, prime))
    for u, v in charts:
        U, V = u*u, k*v*v
        if (a*(U*U+3*V*V) % modulus not in square_int or
            b*(3*U*U+V*V) % modulus not in square_int):
            continue
        beta = (U+V, 2*V)
        gamma = (V+U, 2*U)
        for delta, eta in left:
            if (reduce_mod(mul(delta, beta), modulus) in square_eis and
                reduce_mod(mul(eta, gamma), modulus) in square_eis):
                survivors.add((delta, eta))
    return survivors

def sieve(a, b, k, generators=None):
    pairs = set(product(representatives(a, k, generators), representatives(b, k, generators)))
    counts = [len(pairs)]
    for prime, modulus in MODULI:
        pairs = local_filter(pairs, a, b, k, prime, modulus)
        counts.append(len(pairs))
    return sorted(pairs), counts

# Each line is checked with elliptic Chabauty by final_three.m (power-of-two index obstruction).
CHABAUTY = {
    (19, 1, 19): (((1, 0), (-3, -5)), (Fraction(1, 4), Fraction(4))),
    (95, 1, 57): (((-2, 0), (16, 2)), (Fraction(1), Fraction(225))),
    (182, 7, 13): (((-2, -6), (6, 8)), (Fraction(4),)),
}

def square_rational(x):
    return x>=0 and isqrt(x.numerator)**2==x.numerator and isqrt(x.denominator)**2==x.denominator

@lru_cache(maxsize=128)
def positive_pair(d, a, b):
    k = squarefree(d*a*b)
    remaining, _ = sieve(a, b, k)
    if not remaining:
        return frozenset()
    key = (d, a, b)
    if key not in CHABAUTY:
        for prime, modulus in ((7,49),(61,61),(193,193)):
            remaining=local_filter(remaining,a,b,k,prime,modulus)
            if not remaining:
                return frozenset()
        return None
    pair, ts = CHABAUTY[key]
    if remaining != [pair]:
        raise ArithmeticError("Unexpected surviving Eisenstein class")
    xs=set()
    for t in ts:
        if t == k:
            continue
        x=(k+t)/(k-t)
        if (square_rational(d*(x**6-1)) and square_rational(a*(x*x+x+1)) and
            square_rational(b*(x*x-x+1))):
            xs.add(x)
    return frozenset(xs)

def solve_residual_pair(d, a, b):
    if d==0 or min(a,b)<1:
        raise ValueError("Nonzero twist and positive squareclasses required")
    reflected=a>b
    xs=positive_pair(abs(d), min(a,b), max(a,b))
    if xs is None:
        return None
    values={-x if reflected else x for x in xs}
    if d<0:
        values={1/x for x in values if x}
    return values

RAY_MODULUS=64*81*7
RAY_GENERATOR=(1,7)

@lru_cache(maxsize=128)
def ray_generator(p):
    """A generator proving p belongs to the explicit ray, if one exists."""
    if p<2 or p%RAY_MODULUS!=43 or not pari.isprime(p):
        return None
    g=split_prime(p)
    for root in (g,conjugate(g)):
        for unit in OMEGA_UNITS:
            representative=mul(unit,root)
            if reduce_mod(representative,RAY_MODULUS)==RAY_GENERATOR:
                return representative
    return None

def inert_supported(d):
    """Every odd split prime in the primitive abscissa numerators must divide d or 3."""
    return all(q==3 or q%3==2 for q,_ in factor(abs(d)))
