"""Geometric isogeny and constituent-dimension tests from ordinary Weil reductions."""
from collections import Counter
from functools import lru_cache

import numpy as np
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_factor


@lru_cache(None)
def universal_power(d):
    bound = 4*d*d
    result = 1
    for prime in sp.primerange(2,bound+2):
        q = prime
        while (prime-1)*q <= bound:
            q *= prime
        result *= q
    return result


def trim(coeff):
    vals = list(coeff)
    while len(vals)>1 and not vals[-1]:
        vals.pop()
    return vals


def factors(poly,p):
    _, data = gf_factor(list(poly),p,ZZ)
    return [(tuple(map(int,f)),m) for f,m in data]


def matrix_power(mat,n,p):
    val = np.eye(len(mat),dtype=np.int64)
    while n:
        if n&1:
            val = val@mat%p
        n//=2
        if n:
            mat = mat@mat%p
    return val


@lru_cache(maxsize=32000)
def power_irreducible(poly, p, d):
    deg = len(poly)-1
    exponent = universal_power(d)%(p**deg-1)
    if exponent==0:
        return (1,p-1),deg
    if deg==1:
        return (1, -pow((-poly[1])%p,exponent,p)%p),1
    comp = np.eye(deg,k=-1,dtype=np.int64)
    comp[:,-1] = (-np.array(poly[:0:-1]))%p
    generator = matrix_power(comp, exponent,p)
    # Coefficients are Newton sums of Frobenius-conjugate nonzero roots.
    val = np.eye(deg,dtype=np.int64)
    seq, coeff = [], [1]
    for j in range(1,deg+1):
        val = val@generator%p
        seq.append(int(val.trace())%p)
        coeff.append(-sum(a*b for a,b in zip(coeff,reversed(seq)))*pow(j,-1,p)%p)
    powered = factors(coeff,p)
    if len(powered)!=1:
        raise ArithmeticError('The powered conjugate roots split into different Galois orbits')
    return powered[0]


def fingerprint(coeff, p, d):
    raw = trim(coeff)
    answer = Counter()
    for poly,m in factors(raw,p):
        target, mult = power_irreducible(poly,p,d)
        answer[target] += m*mult
    return tuple(sorted(answer.items()))
