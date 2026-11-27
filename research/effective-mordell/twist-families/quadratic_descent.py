"""Squareclasses of 3*a^2*U^4 + (7+4*s)*V^4 in Z[s], s^2=-2."""
from functools import lru_cache
from itertools import product
from math import isqrt

import sympy as sp
from descent import check, normalized_pairs, square_residues


def multiply(alpha,beta):
    a,b=alpha
    c,d=beta
    return (a*c-2*b*d,a*d+b*c)


@lru_cache(None)
def split_generator(p):
    for v in range(1,isqrt(p//2)+1):
        u=isqrt(p-2*v*v)
        if u*u+2*v*v==p:
            return u,v
    raise ArithmeticError(f'No split generator of norm {p}')


def candidates(d,a):
    factors=sp.factorint(6*d*a)
    values=[(1,0),(-1,0)]
    if factors.get(2,0)%2:
        values=[multiply(z,(0,1)) for z in values]
    for p,k in sorted(factors.items()):
        if p==2:
            continue
        if not k%2:
            # A primitive denominator divisible by 3 can give odd exponents
            # at both factors of 3. At other odd primes their common exponent is even.
            if p==3:
                values += [multiply(z,(3,0)) for z in values]
            continue
        if sp.legendre_symbol(-2,p)==-1:
            return []
        u,v=split_generator(p)
        values=[multiply(z,(u,w)) for z in values for w in (v,-v)]
    return sorted(set(values))


@lru_cache(None)
def squares_ramified(m):
    return {((u*u-2*v*v)%m,2*u*v%m) for u,v in product(range(m),repeat=2)}


def obstruction(a,delta,p,m):
    # Only prime moduli are used for inert p. At split p the two roots give
    # both factors of the residue ring, and at p=2 we enumerate the full ring.
    split=p!=2 and sp.legendre_symbol(-2,p)==1
    if split:
        root=int(sp.sqrt_mod(-2,m))
    for u,v in normalized_pairs(p,m):
        U,V=pow(u,4,m),pow(v,4,m)
        k,h=multiply(delta,((3*a*a*U+7*V)%m,4*V%m))
        if p==2:
            soluble=(k%m,h%m) in squares_ramified(m)
        elif split:
            soluble=all((k+sg*root*h)%m in square_residues(m) for sg in (1,-1))
        else:
            check(p==m,'inert test requires a prime modulus')
            soluble=(k*k+2*h*h)%p in square_residues(p)
        if soluble:
            return False
    return True


def sieve(d,a):
    places=((3,27),(2,16),(2,32),(3,81),(11,121),(17,289),
            (5,5),(7,7),(13,13),(19,19),(23,23),(29,29),(31,31),(41,41),(43,43))
    return [dict(delta=list(z),obstruction=next(
        (m for p,m in places if obstruction(a,z,p,m)),None)) for z in candidates(d,a)]
