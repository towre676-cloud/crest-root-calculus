"""Exact lower density of prime radicands certified by the local Eisenstein tests."""
from collections import Counter
from functools import lru_cache
from fractions import Fraction
from itertools import product
from math import gcd, prod

from eisenstein_quadratic import (conjugate, local_filter, mul, norm, pari,
                                  reduce_mod, split_prime)

PRIMES=(2,3,5,7,13,17,19,97)
FULL_MASK=(1<<128)-1

@lru_cache(maxsize=8)
def squareclass_representatives(q):
    # Unit squares modulo 8 characterize dyadic units that are squares locally,
    # since 1+8O consists of squares. At 3,5,... residue-unit Hensel lifting suffices.
    modulus=8 if q==2 else q
    units={z for z in product(range(modulus),repeat=2) if gcd(norm(z),modulus)==1}
    squares={reduce_mod(mul(z,z),modulus) for z in units}
    reps=[]
    while units:
        z=min(units)
        reps.append(z)
        units.difference_update(reduce_mod(mul(u,z),modulus) for u in squares)
    return tuple(reps)

def rational_local(norm_g,a,b,q,m):
    squares={j*j % m for j in range(m)}
    lines=[(j,1) for j in range(m)]+[(1,j) for j in range(0,m,q)]
    for u,v in lines:
        forms=[a*(u*u+u*v+v*v),b*(u*u-u*v+v*v),norm_g*a*b*(u*u-v*v)]
        if all(z % m in squares for z in forms):
            return True
    return False

def scale(c,pair):
    return tuple(c*v for v in pair)

def local_mask(q,g):
    modulus_z={2:64,3:81}.get(q,q)
    modulus_eis={2:16,3:9}.get(q,q)
    norm_g=norm(g)
    mask=0
    offset=0
    for a,v in product((1,3),repeat=2):
        # Actual ordered squareclasses are (a,pv); k=sf(a*v), independent of p.
        b=norm_g*v
        k=1 if a==v else 3
        signs=(1,-1,2,-2)
        ds=[scale(i,(1,2) if a==3 else (1,0)) for i in signs]
        es=[scale(i,mul((1,2) if v==3 else (1,0),root))
            for i in signs for root in (g,conjugate(g))]
        pairs=list(product(ds,es))
        if rational_local(norm_g,a,b,q,modulus_z):
            survivors=local_filter(pairs,a,b,k,q,modulus_eis)
            for idx,pair in enumerate(pairs):
                if pair in survivors:
                    mask|=1<<(offset+idx)
        offset+=len(pairs)
    if offset!=128:
        raise ArithmeticError("Missing class labels")
    return mask

def certify_prime(p):
    """True proves the four-point lists over both ±p for p outside 2,3 and 7 etc."""
    if not pari.isprime(p) or p in (2,3):
        return False
    if p%3==2:
        return True
    if p in PRIMES:
        return False
    g=split_prime(p)
    mask=FULL_MASK
    for q in PRIMES:
        mask &= local_mask(q,g)
        if not mask:
            return True
    return False

def density():
    live=Counter({FULL_MASK:1})
    receipt=[]
    for q in PRIMES:
        options=Counter(local_mask(q,g) for g in squareclass_representatives(q))
        next_counts=Counter()
        for x,n in live.items():
            for y,m in options.items():
                next_counts[x & y] += n*m
        live=next_counts
        receipt.append((q,sum(live.values()),live[0]))
    total=sum(live.values())
    return receipt,Fraction(total+live[0],2*total)

if __name__=='__main__':
    rows,rate=density()
    expected=[(2,8,2),(3,16,6),(5,32,12),(7,128,82),(13,512,382),
              (17,1024,774),(19,4096,3224),(97,16384,13280)]
    if rows!=expected or rate!=Fraction(927,1024):
        raise ArithmeticError('Prime signature calculation disagrees with the receipt')
    for q,branches,closed in rows:
        print('p-local test:',q,'signatures:',branches,'empty:',closed)
    print('Natural density of the certified prime family:',rate)
