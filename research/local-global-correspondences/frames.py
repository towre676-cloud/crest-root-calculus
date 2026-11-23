"""Rational PGL2 identifications of two smooth binary sextics, including their twist scalar."""
from functools import lru_cache
from math import comb, factorial

import sympy as sp

from frobenius import X
from lift_frames import align
from richelot import gram, match_frames, primitive, transform, vector


def falling(a, n):
    return factorial(a)//factorial(a-n) if a>=n else 0


def transvectant(f, m, g, n, r):
    """Affine polynomial of the order m+n-2r binary transvectant, up to scalar."""
    out=[0]*(m+n-2*r+1)
    for i in range(m+1):
        for j in range(n+1):
            if not (0<=i+j-r<len(out)):
                continue
            v=sum((-1)**k*comb(r,k)*falling(i,r-k)*falling(m-i,k)*
                  falling(j,k)*falling(n-j,r-k) for k in range(r+1))
            out[i+j-r]+=f.nth(i)*g.nth(j)*v
    result=sp.Poly.from_list(out[::-1],X,domain=sp.QQ)
    return result.monic() if result else result


@lru_cache(None)
def frame(f):
    I=transvectant(f,6,f,6,4)
    L=transvectant(f,6,I,4,4)
    M=transvectant(I,4,L,2,2)
    N=transvectant(I,4,M,2,2)
    triple=(L,M,N)
    if not sp.Matrix([vector(q) for q in triple]).det() or any(gram(triple)[i,i]==0 for i in range(3)):
        return None
    return triple


def identify(f,g,primes):
    """Find all matrices modulo scalar, and test the original binary forms exactly."""
    u,v=frame(f),frame(g)
    if (u is None)!=(v is None):
        return [],'covariant-rank exclusion'
    if u is None:
        matrices=align(f,g); method='height-certified split-prime branch triples'
    else:
        matrices=match_frames(u,v,primes); method='quadratic covariants'
    valid=[]
    for m in matrices:
        h=sp.Poly(transform(f,m,6),X,domain=sp.QQ)
        if h.monic()==g.monic():
            valid.append((m,h.LC()/g.LC()))
    return sorted(set(valid)),method
