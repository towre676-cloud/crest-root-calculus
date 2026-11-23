"""Signed Hasse--Witt numerators of the cubic tower by finite-fibre Frobenius."""
from functools import lru_cache
from itertools import product
import numpy as np
import sympy as sp


F = np.array([1, 1, 1, -4], dtype=np.int64)


def trim(a):
    n = np.flatnonzero(a)
    return a[:int(n[-1])+1] if len(n) else np.zeros(1, dtype=np.int64)


def multiply(a, b, p):
    return trim(np.convolve(a, b) % p)


@lru_cache(maxsize=256)
def f_power(p, exponent):
    if exponent == 0:
        return np.array([1], dtype=np.int64)
    sub = f_power(p, exponent//2)
    val = multiply(sub, sub, p)
    return multiply(val, F, p) if exponent % 2 else val


def subtract(a, b, p):
    value = np.zeros(max(len(a), len(b)), dtype=np.int64)
    value[:len(a)] = a
    value[:len(b)] -= b
    return trim(value % p)


@lru_cache(maxsize=512)
def cubic_images(p, index):
    """Columns t^p and t^(2p), suppressing the determined parity factor of w."""
    c = [np.array([1], dtype=np.int64), np.array([0]), np.array([0])]
    linear = np.array([3*index, 3], dtype=np.int64)
    result = []
    for n in range(2*p):
        w_c2 = multiply(c[2], F, p) if n % 2 else c[2]
        c = [w_c2, subtract(c[0], multiply(c[2], linear, p), p), c[1]]
        if n+1 in (p, 2*p):
            result.append(c[1:])
    return result


@lru_cache(maxsize=2048)
def half_column(p, support, exponents, row):
    value = np.array([1],dtype=np.int64)
    for i,e,t in zip(support,exponents,row):
        value = multiply(value,cubic_images(p,i)[e-1][t-1],p)
    return value


def characteristic(matrix, p):
    """Ascending coefficients of det(1-TM) via exact Hessenberg similarity."""
    h = matrix.copy()
    n = len(h)
    for j in range(n-2):
        nz = np.flatnonzero(h[j+1:,j])
        if not len(nz):
            continue
        pivot = j+1+int(nz[0])
        if pivot != j+1:
            h[[pivot,j+1]] = h[[j+1,pivot]]
            h[:,[pivot,j+1]] = h[:,[j+1,pivot]]
        ratios = h[j+2:,j]*pow(int(h[j+1,j]),-1,p)%p
        h[j+2:,j:] = (h[j+2:,j:]-ratios[:,None]*h[j+1,j:])%p
        h[:,j+1] = (h[:,j+1]+h[:,j+2:]@ratios)%p
    # The leading principal minors of a Hessenberg matrix satisfy the
    # first-column cofactor recurrence, written for det(1-TM).
    polys = [np.ones(1,dtype=np.int64)]
    for i in range(n):
        val = np.zeros(i+2,dtype=np.int64)
        val[:-1]=polys[i]
        val[1:]-=h[i,i]*polys[i]
        sub=1
        for j in range(i-1,-1,-1):
            sub=sub*int(h[j+1,j])%p
            if not sub:
                break
            val[i+1-j:] -= (int(h[j,i])*sub % p)*polys[j]
        polys.append(val%p)
    return list(map(int,polys[-1]))


def signed_numerators(p, support):
    """Ascending coefficients of the two reciprocal polynomials modulo p.

    The caller verifies simultaneous good reduction. The prime exceeds the
    factor dimension so Newton division on the lower half is invertible.
    """
    support = tuple(support)
    k = len(support)
    dim = k*2**(k-1)
    if not sp.isprime(p) or p<=max(dim,3) or p>10000 or not 1<=k<=7:
        raise ValueError("Use prime 3,dimension < p <= 10000 and 1 to 7 indices")
    basis = list(product((1,2), repeat=k))
    rad = [-(12*i+1)*pow(12,-1,p) % p for i in support]
    chi = [int(sp.legendre_symbol(v,p)) for v in rad]
    if 0 in chi:
        raise ValueError("Degenerate infinity fibre")
    parity = [sum(e)%2 for e in basis]
    weights = [[(sum(e)+3*((sum(e)+sign)%2))//2 for e in basis] for sign in (0,1)]
    slots=[]
    for weight in weights:
        stops=np.cumsum([0,*weight])
        slots.append([slice(int(a),int(b)) for a,b in zip(stops,stops[1:])])
    matrices=[np.zeros((sum(w),sum(w)),dtype=np.int64) for w in weights]
    max_weight=max(map(max,weights))
    wanted=p*np.arange(1,max_weight+1)[None,:]-np.arange(1,max_weight+1)[:,None]
    split=k//2
    for j, source in enumerate(basis):
        for i, target in enumerate(basis):
            odd = sum((e-t)%2 for e,t in zip(source,target))
            left = half_column(p,support[:split],source[:split],target[:split])
            right = half_column(p,support[split:],source[split:],target[split:])
            base = multiply(left,right,p)
            for d in (0,1):
                sign=d^parity[j]
                extra=d*((p-1)//2)+(d+odd)//2
                coef=multiply(base,f_power(p,extra),p) if extra else base
                indices=wanted[:weights[sign][i],:weights[sign][j]]
                val=np.zeros(indices.shape,dtype=np.int64)
                valid=indices<len(coef)
                val[valid]=coef[indices[valid]]
                matrices[sign][slots[sign][i],slots[sign][j]]=val
    result = []
    for sign,frob in enumerate(matrices):
        coeffs = characteristic(frob,p)[:dim+1]
        # The infinity permutation has order at most two. Its sign summand
        # has dimension 2^(k-1), and its trace determines both multiplicities.
        r = 2**(k-1)
        trace = (int(np.prod([1+c for c in chi]))+
                 (-1)**sign*int(np.prod([1-c for c in chi])))//2
        for scalar, repeated in ((1,(r+trace)//2),(-1,(r-trace)//2)):
            for _ in range(repeated):
                for j in range(1,len(coeffs)):
                    coeffs[j] = (coeffs[j]+scalar*coeffs[j-1]) % p
        result.append(coeffs)
    return result


if __name__ == '__main__':
    import argparse
    import json
    from pathlib import Path
    import sys

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('prime', type=int)
    parser.add_argument('indices', nargs='+', type=int)
    args = parser.parse_args()
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / 'cubic-tower-factors'))
    from tower import check_stage

    check_stage(args.prime, args.indices)
    plus, minus = signed_numerators(args.prime, args.indices)
    print(json.dumps({'prime':args.prime, 'support':args.indices, 'plus':plus, 'minus':minus}, indent=2))
