"""Search target coefficients over finite fields, lift, and certify in characteristic zero."""
from fractions import Fraction as Q
from math import gcd, isqrt
from pathlib import Path
import sys

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
# In the repository this program is a sibling of the previously released jet test.
ELLIPTIC = HERE.parent
sys.path.insert(0, str(ELLIPTIC / "degree-four-jets"))
from reconstruct import CUBICS, S, T, branch_powers, differential, reconstruct


ORDER = 27


def kernel_mod(matrix, p):
    """Nullspace over F_p as a matrix with independent columns."""
    a = np.asarray(matrix, dtype=np.int64).copy() % p
    pivots = []
    for j in range(a.shape[1]):
        i = len(pivots)
        candidates = np.flatnonzero(a[i:,j])
        if not len(candidates):
            continue
        k = i + int(candidates[0])
        a[[i,k]] = a[[k,i]]
        row = a[i]*pow(int(a[i,j]), -1, p) % p
        a = (a-a[:,j:j+1]*row) % p
        a[i] = row
        pivots.append(j)
    free = [j for j in range(a.shape[1]) if j not in pivots]
    coords = np.zeros((a.shape[1],len(free)), dtype=np.int64)
    coords[free] = np.eye(len(free), dtype=np.int64)
    coords[pivots] = -a[:len(pivots),free] % p
    return coords


def residue(x, p):
    q = Q(x)
    return q.numerator*pow(q.denominator,-1,p) % p


def row_reduction(powers, p):
    """Eliminate the ten numerator coefficients and the first two denominator jets."""
    values = np.asarray([[0 if n < i else residue(powers[j][n-i],p) for i,j in CUBICS]
                         for n in range(ORDER)], dtype=np.int64)
    basis = values @ kernel_mod(values[:2],p) % p
    relations = kernel_mod(values[:-2].T,p).T
    if basis.shape[1] != 8 or len(relations) != ORDER-12:
        raise ValueError("Source cubics have exceptional truncated rank at this prime")
    circuit = []
    for k in range(ORDER):
        shifted = np.zeros_like(basis)
        shifted[k:] = basis[:ORDER-k]
        circuit.append((relations @ shifted[2:] % p).ravel())
    return np.asarray(circuit)


def jets(h, a, b, p):
    """Simultaneously compute the local Weierstrass expansions for a batch."""
    h = [residue(h[i],p) for i in range(ORDER)]
    h2 = [sum(h[j]*h[i-j] for j in range(i+1)) % p for i in range(ORDER)]
    one = np.ones_like(a)
    r, square, rhs, deriv = [one], [one], [one], [-2*one]
    for n in range(1,ORDER):
        square.append(sum(r[j]*r[n-j] for j in range(1,n)) % p)
        val = sum(square[j]*r[n-j] for j in range(1,n+1))
        if n >= 4:
            val = val+a*r[n-4]
        if n == 6:
            val = val+b
        rhs.append(val % p)
        known = (sum(deriv[j]*deriv[n-j] for j in range(1,n))
                 - 4*sum(h2[j]*rhs[n-j] for j in range(n+1)))
        rn = (known % p)*pow(4*(n+1),-1,p) % p
        r.append(rn)
        square[n] = (square[n]+2*rn) % p
        rhs[n] = (rhs[n]+3*rn) % p
        deriv.append((n-2)*rn % p)
    return np.asarray(r).T


def singular_minor(matrices, inv, p):
    """Exact batched elimination of the first seven-by-seven minor."""
    a = matrices.copy()
    batch = np.arange(len(a))
    deficient = np.zeros(len(a),dtype=bool)
    for j in range(7):
        rows = a[:,j:,j] != 0
        deficient |= ~np.any(rows,axis=1)
        k = j+np.argmax(rows,axis=1)
        a[batch,j], a[batch,k] = a[batch,k].copy(), a[batch,j].copy()
        scale = a[:,j+1:,j]*inv[a[:,j,j]][:,None] % p
        a[:,j+1:,j+1:] = (a[:,j+1:,j+1:]-scale[:,:,None]*a[:,None,j,j+1:]) % p
    return deficient


def scan_prime(p, h, circuit):
    if p <= ORDER or not sp.isprime(p) or p > 1021:
        raise ValueError("Use a prime between 29 and 1021")
    inv = np.array([0]+[pow(i,-1,p) for i in range(1,p)])
    # Every product and sum in this batched linear map fits in signed int64.
    minor = circuit.reshape(ORDER,-1,8)[:,:7,:7].reshape(ORDER,49)
    answers = []
    for lo in range(0,p*p,4096):
        ab = np.arange(lo,min(lo+4096,p*p))
        a,b = np.divmod(ab,p)
        ok = (4*a**3+27*b**2)%p != 0
        a,b = a[ok],b[ok]
        series = jets(h,a,b,p)
        deficient = singular_minor(((series@minor)%p).reshape(-1,7,7), inv, p)
        for aa,bb,rr in zip(a[deficient],b[deficient],series[deficient]):
            m = (rr@circuit % p).reshape(-1,8)
            if kernel_mod(m,p).shape[1] >= 2:
                answers.append((int(aa),int(bb)))
    return answers


def rational_reconstruction(value, modulus):
    """Find the unique fraction with |num|, den <= floor(sqrt(modulus/2))."""
    n0,n1 = modulus,value
    d0,d1 = 0,1
    bound = isqrt(modulus//2)
    while abs(n1)>bound:
        q,r = divmod(n0,n1)
        n0,n1,d0,d1 = n1,r,d1,d0-q*d1
    if not d1 or abs(d1)>bound or gcd(n1,d1)!=1 or gcd(d1,modulus)!=1:
        return None
    return Q(n1,d1)


def lift_target(chart, u, v, primes, verbose=False):
    terms, powers = branch_powers(chart)
    h = differential(terms,powers,Q(u),Q(v))
    crt = [0,0]
    modulus = 1
    accepted = []
    previous = None
    for p in primes:
        p=int(p)
        try:
            candidates = scan_prime(p,h,row_reduction(powers,p))
        except ValueError as exc:
            if 'invertible' in str(exc) or 'exceptional' in str(exc):
                continue
            raise
        if verbose:
            print(p, candidates[:8], flush=True, file=sys.stderr)
        if len(candidates) != 1:
            continue
        pair = candidates[0]
        inv = pow(modulus,-1,p)
        crt = [v0+modulus*((new-v0)*inv%p) for v0,new in zip(crt,pair)]
        modulus *= p
        accepted.append([p,*pair])
        lifted = tuple(rational_reconstruction(z,modulus) for z in crt)
        if None not in lifted and lifted == previous:
            proof = reconstruct(chart,u,v,*lifted)
            if proof['degree_four']:
                return {'target':[str(x) for x in lifted], 'primes_and_residues':accepted,
                        'reconstruction':proof}
        previous = lifted
    return None
