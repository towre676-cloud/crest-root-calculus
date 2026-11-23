"""Rank of prime-field Frobenius traces on every signed support of stage ten."""
import argparse
from hashlib import sha256
import json
from math import gcd
import os
from pathlib import Path
import sys

# The many thin matrix products run fastest with one BLAS worker.
if __name__ == '__main__':
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
TOWER = HERE.parent / 'cubic-tower-factors'
sys.path.insert(0, str(TOWER))
from tower import trace_sums


def good_primes(count):
    z = sp.Symbol('z')
    f = -4*z**3+z**2+z+1
    earlier = [f]
    bad = 6*sp.discriminant(f, z)
    for i in range(10):
        h = (12*i+1)*z**2+(12*i*i+1)*z+4*i**3+1
        bad *= (12*i+1)*sp.discriminant(h, z)
        for other in earlier:
            bad *= sp.resultant(other, h, z)
        earlier.append(h)
    bad = int(bad)
    primes = []
    for p in sp.primerange(5, 100000):
        if gcd(p, bad) == 1:
            primes.append(p)
            if len(primes) == count:
                return primes
    raise ValueError('Too many primes requested')


def power_mod(vector, exponent, p):
    out = np.ones_like(vector)
    while exponent:
        if exponent & 1:
            out = out*vector % p
        vector = vector*vector % p
        exponent //= 2
    return out


def subset_products(chars):
    out = np.empty((1 << len(chars), chars.shape[1]), dtype=np.float64)
    out[0] = 1
    for i, values in enumerate(chars):
        n = 1 << i
        np.multiply(out[:n], values, out=out[n:2*n])
    return out


def prime_traces(p):
    """Column zero is E, followed by + and - for each increasing nonzero mask."""
    x = np.arange(p, dtype=np.int64)
    sq = x*x % p
    root = np.full(p, -1, dtype=np.int64)
    root[sq] = x
    nonsquare = int(np.flatnonzero(root < 0)[0])
    inv_s = pow(nonsquare, -1, p)
    cube = sq*x % p
    freq = np.array([np.bincount((cube+c*x) % p, minlength=p) for c in (0, 1, nonsquare)])
    typ = np.where(root < 0, 2, 1)
    typ[0] = 0
    scale_root = np.where(typ == 2, root[x*inv_s % p], root)
    scale = power_mod(scale_root, p-4, p)  # inverse cube on nonzero elements
    scale[0] = 1
    f = (((-4*x+1)*x+1)*x+1) % p
    indices = np.arange(10)
    rad_inf = (-(12*indices+1)*pow(12, -1, p)) % p
    chi_inf = np.where(root[rad_inf] >= 0, 1, -1)
    total = []
    for sign, d in ((1, 1), (-1, inv_s)):
        roots = root[f*d % p]
        valid = roots >= 0
        z, w = x[valid], roots[valid]
        weight = np.where(w == 0, 1, 2)
        if sign == 1:
            base_trace = p-int(weight.sum())
        a = 3*(indices[:, None]+z)*d % p
        b = w*d % p
        chars = freq[typ[a], b*scale[a] % p]-1
        chars = np.column_stack((chars, 1+sign*chi_inf))
        weight = np.append(weight, 1)
        # Sum all ten-index subset products via the two five-index tables.
        # Every floating-point multiply and sum is an exact integer: in each
        # entry the sum of absolute values is at most (2p+1) * 2**10 < 2**53.
        left = subset_products(chars[:5])
        right = subset_products(chars[5:])
        total.append(-((right*weight) @ left.T).ravel().astype(np.int64))
    plus, minus = total[0]+total[1], total[0]-total[1]
    if np.any(plus % 2):
        raise ArithmeticError('Nonintegral signed traces')
    return np.r_[base_trace, np.column_stack((plus[1:]//2, minus[1:]//2)).ravel()]


def independent_rows(panel, ell):
    """Select a maximal pivot row set modulo ell in this narrow panel."""
    work = panel.astype(np.int64)
    perm = np.arange(len(work))
    r = 0
    for j in range(work.shape[1]):
        nz = np.flatnonzero(work[r:, j])
        if not len(nz):
            continue
        k = r+int(nz[0])
        work[[r,k]] = work[[k,r]]
        perm[[r,k]] = perm[[k,r]]
        row = work[r,j:]*pow(int(work[r,j]), -1, ell) % ell
        work[r+1:,j:] = (work[r+1:,j:]-work[r+1:,j:j+1]*row) % ell
        r += 1
    return perm, r


def inverse_mod(square, ell):
    n = len(square)
    aug = np.column_stack((square.astype(np.int64), np.eye(n, dtype=np.int64)))
    for j in range(n):
        k = j+int(np.flatnonzero(aug[j:, j])[0])
        aug[[j,k]] = aug[[k,j]]
        row = aug[j]*pow(int(aug[j,j]), -1, ell) % ell
        aug = (aug-aug[:,j:j+1]*row) % ell
        aug[j] = row
    return aug[:,n:]


def certify_full_column_rank(matrix, ell=127, block=32):
    """Exact finite-field Schur complements; reject the first deficient panel."""
    work = (matrix % ell).astype(np.float64)
    cols = work.shape[1]
    while work.shape[1]:
        n = min(block, work.shape[1])
        perm, rank = independent_rows(work[:, :n], ell)
        if rank != n:
            raise ArithmeticError('Trace matrix is deficient modulo the witness prime')
        work = work[perm]
        inv = inverse_mod(work[:n,:n], ell)
        multipliers = (work[n:,:n] @ inv) % ell
        # Each dot product has length at most 32 with residues at most 126,
        # so binary64 performs only exact integer operations here as well.
        work = (work[n:,n:] - multipliers @ work[:n,n:]) % ell
    return cols


def verify():
    primes = good_primes(2055)
    matrix = np.asarray([prime_traces(p) for p in primes], dtype=np.int32)
    for i, p in enumerate(primes[:4]):
        base, signed, _ = trace_sums(p, list(range(10)), 1, 10)
        if matrix[i, 0] != base:
            raise ArithmeticError('Base trace disagrees')
        for support, pair in signed.items():
            mask = sum(1 << index for index in support)
            if list(matrix[i,2*mask-1:2*mask+1]) != list(pair):
                raise ArithmeticError('Independent fibre counter disagrees')
    return {'number_of_good_primes': len(primes), 'first_good_prime': primes[0],
            'last_good_prime': primes[-1], 'columns': matrix.shape[1], 'coefficient_prime': 127,
            'rank': certify_full_column_rank(matrix), 'trace_matrix_int32_little_endian_sha256':
            sha256(matrix.astype('<i4').tobytes(order='C')).hexdigest()}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--compute', action='store_true')
    parser.add_argument('--primes', action='store_true', help='print the primes for the C++ cross-check')
    args = parser.parse_args()
    if args.primes:
        print(*good_primes(2055), sep='\n')
    else:
        result = verify()
        if args.compute:
            print(json.dumps(result, indent=2))
        elif result != json.loads((HERE/'rank.json').read_text(encoding='utf-8')):
            raise ArithmeticError('Rank receipt differs')
        else:
            print('Verified 2,047 independent rational isogeny classes in the stage-ten tower.')
