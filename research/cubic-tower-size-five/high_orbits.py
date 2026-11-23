"""Isolate geometric simple constituents using the published six-index data at 193."""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from hashlib import sha256
from itertools import combinations
import json
from math import comb
import os
from pathlib import Path

if __name__ == '__main__':
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'

import numpy as np
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_sqf_list

from unit_orbits import ResidueRing, degree_blocks, twist_power

HERE = Path(__file__).resolve().parent
PRIME = 193
DIM = 192
EXPECTED_DIGEST = '13562c45b7164e3ba40b0752b9e3f82760cb4f3edf31b39a2a6f4d639e7d3260'


def read(path):
    return json.loads(path.read_text(encoding='utf-8'))


def labels(names, k):
    return {name for name in names if len(name.split('/')[0].split(',')) == k}


def fingerprint(coeff):
    """Return the large powered orbit's degree and three trace moments, if present."""
    p = PRIME
    if len(coeff) != DIM+1 or coeff[0] != 1 or any(type(c) != int or not 0 <= c < p for c in coeff):
        raise ArithmeticError('Invalid six-index numerator')
    last = next(i for i in range(DIM, -1, -1) if coeff[i])
    if last <= DIM//2:
        return None
    for squarefree, multiplicity in gf_sqf_list(coeff[:last+1], p, ZZ)[1]:
        if len(squarefree)-1 <= DIM//2:
            continue
        for n, h in degree_blocks(squarefree, p):
            if n <= DIM//2:
                continue
            if multiplicity != 1 or len(h)-1 != n:
                raise ArithmeticError('The unique large factor has unexpected multiplicity')
            ring = ResidueRing(h, p)
            y = ring.power(ring.x, twist_power(DIM) % (p**n-1))
            phi = ring.frobenius()
            current, previous = y, 0
            for k in sorted(n//prime for prime in sp.factorint(n)):
                for _ in range(k-previous):
                    current = current @ phi % p
                if np.array_equal(current, y):
                    return None
                previous = k
            moments = [n]
            weights = ring.trace_weights()
            current = y
            for j in range(3):
                moments.append(int(current @ weights) % p)
                if j < 2:
                    current = ring.mul(current, y)
            return tuple(moments)
    return None


def verify(workers):
    six = read(HERE/'six-index.json')
    names = [','.join(map(str,s))+'/'+sign for s in combinations(range(10),6) for sign in '+-']
    if six['prime'] != PRIME or six['dimension'] != DIM or set(names) != set(six['unit_polynomials']):
        raise ArithmeticError('The input must contain all 420 signed factors at 193')
    polys = [six['unit_polynomials'][name] for name in names]
    if workers == 1:
        profiles = list(map(fingerprint, polys))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            profiles = list(pool.map(fingerprint, polys, chunksize=2))
    groups = defaultdict(list)
    for name, profile in zip(names, profiles):
        if profile is not None:
            groups[profile].append(name)
    isolated = {matches[0]:profile[0] for profile,matches in groups.items() if len(matches)==1}
    encoded = json.dumps(isolated, sort_keys=True, separators=(',',':')).encode()
    digest = sha256(encoded).hexdigest()
    if digest != EXPECTED_DIGEST:
        raise ArithmeticError('The isolated supports or their certified dimension bounds differ')

    earlier = read(HERE.parent/'cubic-tower-geometric-separation/certificate.json')['simplicity']
    seven = read(HERE.parent/'cubic-tower-stage-seven/certificate.json')['simplicity']
    known_six = labels(earlier,6) | labels(seven,6) | set(six['simple_at_this_prime'])
    known_seven = labels(seven,7)
    five = read(HERE/'certificate.json')
    if len(five['simplicity']) != 2*comb(10,5) or five['unresolved']:
        raise ArithmeticError('The completed five-index classification is required')
    lower_count = 1+sum(2*comb(10,k) for k in range(1,6))
    lower_dimension = 1+sum(2*comb(10,k)*k*(1<<(k-1)) for k in range(1,6))
    known_count = lower_count+len(known_six)+len(known_seven)
    geometric_bound = lower_count+len(known_six | set(isolated))+len(known_seven)
    dimension_bound = (lower_dimension + DIM*len(known_six) + 448*len(known_seven)
                       + sum(n for name,n in isolated.items() if name not in known_six))
    results = {
        'six_index_isolated': len(isolated), 'smallest_certified_dimension': min(isolated.values()),
        'largest_certified_dimension': max(isolated.values()), 'six_index_known_simple': len(known_six),
        'known_simple_among_isolated': len(known_six & set(isolated)),
        'distinct_six_index_simple_classes_lower_bound': len(known_six | set(isolated)),
        'whole_simple_factors_certified': known_count, 'stage_ten_geometric_class_lower_bound': geometric_bound,
        'stage_ten_distinct_constituents_dimension_lower_bound': dimension_bound,
    }
    if tuple(results.values()) != (318,97,192,16,12,322,1293,1599,102691):
        raise ArithmeticError(f'Unexpected geometric totals: {results}')
    return results, isolated


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--workers', type=int, default=1)
    parser.add_argument('--list', action='store_true', help='list the dimension bound for each isolated support')
    args = parser.parse_args()
    summary, supports = verify(args.workers)
    if args.list:
        summary['support_dimension_lower_bounds'] = dict(sorted(supports.items()))
    print(json.dumps(summary, indent=2))
