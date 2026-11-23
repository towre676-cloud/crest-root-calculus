"""Verify the 504 five-index geometric certificates, with optional full matrix recount."""
import argparse
from collections import defaultdict
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_diff, gf_factor, gf_mul, gf_pow_mod, gf_rem

HERE = Path(__file__).resolve().parent
for folder in ('cubic-tower-factors','cubic-tower-cartier','cubic-tower-geometric-separation'):
    sys.path.insert(0,str(HERE.parent/folder))

from affine_frobenius import signed_numerators
from cartier import Cartier
from geometric_dimensions import root_of_unity_exponent, subdegrees
from tower import check_stage, trace_from_factor, trace_sums
from unit_orbits import observations

DIM = 80
KEYS = {','.join(map(str,s))+suffix:(s,sign) for s in combinations(range(10),5)
        for sign,suffix in enumerate(('/+','/-'))}


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def unique_dict(pairs):
    result = {}
    for key,value in pairs:
        if key in result:
            raise ValueError(f'Duplicate certificate key {key}')
        result[key] = value
    return result


def good_subsets(p):
    """Test the branch conditions over F_p independently of integer-resultant selection."""
    z = sp.Symbol('z')
    base = sp.Poly(-4*z**3+z**2+z+1,z,modulus=p)
    require(base.is_sqf,'Bad elliptic base prime')
    polys = [sp.Poly((12*i+1)*z*z+(12*i*i+1)*z+4*i**3+1,z,modulus=p) for i in range(10)]
    good = {i for i,q in enumerate(polys) if q.degree()==2 and q.is_sqf and q.gcd(base).degree()==0}
    disjoint = {pair for pair in combinations(sorted(good),2)
                if polys[pair[0]].gcd(polys[pair[1]]).degree()==0}
    return {s for s in combinations(range(10),5)
            if good.issuperset(s) and all(pair in disjoint for pair in combinations(s,2))}


def load_observations(receipt):
    require(receipt['dimension']==DIM and type(receipt['dimension']) is int, 'Wrong signed dimension')
    require(set(receipt['simplicity'])==set(KEYS), 'Missing or extra signed supports')
    data = {}
    for label,records in receipt['observations'].items():
        p = int(label)
        require(label==str(p) and DIM<p<=10000 and sp.isprime(p), 'Invalid witness prime')
        good = good_subsets(p)
        require(records and set(records)<=set(KEYS), 'Missing or unknown local coefficient records')
        for key,coeff in records.items():
            require(KEYS[key][0] in good, f'Bad reduction for {key} at {p}')
            require(isinstance(coeff,list) and len(coeff)==DIM+1 and coeff[0]==1 and
                    all(type(c) is int and 0<=c<p for c in coeff), 'Invalid mod-p coefficients')
            data[p,key] = coeff
    require(data,'Empty certificate')
    return data


def orbit_result(job):
    pair,coeff = job
    p,key = pair
    dims,signature = observations(coeff,p,DIM)
    return pair,(None if dims is None else frozenset(dims),tuple(signature))


def recounted(job):
    p,s,records = job
    rebuilt = signed_numerators(p,s)
    for sign,coeff in records.items():
        require(rebuilt[sign]==coeff, f'Local coefficient mismatch at {p}, {s}, sign {sign}')
    return len(records)


def parallel(function, jobs, workers, chunksize):
    if workers==1:
        yield from map(function,jobs)
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            yield from pool.map(function,jobs,chunksize=chunksize)


def ordered_primes(primes):
    return (isinstance(primes,list) and all(type(p) is int and DIM<p<=10000 for p in primes)
            and primes == sorted(set(primes)))


def certificates(receipt, data, results, *, complete=True):
    used = set()
    outstanding = {}
    for key,primes in receipt['simplicity'].items():
        require(ordered_primes(primes),'Simplicity witness primes must be increasing and distinct')
        dims = set(range(1,DIM//2+1))
        for p in primes:
            require((p,key) in data, f'Missing simplicity observation at {p}, {key}')
            used.add((p,key))
            permitted,_ = results[p,key]
            require(permitted is not None, f'Nonordinary simplicity witness at {p}, {key}')
            narrowed = dims.intersection(permitted)
            require(narrowed != dims, f'Redundant simplicity witness at {p}, {key}')
            dims = narrowed
        if dims:
            outstanding[key] = sorted(dims)
    require(receipt['unresolved']==outstanding,'Saved unresolved dimensions differ from the replay')
    if complete:
        require(not outstanding,'Some geometric proper dimensions survive')

    groups = [sorted(KEYS)]
    require(ordered_primes(receipt['separation']), 'Separation primes must be increasing and distinct')
    for p in receipt['separation']:
        check_stage(p,tuple(range(10)))
        split = []
        for group in groups:
            if len(group)==1:
                split.append(group)
                continue
            buckets = defaultdict(list)
            for key in group:
                require((p,key) in data, f'Missing comparison observation at {p}, {key}')
                used.add((p,key))
                buckets[results[p,key][1]].append(key)
            split.extend(buckets.values())
        groups = split
    require(used==set(data),'Stored observations do not equal the witness inputs')
    unresolved_pairs = sum(len(g)*(len(g)-1)//2 for g in groups)
    if complete:
        require(unresolved_pairs==0 and len(groups)==len(KEYS),'Some geometric pairs remain unresolved')
    return len(outstanding),unresolved_pairs


def independent_signature(p, coeff):
    """Use irreducible factorization and a Lagrange-interpolation trace formula."""
    last = max(i for i,c in enumerate(coeff) if c)
    traces = [0,0,0]
    exponent = root_of_unity_exponent(2*DIM)
    for factor,multiplicity in gf_factor(coeff[:last+1],p,ZZ)[1]:
        n = len(factor)-1
        y = gf_pow_mod([1,0],exponent%(p**n-1),factor,p,ZZ)
        deriv = gf_diff(factor,p,ZZ)
        power = y
        for j in range(3):
            # For separable monic f, Tr(a) is the x^(n-1) coefficient
            # of (a f') modulo f, by interpolation at the n roots of f.
            residue = gf_rem(gf_mul(power,deriv,p,ZZ),factor,p,ZZ)
            if len(residue)==n:
                traces[j] = (traces[j]+multiplicity*int(residue[0])) % p
            power = gf_rem(gf_mul(power,y,p,ZZ),factor,p,ZZ)
    return (last,*traces)


def sample_pairs(data):
    ordered = sorted(data)
    # Include the first and last prime, and both unit-rank cases (ordinary and
    # nonordinary cases) when present. Each selected pair uses the stored data.
    ordinary = next((k for k in ordered if data[k][-1]),None)
    deficient = next((k for k in ordered if not data[k][-1]),None)
    return sorted({k for k in (ordered[0],ordered[-1],ordinary,deficient) if k is not None})


def trace_characteristic(matrix,p):
    # Faddeev--LeVerrier for det(1-TM), valid here because p exceeds DIM.
    n = len(matrix)
    work = np.eye(n,dtype=np.int64)
    diagonal = np.diag_indices(n)
    coeffs = [1]
    for j in range(1,n+1):
        work = matrix @ work % p
        c = -int(work.trace())*pow(j,-1,p) % p
        coeffs.append(c)
        work[diagonal] = (work[diagonal]+c)%p
    return coeffs


def check_samples(data, results):
    for index,(p,key) in enumerate(sample_pairs(data)):
        s,sign = KEYS[key]
        check_stage(p,s)
        coeff = data[p,key]
        rebuilt = signed_numerators(p,s)
        require(rebuilt[sign]==coeff,f'Fresh matrix failed for {key} at {p}')
        if index==0:
            cp,cs = (p,s) if p<=997 else (83,tuple(range(5)))
            check_stage(cp,cs)
            direct = [trace_characteristic(matrix,cp) for matrix in Cartier(cp,5).matrices(cs)]
            baseline = rebuilt if (cp,cs)==(p,s) else signed_numerators(cp,cs)
            require(direct==baseline,'Independent Cartier matrices give different complete polynomials')
        permitted,signature = results[p,key]
        require(signature==independent_signature(p,coeff),'Independent powered traces differ')
        previous = subdegrees(p,tuple(coeff))
        require(permitted==(None if previous is None else frozenset(previous)),
                'Independent geometric subdimensions differ')
        # Direct fibre counts on the tower test actual Frobenius traces in both
        # signs and, at the first affordable prime, in the quadratic extension.
        extensions = (1,2) if index==0 and p*p<=300000 else (1,)
        for power in extensions:
            _,signed,_ = trace_sums(p,s,power,5)
            for e in (0,1):
                require((trace_from_factor(rebuilt[e],power)-signed[s][e])%p==0,
                        'Independent fibre-count trace differs')


def run(path,workers=1,recount=False):
    require(type(workers) is int and workers>=1,'Workers must be positive')
    receipt = json.loads(Path(path).read_text(encoding='utf-8'),object_pairs_hook=unique_dict)
    data = load_observations(receipt)
    results = dict(parallel(orbit_result,sorted(data.items()),workers,8))
    certificates(receipt,data,results)
    check_samples(data,results)
    if recount:
        jobs = defaultdict(dict)
        for (p,key),coeff in data.items():
            s,sign = KEYS[key]
            jobs[p,s][sign] = coeff
        total = sum(parallel(recounted,[(p,s,records) for (p,s),records in sorted(jobs.items())],workers,4))
        require(total==len(data),'Local matrix recount missed a coefficient record')
    return {'signed_supports':len(KEYS),'geometrically_separated_pairs':len(KEYS)*(len(KEYS)-1)//2,
            'stored_local_polynomials':len(data),'complete_matrix_recount':recount}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--certificate',type=Path,default=HERE/'certificate.json')
    parser.add_argument('--workers',type=int,default=1)
    parser.add_argument('--recount',action='store_true')
    options=parser.parse_args()
    print(json.dumps(run(options.certificate,options.workers,options.recount),indent=2))
