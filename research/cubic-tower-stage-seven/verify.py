"""Prove the higher constituents needed for the complete seventh-stage decomposition."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
import os
from pathlib import Path
import sys

if __name__ == '__main__':
    os.environ['OPENBLAS_NUM_THREADS'] = '1'
    os.environ['MKL_NUM_THREADS'] = '1'

HERE=Path(__file__).resolve().parent
SIEVE=HERE.parent/'cubic-tower-size-five'
sys.path.insert(0,str(SIEVE))
sys.path.insert(0,str(HERE.parent/'cubic-tower-factors'))
from block_frobenius import signed_blocks
from unit_orbits import observations
from tower import check_stage, trace_from_factor, trace_sums


def label(S,e):
    return ','.join(map(str,S))+'/'+('+-'[e])


def parse(name):
    s,sign=name.split('/')
    return tuple(map(int,s.split(','))),'+-'.index(sign)


def inspect(job):
    p,S,e,c=job
    d=len(S)*(1<<(len(S)-1))
    if not len(c)==d+1 or any(type(x)!=int or not 0<=x<p for x in c):
        raise ArithmeticError('Noncanonical reported unit coefficients')
    if c[0]!=1:
        raise ArithmeticError('Invalid constant coefficient')
    return observations(c,p,d)


def recount(job):
    p,S,rows=job
    vals=signed_blocks(p,S)
    if any(vals[sign]!=coeff for sign,coeff in rows.items()):
        raise ArithmeticError(f'Unit numerator differs at {p} on {S}')


def run(data,recheck,workers):
    universe=tuple(range(7));six=list(combinations(universe,6))
    target={label(S,e) for S in six[1:]+[universe] for e in (0,1)}
    if set(data['simplicity'])!=target:
        raise ArithmeticError('Unexpected scope of the higher simplicity proof')
    rows=data['polynomials']
    used={}
    def get(p,name):
        S,e=parse(name)
        check_stage(p,S)
        coeff=rows[str(p)][name]
        used[p,name]=(p,S,e,coeff)
        return coeff
    for name,primes in data['simplicity'].items():
        for p in primes: get(p,name)
    seven_prime=data['separation_prime']
    for e in (0,1):
        get(seven_prime,label(universe,e))
    if {(int(p),n) for p,t in rows.items() for n in t}!=set(used):
        raise ArithmeticError('Unused witness polynomials')
    with ProcessPoolExecutor(max_workers=workers) as pool:
        results=dict(zip(used,pool.map(inspect,used.values(),chunksize=1)))
        for name,primes in data['simplicity'].items():
            S,e=parse(name);d=len(S)*(1<<(len(S)-1));remain=set(range(1,1+d//2))
            for p in primes:
                options=results[p,name][0]
                if options is None: raise ArithmeticError('Nonordinary dimension witness')
                remain.intersection_update(options)
            if remain: raise ArithmeticError(f'Proper dimensions survive for {name}: {sorted(remain)}')
        if results[seven_prime,label(universe,0)][1]==results[seven_prime,label(universe,1)][1]:
            raise ArithmeticError('The two seven-index factors were not separated')
        groups={}
        for (p,name),(_,S,e,coeff) in used.items():groups.setdefault((p,S),{})[e]=coeff
        jobs=[(p,S,vals) for (p,S),vals in groups.items()]
        if not recheck:
            jobs=[next(j for j in jobs if j[0]==seven_prime and len(j[1])==7)]
        list(pool.map(recount,jobs,chunksize=1))
        # The all-stage-ten 193 receipt contains these fourteen first-seven factors.
        external=json.loads((SIEVE/'six-index.json').read_text())
        p=external['prime']; six_jobs=[]
        for S in six:
            check_stage(p,S)
            for e in (0,1):
                six_jobs.append((p,S,e,external['unit_polynomials'][label(S,e)]))
        signatures={tuple(fp) for opt,fp in pool.map(inspect,six_jobs)}
        if len(signatures)!=14:raise ArithmeticError('Some first-seven six-index factors coincide')
    for n in (1,2):
        literal=trace_sums(seven_prime,universe,n,len(universe))[1][universe]
        for e,t in enumerate(literal):
            if (trace_from_factor(get(seven_prime,label(universe,e)),n)-t)%seven_prime:
                raise ArithmeticError('Direct first-seven signed point count failed')
    print(f'Verified all 14 additional geometric simplicity proofs and the 255-factor seventh stage. '
          f'Recomputed {len(jobs)}/{len(groups)} prime/support matrices.')


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--recount',action='store_true')
    parser.add_argument('--workers',type=int,default=1)
    args=parser.parse_args()
    run(json.loads((HERE/'certificate.json').read_text()),args.recount,args.workers)
