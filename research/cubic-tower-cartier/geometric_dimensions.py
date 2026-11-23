"""Certify geometric simplicity by intersecting ordinary dimension restrictions."""
import argparse
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.galoistools import gf_factor, gf_pow_mod
from sympy.polys.domains import ZZ

from cartier import Cartier, numerators

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'cubic-tower-factors'))
from tower import check_stage


@lru_cache(None)
def root_of_unity_exponent(degree):
    """LCM of all n with phi(n) <= degree^2."""
    bound=degree**2
    result=1
    for prime in sp.primerange(2,bound+2):
        power=prime
        while power*(prime-1)<=bound:
            power*=prime
        result*=power
    return result


@lru_cache(None)
def subdegrees(p,coeffs):
    if not coeffs[-1]:
        return None
    dim=len(coeffs)-1
    M=root_of_unity_exponent(2*dim)
    _,factors=gf_factor(list(coeffs),p,ZZ)
    orbit_sizes=[]
    for f,mult in factors:
        n=len(f)-1
        h=p**n-1
        remaining=n
        for ell in sp.factorint(n):
            while remaining%ell==0:
                d=remaining//ell
                exponent=(M*(p**d-1))%h
                if gf_pow_mod([1,0],exponent,f,p,ZZ)!=[1]:
                    break
                remaining=d
        orbit_sizes.extend([remaining]*(mult*n//remaining))
    allowed=1
    for d in orbit_sizes:
        allowed |= allowed<<d
    return tuple(i for i in range(1,dim//2+1) if (allowed>>i)&1)


def collect(max_prime=251):
    supports=list(combinations(range(10),3))
    remaining={(s,e):set(range(1,7)) for s in supports for e in (0,1)}
    records={','.join(map(str,s))+('/-' if e else '/+'):[] for s,e in remaining}
    count=0
    for p in sp.primerange(11,max_prime+1):
        if not any(remaining.values()):
            break
        engine=Cartier(p,3)
        processed=0
        for s in supports:
            if not (remaining[s,0] or remaining[s,1]):
                continue
            try:
                check_stage(p,s)
            except ValueError:
                continue
            coeffs=numerators(engine.matrices(s),p)
            processed+=1
            for e,co in enumerate(coeffs):
                previous=remaining[s,e]
                if not previous:continue
                options=subdegrees(p,tuple(co))
                if options is None:continue
                reduced=previous.intersection(options)
                if reduced!=previous:
                    remaining[s,e]=reduced
                    records[','.join(map(str,s))+('/-' if e else '/+')].append(
                        {'p':p,'unit_polynomial':co,'permitted_half_dimensions':list(options)})
        count+=processed
        print(p,processed, sum(not x for x in remaining.values()),count,flush=True,file=sys.stderr)
    return records


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--search',action='store_true')
    parser.add_argument('--max-prime',type=int,default=251)
    options=parser.parse_args()
    if options.search:
        print(json.dumps(collect(options.max_prime),indent=2))
