"""Find and replay geometric certificates for every five-index support at stage ten."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys

import sympy as sp

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"cubic-tower-geometric-separation"))
sys.path.insert(0,str(HERE.parent/"cubic-tower-factors"))
from affine_frobenius import signed_numerators
from tower import check_stage
from unit_orbits import observations

SUPPORTS=list(combinations(range(10),5))
DIM=80


@lru_cache(None)
def excluded_prime_products():
    z=sp.Symbol("z")
    f=-4*z**3+z**2+z+1
    h=[(12*i+1)*z**2+(12*i*i+1)*z+4*i**3+1 for i in range(10)]
    singles=[int(6*sp.discriminant(f,z)*(12*i+1)*sp.discriminant(h[i],z)*sp.resultant(f,h[i],z))
             for i in range(10)]
    pairs={(i,j):int(sp.resultant(h[i],h[j],z)) for i,j in combinations(range(10),2)}
    return {S:sp.prod(singles[i] for i in S)*sp.prod(pairs[e] for e in combinations(S,2)) for S in SUPPORTS}


def evaluate(job):
    p,index,signs=job
    polys=signed_numerators(p,SUPPORTS[index])
    return index, {sign:(polys[sign],*observations(polys[sign],p,DIM)) for sign in signs}


def name(i,sign):
    return ",".join(map(str,SUPPORTS[i]))+("/+" if sign==0 else "/-")


def build(workers,output):
    possible={name(i,e):set(range(1,DIM//2+1)) for i in range(len(SUPPORTS)) for e in (0,1)}
    groups=[sorted(possible)]
    receipt={"dimension":DIM,"simplicity":{key:[] for key in possible},"separation":[],"observations":{}}
    goods=excluded_prime_products()
    with ProcessPoolExecutor(max_workers=workers) as pool:
        for p in sp.primerange(DIM+1,4000):
            need_separate={key for g in groups if len(g)>1 for key in g}
            compare=need_separate if all(v%p for v in goods.values()) else set()
            jobs=[]
            for i,S in enumerate(SUPPORTS):
                signs=tuple(e for e in (0,1) if name(i,e) in compare or possible[name(i,e)])
                if signs and goods[S]%p: jobs.append((p,i,signs))
            fp={}
            for i,res in pool.map(evaluate,jobs,chunksize=4):
                for e,(coeff,allowed,fingerprint) in res.items():
                    key=name(i,e)
                    store=False
                    if key in compare:
                        fp[key]=tuple(fingerprint)
                        store=True
                    if allowed is not None and possible[key]:
                        new=possible[key].intersection(allowed)
                        if new!=possible[key]:
                            receipt["simplicity"][key].append(p)
                            possible[key]=new
                            store=True
                    if store:
                        receipt["observations"].setdefault(str(p),{})[key]=coeff
            if compare:
                refined=[]
                for grp in groups:
                    if len(grp)==1:
                        refined.append(grp)
                        continue
                    buckets={}
                    for key in grp:
                        buckets.setdefault(fp[key],[]).append(key)
                    refined.extend(buckets.values())
                groups=refined
                receipt["separation"].append(p)
            unresolved=sum(bool(v) for v in possible.values())
            print(p,len(jobs),unresolved,sum(len(g)*(len(g)-1)//2 for g in groups),flush=True)
            receipt["unresolved"]={key:sorted(s) for key,s in possible.items() if s}
            tmp=output.with_suffix(".tmp")
            tmp.write_text(json.dumps(receipt,separators=(",",":"))+"\n")
            tmp.replace(output)
            if unresolved==0 and all(len(g)==1 for g in groups):
                return
    raise RuntimeError("Some classes still need witnesses")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers",type=int,default=2)
    parser.add_argument("--output",type=Path,required=True)
    args=parser.parse_args()
    build(args.workers,args.output)
