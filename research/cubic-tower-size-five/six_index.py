"""Separate the six-index signed factors of stage ten using powered unit-root traces."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from itertools import combinations
import json
import os
from pathlib import Path

if __name__=="__main__":
    os.environ["OPENBLAS_NUM_THREADS"]="1"
    os.environ["MKL_NUM_THREADS"]="1"

from five import check_stage, signed_numerators
from block_frobenius import signed_blocks
from unit_orbits import observations

HERE=Path(__file__).resolve().parent
SUPPORTS=list(combinations(range(10),6))
PRIME=193
DIM=192


def label(S,e):
    return ",".join(map(str,S))+("/+" if e==0 else "/-")


def local(S):
    check_stage(PRIME,S)
    return {label(S,e):c for e,c in enumerate(signed_numerators(PRIME,S))}


def blocked(S):
    check_stage(PRIME,S)
    return {label(S,e):c for e,c in enumerate(signed_blocks(PRIME,S))}


def inspect(pair):
    name,poly=pair
    if len(poly)!=DIM+1 or poly[0]!=1 or any(type(c)!=int or not 0<=c<PRIME for c in poly):
        raise ArithmeticError("Noncanonical unit polynomial")
    allowed,signature=observations(poly,PRIME,DIM)
    return name,signature,allowed==[]


def run(path,build,recount,workers):
    if build:
        coeffs={}
    else:
        data=json.loads(path.read_text())
        if data["prime"]!=PRIME or data["dimension"]!=DIM:
            raise ArithmeticError("Unexpected six-index certificate parameters")
        coeffs=data["unit_polynomials"]
    for S in SUPPORTS:
        check_stage(PRIME,S)
    with ProcessPoolExecutor(max_workers=workers) as pool:
        jobs=SUPPORTS if build or recount else [SUPPORTS[0]]
        for group in pool.map(local if build else blocked,jobs,chunksize=2):
            for key,row in group.items():
                if build: coeffs[key]=row
                elif row!=coeffs[key]: raise ArithmeticError("Local polynomial differs for "+key)
        if set(coeffs)!={label(S,e) for S in SUPPORTS for e in (0,1)}:
            raise ArithmeticError("Incomplete six-index support coverage")
        groups={}
        simple=[]
        for name,signature,proved in pool.map(inspect,coeffs.items(),chunksize=3):
            groups.setdefault(tuple(signature),[]).append(name)
            if proved: simple.append(name)
    collisions=sorted(sorted(g) for g in groups.values() if len(g)>1)
    summary={"separated_pairs":420*419//2-sum(len(g)*(len(g)-1)//2 for g in collisions),
             "unseparated_groups":collisions,"simple_at_this_prime":sorted(simple)}
    if build:
        path.write_text(json.dumps({"prime":PRIME,"dimension":DIM,**summary,"unit_polynomials":coeffs},
                                    separators=(",",":"))+"\n")
    elif any(summary[k]!=data[k] for k in summary):
        raise ArithmeticError("Six-index geometric result differs")
    print(json.dumps(summary,indent=2))


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--build",action="store_true")
    parser.add_argument("--recount",action="store_true")
    parser.add_argument("--workers",type=int,default=1)
    parser.add_argument("--file",type=Path,default=HERE/"six-index.json")
    args=parser.parse_args()
    run(args.file,args.build,args.recount,args.workers)
