"""Reconstruct finite one-dimensional site patterns and select a local observable."""
import argparse
from bisect import bisect_right
from collections import Counter
from itertools import combinations
import json
from math import comb, isqrt
from pathlib import Path

ROOT=Path(__file__).resolve().parent


def distances(points):
    return sorted(b-a for a,b in combinations(points,2))


def canonical(points):
    forward=tuple(sorted(points))
    reverse=tuple(forward[-1]-x for x in reversed(forward))
    return min(forward,reverse)


def local_count(points, k, width):
    """Count k-site subsets whose complete span is at most width."""
    total=0
    for i,x in enumerate(points):
        others=bisect_right(points,x+width)-i-1
        if others>=k-1:
            total+=comb(others,k-1)
    return total


def integer(value, label, low, high):
    if type(value) is not int or not low<=value<=high:
        raise ValueError(f"{label} must be an integer in [{low}, {high}]")
    return value


def parse_probe(probe,n,observed):
    if not isinstance(probe,dict):
        raise ValueError("each probe must be an object")
    k=integer(probe.get("sites"),"sites",2,n)
    width=integer(probe.get("span"),"span",0,4096)
    key="count" if observed else "cost"
    val=integer(probe.get(key,1 if not observed else None),key,0,10**9)
    return k,width,val


class BudgetReached(Exception):
    pass


def reconstruct(y,n,e,budget):
    seen=set()
    result=set()
    visited=0
    def visit():
        nonlocal visited
        if visited>=budget or len(result)>=1024:
            raise BudgetReached
        visited+=1
    def exact(points,remaining):
        # Remaining distances are uniquely determined by the placed point set.
        key=canonical(points)
        if key in seen:
            return
        seen.add(key)
        visit()
        if not remaining:
            if len(points)==n:
                result.add(key)
            return
        far=max(remaining)
        for x in sorted({far,y[-1]-far}):
            if x in points:
                continue
            used=Counter(abs(x-v) for v in points)
            if all(remaining[d]>=k for d,k in used.items()):
                exact((*points,x),remaining-used)
    def fits_partial(measured):
        # Greedy earliest available match is optimal for equal-width intervals.
        j=0
        for d in measured:
            while j<len(y) and y[j]<d-e:
                j+=1
            if j==len(y) or y[j]>d+e:
                return False
            j+=1
        return True
    def approximate(interior,span):
        visit()
        points=(0,)+interior+(span,)
        if not fits_partial(distances(points)):
            return
        if len(points)==n:
            result.add(canonical(points))
            return
        next_x=interior[-1]+1 if interior else 1
        for x in range(next_x,span-(n-len(points))+1):
            approximate(interior+(x,),span)
    complete=True
    try:
        if e==0:
            remaining=Counter(y)
            remaining.subtract([y[-1]])
            exact((0,y[-1]),+remaining)
        else:
            for span in range(max(n-1,y[-1]-e),y[-1]+e+1):
                approximate((),span)
    except BudgetReached:
        complete=False
    return sorted(result),visited,complete


def solve(data):
    if not isinstance(data,dict):
        raise ValueError("input must be an object")
    if data.get("mode")=="cubic_powder":
        from powder import solve as index_powder
        return index_powder(data)
    raw=data.get("distances")
    if not isinstance(raw,list) or not 1<=len(raw)<=190:
        raise ValueError("supply all pair distances for 2 to 20 distinct sites")
    y=sorted(integer(v,"distance",0,2048) for v in raw)
    n=(1+isqrt(1+8*len(y)))//2
    if comb(n,2)!=len(y):
        raise ValueError("distance count must equal n*(n-1)/2")
    e=integer(data.get("distance_error",0),"distance_error",0,64)
    if not e and y[0]==0:
        raise ValueError("distinct sites have positive exact distances")
    budget=integer(data.get("node_budget",50000),"node_budget",1,250000)
    obs=data.get("observations",[])
    if not isinstance(obs,list) or len(obs)>200:
        raise ValueError("use a list of up to 200 observations")
    observations=[parse_probe(p,n,True) for p in obs]
    probes=data.get("probes")
    if probes is not None:
        if not isinstance(probes,list) or len(probes)>200:
            raise ValueError("use a list of up to 200 available probes")
        probes=[parse_probe(p,n,False) for p in probes]
    models,nodes,complete=reconstruct(y,n,e,budget)
    before=len(models)
    models=[p for p in models if all(local_count(p,k,w)==v for k,w,v in observations)]
    if not complete:
        status="search_incomplete"
    elif not models:
        status="incompatible"
    elif len(models)==1:
        status="unique_up_to_translation_and_reflection"
    else:
        status="ambiguous"
    best=None
    if status=="ambiguous":
        if probes is None:
            widths=sorted({d for p in models for d in distances(p)})
            probes=[(k,w,1) for k in range(3,min(5,n)+1) for w in widths]
        ranking=None
        for k,w,cost in probes:
            values=[local_count(p,k,w) for p in models]
            bins=Counter(values)
            if len(bins)<=1:
                continue
            score=(max(bins.values()),cost,sum(t*t for t in bins.values()),k,w)
            if ranking is None or score<ranking:
                ranking=score
                best={"sites":k,"span":w,"cost":cost,"worst_remaining":score[0],
                      "predicted_count_per_candidate":values,
                      "question":f"How many groups of {k} sites span at most {w} grid units?"}
    return {"status":status,"complete":complete,"site_count":n,"nodes_visited":nodes,
            "classes_before_observations":before,"candidates":[list(p) for p in models],
            "next_probe":best,"model":"distinct 1-D integer sites, equal weights, all pair separations measured"}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",nargs="?",type=Path,default=ROOT/"example.json")
    args=parser.parse_args()
    try: print(json.dumps(solve(json.loads(args.input.read_text(encoding="utf-8"))),indent=2))
    except (ValueError,OSError) as error: parser.exit(2,f"Input error: {error}\n")
