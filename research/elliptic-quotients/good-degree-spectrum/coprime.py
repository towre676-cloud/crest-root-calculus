"""Certify separation using only degrees coprime to the theta level."""
import argparse
from collections import defaultdict
import json
from math import gcd
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "inverse-spectrum"))
from census import CLASS_NUMBER_ONE, residue_field, scalar_orbits, square_shells


def primitive_shells(field, bound):
    return [[(u,x,y) for u,x,y in shell if (u,x,y) > (0,0,0) and gcd(u,x,y)==1]
            for shell in square_shells(field,bound)]


def counts_on_shell(centers, traces, shell, p):
    """Count using affine lines in the trace coordinates or direct evaluation."""
    if len(centers) <= p:
        return [2*sum((u+s*x+t*y)%p==0 for u,x,y in shell) for s,t in (traces[a] for a in centers)]
    counts = defaultdict(int)
    lines = defaultdict(int)
    for u,x,y in shell:
        u,x,y = u%p,x%p,y%p
        if y:
            inv = pow(y,-1,p)
            lines[(1,-x*inv%p,-u*inv%p)] += 2
        elif x:
            lines[(0,0,-u*pow(x,-1,p)%p)] += 2
        # If x=y=0 modulo p, a degree prime to p forces u nonzero.
    for (vary_s,slope,offset),weight in lines.items():
        if vary_s:
            for s in range(p):
                counts[s,(slope*s+offset)%p] += weight
        else:
            for t in range(p):
                counts[offset,t] += weight
    return [counts[traces[a]] for a in centers]


def certify(D, field, shells):
    p = field.p
    centers = list(scalar_orbits(field))
    traces = {a:(field.trace(a),field.trace(field.mul(a,(0,1)))) for a in centers}
    groups = [centers] if len(centers)>1 else []
    last = None
    for degree,shell in enumerate(shells,1):
        if not groups:
            break
        if gcd(degree,2*D*p)!=1:
            continue
        flat = [a for group in groups for a in group]
        counts = dict(zip(flat,counts_on_shell(flat,traces,shell,p)))
        refined = []
        for group in groups:
            by_count = defaultdict(list)
            for a in group:
                by_count[counts[a]].append(a)
            if len(by_count)>1:
                lower,higher = sorted(by_count)[:2]
                last = {"degree":degree,"centers":[by_count[lower][0],by_count[higher][0]],"counts":[lower,higher]}
            refined.extend(g for g in by_count.values() if len(g)>1)
        groups = refined
    return {"D":D,"ell":p,"forms":len(centers),"last_witness":last,"unresolved":groups}


def compare_receipt(output, receipt):
    rendered = json.dumps(output, indent=2)
    if receipt is not None and json.loads(rendered) != json.loads(receipt.read_text()):
        raise ArithmeticError("The recomputed census differs from the receipt.")
    print(rendered)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-limit",type=int,default=199)
    parser.add_argument("--degree-limit",type=int,default=512)
    parser.add_argument("--discriminant",type=int,choices=CLASS_NUMBER_ONE,action="append")
    parser.add_argument("--verify",type=Path,help="Compare with the saved receipt.")
    args=parser.parse_args()
    if args.prime_limit<3 or args.degree_limit<1:
        parser.error("Use a gluing prime bound of at least 3 and a positive degree limit.")
    output=[]
    for D in args.discriminant or CLASS_NUMBER_ONE:
        shells=None
        for p in range(3,args.prime_limit+1,2):
            try:
                field=residue_field(D,p)
            except ValueError:
                continue
            if shells is None:
                shells=primitive_shells(field,args.degree_limit)
            result=certify(D,field,shells)
            output.append(result)
            witness=result["last_witness"]
            print(D,p,result["forms"],witness["degree"] if witness else None,len(result["unresolved"]),
                  file=sys.stderr,flush=True)
    compare_receipt(output,args.verify)


if __name__ == "__main__":
    main()
