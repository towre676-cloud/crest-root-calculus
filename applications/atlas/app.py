"""Compute an elliptic-quotient certificate atlas for smooth plane quartics."""
import argparse
from fractions import Fraction
import json
from itertools import product
from pathlib import Path
import sys
from algebra import (characteristic, difference, gcd, mons, parse_form, point_count,
                     simplicity_witness, smoothness, trim)

ROOT = Path(__file__).resolve().parent


def poly_mul(a,b):
    out = [0]*(len(a)+len(b)-1)
    for i,v in enumerate(a):
        for j,w in enumerate(b):
            out[i+j] += v*w
    return out


def poly_sum(*terms):
    out = [0]*max(len(a) for _,a in terms)
    for scalar,a in terms:
        for i,v in enumerate(a):
            out[i] += scalar*v
    return trim(out)


def evaluate(form, xyz):
    x,y,z=xyz
    return sum(v*x**a*y**b*z**c for (a,b,c),v in form.items())


def cyclic_map(form):
    if any(v!=form.get((b,c,a),0) for (a,b,c),v in form.items()):
        return None
    def sym(a,b,c):
        return form.get((a,b,c),0)+form.get((b,a,c),0)
    A=sym(4,0,0)
    B=sym(3,1,0)-4*A
    C=sym(2,2,0)-6*A-2*B
    D=sym(2,1,1)-12*A-5*B-2*C
    K=form.get((1,3,0),0)-form.get((3,1,0),0)
    if D==0:
        return None
    # Independent polynomial-identity test. Bounds on each variable are 4,
    # so the 5^3 grid proves the identity, rather than merely sampling it.
    for x,y,z in product(range(5), repeat=3):
        s=x+y+z
        t=x*y+y*z+z*x
        u=x*y*z
        v=(x-y)*(y-z)*(z-x)
        if 2*evaluate(form,(x,y,z)) != A*s**4+B*s*s*t+C*t*t+D*s*u+K*s*v:
            raise ArithmeticError("invariant identity failed")
    N=[A,B,C]
    L=D*D+27*K*K
    S=poly_sum((K,[4*D,-18*D]), (-54*K,N))
    T=poly_sum((D*D,[0,0,1,-4]),(1,poly_mul([4*D,-18*D],N)),(-27,poly_mul(N,N)))
    target=poly_sum((1,poly_mul(S,S)),(4*L,T))
    if len(target) not in (4,5):
        return None
    for ell in (5,7,11,13,17,19,23,29,31,37,41,43):
        mod = trim([a%ell for a in target])
        deriv = trim([(i*a)%ell for i,a in enumerate(target) if i])
        if len(mod)==len(target) and len(gcd(mod,deriv,ell))==1:
            break
    else:
        return None
    origin=None
    for x,y,z in product([0,1,-1,2,-2,3,-3],repeat=3):
        s=x+y+z
        if not s or evaluate(form,(x,y,z)):
            continue
        U=Fraction(x*y+y*z+z*x,s*s)
        V=Fraction((x-y)*(y-z)*(z-x),s**3)
        Y=2*L*V-sum(v*U**j for j,v in enumerate(S))
        if Y*Y!=sum(v*U**j for j,v in enumerate(target)):
            raise ArithmeticError("point failed target equation")
        origin={"source":[x,y,z],"target":[str(U),str(Y)]}
        break
    return {"degree":3,"field":"Q","invariants_A_B_C_D_K":[A,B,C,D,K],"L":L,
            "S_coefficients_ascending":S,"target_Y_squared_coefficients_ascending":target,
            "coordinates":"s=x+y+z; U=(xy+yz+zx)/s^2; V=(x-y)(y-z)(z-x)/s^3; Y=2*L*V-S(U)",
            "squarefree_target_mod":ell,"rational_origin":origin}


def classify(entry, primes):
    if not isinstance(entry,dict) or not isinstance(entry.get("id"),str):
        raise ValueError("each curve needs id and quartic")
    form=parse_form(entry.get("quartic"))
    results=[]
    for p in primes:
        smooth=smoothness(form,p)
        if not smooth:
            results.append({"p":p,"status":"smoothness_not_certified"})
            continue
        quotient=cyclic_map(form)
        if quotient:
            status="degree_three_elliptic_map_over_Q" if quotient["rational_origin"] else "degree_three_genus_one_map_over_Q"
            return {"id":entry["id"],"status":status,"smoothness":{"p":p,**smooth},"map":quotient}
        counts=[point_count(form,p,n) for n in (1,2,3)]
        char=characteristic(p,counts)
        witness=simplicity_witness(char)
        results.append({"p":p,"smoothness":smooth,"point_counts":counts,
                        "frobenius_ascending":char,"simplicity_certificate":witness})
        if witness:
            return {"id":entry["id"],"status":"no_elliptic_map_over_Qbar","reductions":results}
    return {"id":entry["id"],"status":"unresolved","reductions":results}


def solve(data):
    if isinstance(data,dict) and data.get("engine")=="part_xxx":
        from original.runner import solve as run_original
        return run_original(data)
    if not isinstance(data,dict): raise ValueError("input must be a JSON object")
    entries=data.get("curves")
    primes=data.get("primes",[2,3,5])
    if not isinstance(entries,list) or not 1<=len(entries)<=32:
        raise ValueError("provide 1–32 curves")
    if not isinstance(primes,list) or not primes or len(set(primes))!=len(primes) or any(type(p) is not int or p not in (2,3,5) for p in primes):
        raise ValueError("primes must list distinct choices among 2, 3, 5")
    if len({e.get("id") for e in entries if isinstance(e,dict)}) != len(entries):
        raise ValueError("curve ids must be distinct strings")
    results=[classify(entry,primes) for entry in entries]
    return {"counts":{s:sum(r["status"]==s for r in results) for s in sorted({r["status"] for r in results})},"curves":results}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",nargs="?",type=Path,default=ROOT/"example.json")
    args=parser.parse_args()
    try: print(json.dumps(solve(json.loads(args.input.read_text(encoding="utf-8"))),indent=2))
    except (OSError,ValueError,KeyError,TypeError) as error: parser.exit(2,f"Input error: {error}\n")
