"""Certify the index and rational map degrees of a smooth plane quartic."""
import argparse
from itertools import combinations
import json
from math import gcd, isqrt, prod
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "involution-atlas"))
from classify import classify, parse_quartic

XYZ = sp.symbols("X Y Z")
ODD_ISOGENY_PRIMES = (3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163)


def expr(form):
    return sum(c*prod(x**i for x, i in zip(XYZ, m)) for m,c in form.items())


def integral_form(f):
    poly = sp.Poly(f, *XYZ, domain=sp.QQ).clear_denoms(convert=True)[1].primitive()[1]
    return {m:int(c) for m,c in poly.terms()}


def normalized(form):
    return {m:c//gcd(*form.values()) for m,c in form.items()}


def is_even(form):
    return all(all(i%2 == 0 for i in m) for m in form)


def smooth(form):
    f = expr(form)
    derivatives = [f.diff(x) for x in XYZ]
    for i, x in enumerate(XYZ):
        rest = XYZ[:i]+XYZ[i+1:]
        if list(sp.groebner([g.subs(x,1) for g in derivatives], *rest,domain=sp.QQ)) != [1]:
            raise ValueError("The source is singular over the algebraic closure")


def valuation(a,p):
    n=0
    while a%p == 0:
        a//=p
        n+=1
    return n,a


def pure(form):
    return [form.get(tuple(4 if k==i else 0 for k in range(3)),0) for i in range(3)]


def index_witness(form):
    diagonal=pure(form)
    if not all(diagonal):
        return None
    for p in sorted(sp.factorint(abs(prod(diagonal)))):
        if p==2:
            continue
        factors=[valuation(v,p) for v in diagonal]
        for shift in range(4):
            scales=[(v+shift)//4 for v,u in factors]
            vals=[(v+shift)%4 for v,u in factors]
            if len(set(vals)) != 3:
                continue
            mixed=[valuation(c,p)[0]+shift-sum(i*k for i,k in zip(m,scales))
                   for m,c in form.items() if 4 not in m]
            if mixed and min(mixed) <= max(vals):
                continue
            i,j=next((i,j) for i,j in combinations(range(3),2) if vals[i]%2==vals[j]%2)
            residue=-factors[i][1]*pow(factors[j][1],-1,p)%p
            if pow(residue,(p-1)//2,p) == p-1:
                return {"prime":int(p),"global_valuation_shift":shift,"coordinate_scales":scales,
                        "diagonal_valuations":vals,"competing_indices":[i,j],
                        "nonsquare_unit_ratio":int(residue),
                        "least_mixed_valuation":min(mixed) if mixed else None}
    return None


def even_frame(form):
    if is_even(form):
        return form,sp.eye(3)
    f=expr(form)
    factors=sp.factor_list(sp.hessian(f,XYZ).det())[1]
    if len(factors)==3 and all(sp.Poly(h,*XYZ).total_degree()==1 and m==2 for h,m in factors):
        axes=sp.Matrix([[h.diff(x) for x in XYZ] for h,m in factors])
        if axes.det():
            change=axes.inv()
            recovered=integral_form(f.subs(dict(zip(XYZ,change*sp.Matrix(XYZ))),simultaneous=True))
            if is_even(recovered):
                return recovered,change
    atlas=classify({"terms":[[*m,c] for m,c in form.items()]})
    centers=[sp.Matrix([sp.Rational(x) for x in r["center"]]) for r in atlas["orbits"]
             if r["orbit_size"]==1]
    for triple in combinations(centers,3):
        change=sp.Matrix.hstack(*triple)
        if not change.det():
            continue
        recovered=integral_form(f.subs(dict(zip(XYZ,change*sp.Matrix(XYZ))),simultaneous=True))
        if is_even(recovered):
            return recovered,change
    return None,None


def sign_models(form):
    # On the i-th quotient, w = 2a_i X_i^2 + m_ij X_j^2 + m_ik X_k^2.
    diag=pure(form)
    mixed=lambda i,j: form.get(tuple(2 if k in (i,j) else 0 for k in range(3)),0)
    rows=[]
    for i in range(3):
        j,k=[v for v in range(3) if v!=i]
        a=diag[i]
        b,c=mixed(i,j),mixed(i,k)
        A=b*b-4*a*diag[j]
        B=2*b*c-4*a*mixed(j,k)
        C=c*c-4*a*diag[k]
        if not A*C*(B*B-4*A*C):
            raise ArithmeticError("Repeated branch point on the sign quotient")
        rows.append({"axis":i,"binary_quartic":[A,B,C],"w_coefficients":[2*a,b,c],
                     "jacobian":[-2*B,B*B-4*A*C],"degree_four_target":[B,A*C]})
    return rows


def curve_trace(model,p):
    a,b=model
    if (2*b*(a*a-4*b))%p==0:
        return None
    squares={x*x%p for x in range(1,p)}
    return -sum(0 if k==0 else 1 if k in squares else -1
                for k in (x*(x*x+a*x+b)%p for x in range(p)))


def nonsquare(n):
    return n<0 or isqrt(n)**2 != n


def degree_sieve(models,bound):
    remaining=[set(ODD_ISOGENY_PRIMES) for m in models]
    separation={}
    traces={}
    for p in sp.primerange(3,bound+1):
        values=[curve_trace(m["jacobian"],p) for m in models]
        used=False
        for i,t in enumerate(values):
            if t is None:
                continue
            excluded={l for l in remaining[i] if l!=p and pow(t*t-4*p,(l-1)//2,l)==l-1}
            if excluded:
                used=True
                remaining[i]-=excluded
            for j in range(i):
                if values[j] is not None and t!=values[j] and (j,i) not in separation:
                    separation[j,i]=p
                    used=True
        if used:
            traces[p]=values
        if len(separation)==3 and not any(remaining):
            break
    two_chains=[]
    for model in models:
        two_chains.append(all(nonsquare(a*a-4*b) for a,b in
                              (model["jacobian"],model["degree_four_target"])))
    return {"frobenius_traces":traces,"pair_witness_primes":{f"{i},{j}":p for (i,j),p in separation.items()},
            "surviving_odd_prime_degrees":[sorted(s) for s in remaining],
            "unique_rational_two_kernel_at_both_ends":two_chains}


def cubic_orbit(form):
    # A unit determinant over Z/4 certifies independent fourth-root classes over Q(i).
    if any(4 not in m for m in form):
        return None
    a,b,c=pure(form)
    primes=[p for p in sp.factorint(abs(a*b*c)) if p!=2]
    rows={p:[valuation(u,p)[0]-valuation(a,p)[0] for u in (b,c)] for p in primes}
    for p,q in combinations(primes,2):
        u,v=rows[p]; w,z=rows[q]
        if (u*z-v*w)%2:
            return {"odd_primes":[int(p),int(q)],"valuation_rows":[rows[p],rows[q]],
                    "geometric_degree_three_classes":16,"galois_orbit_size":16,
                    "least_degree_of_a_direction_field":16}
    return None


def analyze(data,bound=113,recover=False):
    form=normalized(parse_quartic(data))
    smooth(form)
    out={"input_index_witness":index_witness(form)}
    if is_even(form) or recover:
        even,change=even_frame(form)
    else:
        even,change=None,None
    wit=out["input_index_witness"]
    if even:
        models=sign_models(even)
        w=index_witness(even)
        wit=wit or w
        out.update({"old_coordinates_in_new_frame":[[str(a) for a in row] for row in change.tolist()],
                    "even_frame_terms":[[*m,c] for m,c in sorted(even.items(),reverse=True)],
                    "frame_index_witness":w,"sign_quotients":models,
                    "fermat_cubic_direction_orbit":cubic_orbit(even)})
        sieve=degree_sieve(models,bound)
        out["isogeny_sieve"]=sieve
        if wit and len(sieve["pair_witness_primes"])==3:
            out["rational_optimal_genus_one_degrees"]=[2,2,2]
            if not any(sieve["surviving_odd_prime_degrees"]) and all(sieve["unique_rational_two_kernel_at_both_ends"]):
                out["all_rational_genus_one_degrees"]=["2*n^2","4*n^2"]
                out["all_rational_elliptic_degrees"]=["4*n^2","8*n^2"]
        if wit:
            out["least_rational_genus_one_degree"]=2
            out["least_rational_elliptic_degree"]=4
    if wit:
        out["source_index"]=4
        out["rational_genus_one_map_degrees_are_even"]=True
        out["rational_elliptic_map_degrees_are_divisible_by_four"]=True
    return out


if __name__=="__main__":
    p=argparse.ArgumentParser(description=__doc__)
    p.add_argument("input",type=Path)
    p.add_argument("--prime-bound",type=int,default=113)
    p.add_argument("--recover-frame",action="store_true",help="recover rational sign axes in mixed coordinates")
    args=p.parse_args()
    print(json.dumps(analyze(json.loads(args.input.read_text()),args.prime_bound,args.recover_frame),indent=2))
