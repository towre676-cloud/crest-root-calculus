"""Build the three projective-map incidence charts for a smooth rational plane quartic."""
import argparse
from collections import defaultdict
from functools import lru_cache
from itertools import combinations
import json
from pathlib import Path
import sys

import sympy as sp

U,V,W = sp.symbols("u v w")


def powers(d):
    return [(i,j,d-i-j) for i in range(d+1) for j in range(d-i+1)]


def source(data):
    terms = data["terms"] if isinstance(data, dict) else data
    F = defaultdict(lambda: sp.S.Zero)
    for i,j,k,c in terms:
        if any(type(m) != int or m < 0 for m in (i,j,k)) or i+j+k != 4:
            raise ValueError("Expected homogeneous quartic terms")
        F[i,j,k] += sp.Rational(c)
    F = {m:c for m,c in F.items() if c}
    f = evaluate(F, (U,V,W))
    for t in (U,V,W):
        others = [s for s in (U,V,W) if s != t]
        if list(sp.groebner([f.diff(s).subs(t,1) for s in (U,V,W)], *others, domain=sp.QQ)) != [1]:
            raise ValueError("Expected a smooth quartic")
    return F


def evaluate(form, point):
    x,y,z = point
    return sum(c*x**i*y**j*z**k for (i,j,k),c in form.items())


def add(*forms):
    out = defaultdict(list)
    for f in forms:
        for m,c in f.items():
            out[m].append(c)
    return {m:sum(terms) for m,terms in out.items()}


def scale(f, scalar):
    return {m:scalar*c for m,c in f.items()}


def mul(f,g):
    out = defaultdict(list)
    for (i,j,k),c in f.items():
        for (p,q,r),d in g.items():
            out[i+p,j+q,k+r].append(c*d)
    return {m:sum(terms) for m,terms in out.items()}


def deriv(f, index):
    out = {}
    for m,c in f.items():
        if m[index]:
            n = list(m)
            n[index] -= 1
            out[tuple(n)] = m[index]*c
    return out


def reducer(F):
    pivot = max(F)
    @lru_cache(None)
    def monomial(m):
        if any(i < j for i,j in zip(m,pivot)):
            return {m:sp.S.One}
        shift = tuple(i-j for i,j in zip(m,pivot))
        return add(*(scale(monomial(tuple(a+b for a,b in zip(n,shift))), -c/F[pivot])
                     for n,c in F.items() if n != pivot))
    def reduce_form(form):
        # Group whole coefficient expressions before SymPy combines them.
        return add(*(scale(monomial(m),c) for m,c in form.items()))
    return reduce_form


def identities(F, forms, a, b, ell, lam):
    A,B,C = forms
    red = reducer(F)
    def D(p):
        return add(mul(deriv(F,1),deriv(p,0)),scale(mul(deriv(F,0),deriv(p,1)),-1))
    AA, BB, CC = mul(A,A), mul(B,B), mul(C,C)
    cubic = add(mul(BB,C),scale(mul(AA,A),-1),scale(mul(A,CC),-a),scale(mul(CC,C),-b))
    W_ell = {(1,0,1):ell[0], (0,1,1):ell[1], (0,0,2):ell[2]}
    diff = add(mul(C,D(A)),scale(mul(A,D(C)),-1),scale(mul(mul(B,C),W_ell),-2*lam))
    return [*red(cubic).values(), *red(diff).values()]


def axis_unit(orbit, ell, separator):
    # Every conjugate axis must differ from ell. A single generic separator works
    # for all conjugates and orbits, over an infinite field.
    T = sp.Symbol("T")
    q = sp.sympify(orbit["field_polynomial"], locals={"T":T})
    axis = sp.Matrix([sp.sympify(t, locals={"T":T}) for t in orbit["axis"]])
    cross = sp.Matrix(ell).cross(axis)
    test = sum(c*separator**i for i,c in enumerate(cross))
    return sp.resultant(q,test,T)


def chart(F, pivot, orbits=()):
    """Disjoint chart whose last nonzero differential coordinate has index pivot."""
    kept = sp.symbols(f"e0:{pivot}")
    ell = (*kept, sp.S.One, *([sp.S.Zero]*(2-pivot)))
    nuisance = []
    def fresh(name,count):
        values = sp.symbols(f"{name}0:{count}")
        nuisance.extend(values)
        return values
    basis = [m for m in powers(4) if m != max(F)]
    forms = [dict(zip(basis,fresh(prefix,len(basis)))) for prefix in ("aa","bb","cc")]
    a,b,lam,xt,yt = fresh("t",5)
    units = [4*a**3+27*b**2, lam, yt]
    if orbits:
        sep, = fresh("s",1)
        units.extend(axis_unit(orbit,ell,sep) for orbit in orbits)
        if any(t == 0 for t in units):
            return nuisance, kept, [sp.S.One]
    eqs = identities(F,forms,a,b,ell,lam)
    H = fresh("h",3)
    L = fresh("l",3)
    points = [fresh(f"p{i}c",3) for i in range(8)]
    values = []
    for i,point in enumerate(points):
        dot = lambda line: sum(s*t for s,t in zip(line,point))
        eqs.extend((evaluate(F,point),dot(H)-1))
        values.append(dot(L))
        A,B,C = (evaluate(f,point) for f in forms)
        if i < 4:
            eqs.extend((A,B,C))
        else:
            eqs.extend((A-xt*C,B-yt*C))
            units.append(C)
    for group in (values[:4],values[4:]):
        units.extend(x-y for x,y in combinations(group,2))
    eqs.extend(inv*f-1 for inv,f in zip(fresh("inv",len(units)),units))
    return nuisance, kept, [e for e in eqs if e != 0]


def write_singular(path, nuisance, kept, eqs):
    names = ",".join(map(str,(*nuisance,*kept)))
    order = f"(dp({len(nuisance)}),lp({len(kept)}))" if kept else "dp"
    with path.open("w") as out:
        out.write(f"ring R=0,({names}),{order};\nideal I=\n")
        out.write(",\n".join(str(e).replace("**","^") for e in eqs))
        out.write(";\n")
        if kept:
            out.write("ideal J=eliminate(I,"+"*".join(map(str,nuisance))+");\nJ;\n")
        else:
            out.write('ideal G=std(I);\nif (reduce(1,G)==0) {print("empty");} else {print("present");}\n')
        out.write("quit;\n")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quartic",type=Path)
    parser.add_argument("pivot",type=int,choices=range(3))
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--singular",type=Path,help="write an exact elimination job")
    group.add_argument("--solve",action="store_true",help="run an unrestricted lexicographic elimination with SymPy")
    args=parser.parse_args()
    data = json.loads(args.quartic.read_text())
    F = source(data)
    sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "involution-atlas"))
    from classify import classify
    nuisance, kept, eqs = chart(F,args.pivot,classify(data)["orbits"])
    if args.singular:
        write_singular(args.singular,nuisance,kept,eqs)
        print(f"{len(eqs)} equations, {len(nuisance)} auxiliary variables, {len(kept)} line coordinates")
    else:
        G = sp.groebner(eqs,*nuisance,*kept,domain=sp.QQ,order="lex")
        equations = [p for p in G if p.free_symbols <= set(kept)]
        print(json.dumps({"pivot":args.pivot,"elimination_generators":list(map(str,equations))},indent=2))


if __name__ == "__main__":
    main()
