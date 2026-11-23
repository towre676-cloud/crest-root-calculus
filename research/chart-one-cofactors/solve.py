"""Compute a finite triangular certificate for a marked transverse chart."""
import argparse
import json
from pathlib import Path

import sympy as sp
from sympy.polys.agca.extensions import MonogenicFiniteExtension
from sympy.polys.rings import ring

from chart import build, C, S, R, c, s


def coefficients(p):
    return [str(p.nth(i)) for i in range(p.degree()+1)] if p else []


def determinant3(m):
    a,b,c = m[0]
    d,e,f = m[1]
    g,h,i = m[2]
    return a*(e*i-f*h)-b*(d*i-f*g)+c*(d*h-e*g)


def norm_with_cofactor(q,row):
    """Compute e,u with row*u=e mod q via four 3-by-3 determinants."""
    K = q.ring.domain
    U,x = ring("s",K)
    columns = []
    v = row
    for _ in range(4):
        col = [U.zero for _ in range(4)]
        for (i,j),a in v.items():
            col[i][j,] = a
        columns.append(col)
        v = (v*q.ring.gens[0]).rem([q])
    adj = [(-1)**j*determinant3([col[1:] for i,col in enumerate(columns) if i!=j])
           for j in range(4)]
    norm = sum((col[0]*a for col,a in zip(columns,adj)),U.zero)
    cofactor = q.ring.from_dict({(i,j):a for i,u in enumerate(adj) for (j,),a in u.items()})
    observed = (row*cofactor).rem([q])
    if observed != q.ring.from_dict({(0,j):a for (j,),a in norm.items()}):
        raise ArithmeticError("The determinant cofactor does not reproduce the norm.")
    return sp.Poly.from_dict(norm.to_dict(),(S,),domain=K),cofactor


def select_norms(chart):
    """Some combination of the four residuals has nonzero norm on the generic line."""
    slots = []
    for j in (1,3,0,2):
        if chart.rows[j]:
            recipe = [int(i==j) for i in range(4)]
            slots.append((chart.rows[j],recipe))
    # A generic scalar polynomial has at most 3 exceptional scalars per component of q.
    for k in range(13):
        recipe = [k**i for i in range(4)]
        slots.append((sum((a*r for a,r in zip(recipe,chart.rows)),chart.q.ring.zero),recipe))
    used = []
    for row,recipe in slots:
        if not row or any(row==r[0] for r in used):
            continue
        n,u = norm_with_cofactor(chart.q,row)
        if n:
            used.append((row,recipe,n,u))
            if len(used)==2 or n.degree()==0:
                break
    if not used:
        raise ArithmeticError("Every tested generic norm vanished.")
    return used


class ResidueField(MonogenicFiniteExtension):
    """QQ[s]/(f) for a verified irreducible f, using field exact division."""
    def exquo(self,a,b):
        return a/b


def specialize(p, factor, field):
    """Reduce the coefficients in s before entering the number field."""
    parts = {}
    for (i,j),v in p.items():
        parts.setdefault(i,{})[j,] = v
    return sp.Poly.from_dict({(i,):sp.Poly.from_dict(a,(S,),domain=sp.QQ).rem(factor).as_expr()
                              for i,a in parts.items()},(C,),domain=field)


def lifted(p):
    """Canonical c-polynomial lift, with s-degree below the number field degree."""
    vals = {}
    for (i,),co in p.rep.to_dict().items():
        for (j,),v in co.rep.to_dict().items():
            vals[i,j] = v
    return R.from_dict(vals)


def sparse(p):
    return [[i,j,str(a)] for (i,j),a in sorted(p.items())]


def profiled_gcd(polynomials):
    """Return d and the explicit Bezout row d = sum u_j f_j over the number field."""
    unit = polynomials[0].one
    zero = polynomials[0].zero
    d = polynomials[0]
    row = [unit]
    for p in polynomials[1:]:
        if d.degree()==0 or not p:
            row.append(zero)
            continue
        # Divide before gcdex to limit intermediate coefficient growth.
        quo, rem = p.div(d)
        if not rem:
            row.append(zero)
            continue
        a,b,new = d.gcdex(rem)
        row = [(a-b*quo)*u for u in row] + [b]
        d = new
    denom = d.LC()
    if denom!=1:
        k = sp.Poly(denom,C,domain=d.domain)
        row = [u.exquo(k) for u in row]
        d = d.monic()
    if sum((u*p for u,p in zip(row,polynomials)),zero)!=d:
        raise ArithmeticError("Fibre Bezout identity failed.")
    return d,row


def good_prime_for_gcd(norms,g):
    """An invertible Sylvester matrix of the quotients certifies the exact gcd."""
    quot = [n.exquo(g).clear_denoms(convert=True)[1].primitive()[1] for n in norms]
    p = 991
    while True:
        p = int(sp.nextprime(p))
        polys = [sp.Poly(a.as_expr(),S,modulus=p) for a in quot]
        if any(a.degree()!=b.degree() for a,b in zip(quot,polys)):
            continue
        rem = polys[0]
        for a in polys[1:]:
            rem = rem.gcd(a)
        if rem.degree()==0:
            return int(p)


def quotient_targets(chart, factor, field, good, specialization=None):
    """Recover b and all coefficients of the original smooth cubic w²=K(z)."""
    if good.degree()==0:
        return None
    red = lambda a:a.rem(good)
    mul = lambda a,b:red(a*b)
    spec = specialization or (lambda p:specialize(p,factor,field))
    h = red(spec(chart.h))
    invh = h.invert(good)
    denom = red(spec(chart.denominator))
    b = mul(spec(chart.numerator),denom.invert(good))
    table = {k:red(spec(v)) for k,v in chart.table.items()}
    bp = [good.one]
    for _ in range(6):
        bp.append(mul(bp[-1],b))
    from math import comb
    def g(i,l):
        w = 6-i-2*l
        return sum((comb(l+w-j,l)*mul(table[j,l+w-j],bp[w-j])
                    for j in range(w+1) if j<=4 and l+w-j<=3),good.zero)
    half = sp.Rational(1,2)
    e2 = mul(invh,g(5,0))*half
    e1 = (mul(invh,g(4,0))-mul(e2,e2))*half
    f1 = mul(invh,g(4,1))*half
    e0 = mul(invh,g(3,0))*half-mul(e1,e2)
    f0 = mul(invh,g(3,1))*half-mul(f1,e2)
    square = {(3,0):good.one,(2,0):e2,(1,0):e1,(1,1):f1,(0,0):e0,(0,1):f0}
    K = []
    for i in range(7):
        for j in range(4):
            val = sum((mul(a,b) for (k,l),a in square.items() for (u,v),b in square.items()
                       if k+u==i and l+v==j),good.zero)-mul(invh,g(i,j))
            if i and val:
                raise ArithmeticError("The reconstructed square fails the full original equation.")
            if i==0:
                K.append(val)
    a,b1,c1,d = K
    disc = red(b1*b1*c1*c1-4*a*c1**3-4*d*b1**3-27*a*a*d*d+18*a*b1*c1*d)
    if disc.gcd(good).degree()!=0 or d.gcd(good).degree()!=0:
        raise ArithmeticError("A reconstructed target is singular.")
    out = lambda p:sparse(lifted(p))
    return {"b":out(b),"square":[out(v) for v in (e2,e1,f1,e0,f0)],
            "target_constant_first":[out(v) for v in K]}


def fibre(chart,f,multiplicity):
    field = ResidueField(f)
    q,h = [specialize(p,f,field) for p in [chart.q,chart.h]]
    rows = [specialize(p,f,field) for p in chart.rows]
    d,bezout = profiled_gcd([q,*rows])
    good = d.exquo(d.gcd(h**4)).monic()
    if good.gcd(h).degree()!=0 or good.gcd(good.diff()).degree()!=0:
        raise ArithmeticError("Boundary removal did not produce a reduced open fibre.")
    return {"s_factor":coefficients(f),"s_multiplicity_in_eliminant":multiplicity,
            "fibre_gcd":sparse(lifted(d)),"unit_fibre":sparse(lifted(good)),
            "geometric_marks":f.degree()*good.degree(),
            "parameters":quotient_targets(chart,f,field,good)}, [lifted(u) for u in bezout]


def describe(data,point,details=False):
    chart = build(data,point)
    selected = select_norms(chart)
    norms = [n for r,p,n,u in selected]
    g = norms[0].monic()
    for n in norms[1:]:
        g = g.gcd(n).monic()
    prime = good_prime_for_gcd(norms,g)
    fibres,bezout = [],[]
    for f,m in g.factor_list()[1]:
        record,row = fibre(chart,f.monic(),m)
        fibres.append(record)
        bezout.append(row)
    result = {"frame":chart.frame,"norm_rows":[p for _,p,_,_ in selected],
              "norm_degrees":[int(n.degree()) for n in norms],"gcd_degree":g.degree(),
              "coprimality_prime":prime,"fibres":fibres,
              "geometric_marks":sum(r["geometric_marks"] for r in fibres)}
    if details:
        result["norm_cofactors"] = [sparse(u) for _,_,_,u in selected]
        result["fibre_cofactors"] = [[sparse(u) for u in row] for row in bezout]
    return result


if __name__=="__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("quartic",type=Path)
    ap.add_argument("point",nargs=3,metavar="coordinate")
    ap.add_argument("--cofactors",action="store_true",help="include individual polynomial cofactor rows")
    ap.add_argument("--output",type=Path)
    args = ap.parse_args()
    result = describe(json.loads(args.quartic.read_text()),args.point,args.cofactors)
    payload = json.dumps(result,indent=2)+"\n"
    if args.output:
        args.output.write_text(payload)
    else:
        print(payload,end="")
