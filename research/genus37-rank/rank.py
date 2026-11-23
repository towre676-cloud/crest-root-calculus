"""Prove Mordell-Weil rank lower bounds for the character factors of the genus-37 example."""
from fractions import Fraction
from itertools import combinations
from math import comb
from pathlib import Path
import json

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcd, gf_diff, gf_monic, gf_rem

from mumford import Jacobian, IDENTITY

HERE=Path(__file__).resolve().parent
CUBICS=[(15,-34,1),(10,-19,4),(11,3,-14),(11,-11,16)]
x=sp.Symbol("x")


def polynomial(S):
    return sp.Poly(sp.prod(x**3+CUBICS[i][0]*x**2+CUBICS[i][1]*x+CUBICS[i][2] for i in S),x)


def rational_poly(row,key):
    return sp.Poly.from_list([sp.Rational(c) for c in row[key]],x,domain=sp.QQ)


def validate(F,D):
    if D.get("infinity"):
        if F.degree()%2 or F.LC()!=1:
            raise ArithmeticError("Two distinct rational infinities required")
        return
    u,v=rational_poly(D,"u"),rational_poly(D,"v")
    if not 1 <= u.degree() <= (F.degree()-1)//2 or u.LC()!=1 or v.degree()>=u.degree():
        raise ArithmeticError("Invalid effective divisor")
    if not u.gcd(u.diff()).degree()==0 or not (v*v-F).rem(u).is_zero:
        raise ArithmeticError("Invalid global divisor relation")


def modq(q,p):
    t=Fraction(q)
    return t.numerator*pow(t.denominator,-1,p)%p


def chart(F,p,divs):
    coeff=[int(c)%p for c in F.all_coeffs()]
    if not sp.isprime(p) or p==2 or len(gf_gcd(coeff,gf_diff(coeff,p,ZZ),p,ZZ))>1:
        raise ValueError("Bad prime")
    g=(F.degree()-1)//2
    if F.degree()%2:
        J=Jacobian(coeff,p)
        pts=[]
        for D in divs:
            us=[modq(a,p) for a in D["u"]]
            vs=gf_rem([modq(b,p) for b in D["v"]],us,p,ZZ)
            pts.append((tuple(us),tuple(vs)))
    else:
        polys=[None if D.get("infinity") else ([modq(a,p) for a in D["u"]],
                                               [modq(b,p) for b in D["v"]]) for D in divs]
        def val(c,t):
            result=0
            for a in c:result=(result*t+a)%p
            return result
        root=next((r for r in range(p) if val(coeff,r)==0 and all(
            d is None or val(d[0],r)!=0 for d in polys)),None)
        if root is None:
            raise ValueError("No available rational branch point for this chart")
        def shifted(c):
            return [sum(a*comb(i,k)*pow(root,i-k,p) for i,a in enumerate(reversed(c)) if i>=k)%p
                    for k in range(len(c))]
        shift=shifted(coeff)
        L=shift[1]
        # X=L/(x-r), Y=y*X^(g+1)/L.
        Fnew=[v*pow(L,k-2,p)%p for k,v in enumerate(shift) if k>0]
        J=Jacobian(Fnew,p)
        pts=[]
        for d in polys:
            if d is None:
                pts.append(((1,0),(pow(L,g,p),)))
                continue
            u,v=d
            unew=gf_monic([c*pow(L,k,p)%p for k,c in enumerate(shifted(u))],p,ZZ)[1]
            vshift=shifted(v)
            vnew=[c*pow(L,k-1,p)%p for k,c in enumerate(vshift)]
            vnew+=[0]*(g+2-len(vnew))
            pts.append((tuple(unew),tuple(gf_rem(vnew,unew,p,ZZ))))
    for D in pts:J.check(D)
    # Twice this local class is globally defined as D-iota(D).
    return J,[J.mul(2,D) for D in pts]


def row_rank(rows,p):
    rows=[list(r) for r in rows];rank=0
    for col in range(len(rows[0]) if rows else 0):
        j=next((j for j in range(rank,len(rows)) if rows[j][col]%p),None)
        if j is None:continue
        rows[rank],rows[j]=rows[j],rows[rank]
        scale=pow(rows[rank][col],-1,p)
        rows[rank]=[v*scale%p for v in rows[rank]]
        for j in range(rank+1,len(rows)):
            m=rows[j][col]
            rows[j]=[(a-m*b)%p for a,b in zip(rows[j],rows[rank])]
        rank+=1
    return rank


def verify_sector(S,data):
    F=polynomial(S)
    divs=data["divisors"]
    for D in divs:validate(F,D)
    controls=data["controls"]
    if len({p for p,m in controls})!=2:
        raise ArithmeticError("Two distinct control primes are required")
    B=1
    for p,m in controls:
        if m<=0:raise ArithmeticError("Invalid annihilator")
        J,Ds=chart(F,p,divs)
        if any(J.mul(m,D)!=IDENTITY for D in Ds):
            raise ArithmeticError("Control multiplier failed")
        B*=m
    ell=data["ell"]
    if not sp.isprime(ell):
        raise ArithmeticError("Expected a prime coefficient field")
    rows=[]
    for p,m,expected in data["witnesses"]:
        if m<=0 or p==ell or B%ell==0:raise ArithmeticError("Invalid witness prime or multiplier")
        J,Ds=chart(F,p,divs)
        images=[J.mul(m,D) for D in Ds]
        U=next((D for D in images if D!=IDENTITY),None)
        if U is None or J.mul(ell,U)!=IDENTITY:
            raise ArithmeticError("Missing cyclic witness")
        multiples={J.mul(i,U):i for i in range(ell)}
        row=[multiples[D] for D in images]
        if row!=expected:raise ArithmeticError("Incorrect reduction relation")
        rows.append(row)
    rank=row_rank(rows,ell)
    if rank!=len(divs) or rank!=data["rank_lower_bound"]:
        raise ArithmeticError("Reduction rows fail to prove the stated rank")
    return rank


def main():
    d=json.loads((HERE/"certificate.json").read_text())
    ranks={key:verify_sector(tuple(map(int,key)),v) for key,v in d["sectors"].items()}
    print(json.dumps({"bounds":ranks,"total_rank_lower_bound":sum(ranks.values())},indent=2))


if __name__=="__main__":
    main()
