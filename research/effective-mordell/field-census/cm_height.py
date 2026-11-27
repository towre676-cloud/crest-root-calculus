"""Certify an abscissa height window over K and optionally exhaust its polynomial box."""
import argparse
from math import comb, gcd
import json
from cypari2 import Pari
from rank_one_height import verify_identities
from classgroup_descent import classify as selmer

p=Pari()
w,x=p('w'),p('x')

def rank_bound(f):
    degree=int(f.poldegree())
    chars={int(p.core(p.poldisc(row[0]))) for row in p.nfsubfields(f,2)} if degree>1 else set()
    if len(chars)+1==degree:
        ranks={str(d):list(map(int,p.ellinit([0,-d**3]).ellrank()[:2])) for d in sorted({1,*chars})}
        return sum(interval[1] for interval in ranks.values()),dict(rational_character_ranks=ranks)
    cert=selmer(list(map(int,f.Vecrev())))
    return cert['rank_upper_bounds']['-1'],dict(selmer_sizes=cert['selmer_sizes']['-1'])

def certificate(co):
    if not isinstance(co,list) or len(co)<2 or co[-1]!=1 or any(type(c)!=int for c in co):
        raise ValueError('Supply ascending integral coefficients of a monic field polynomial.')
    f=p.Polrev(co,w)
    if not p.polisirreducible(f):raise ValueError('Field polynomial must be irreducible.')
    evidence=[]
    for cyclo,r_max,box in [(4,1,5),(12,2,5)]:
        # Cyclotomic fields are normal; all choices give Q-isomorphic composita.
        group=p.polcompositum(f,p.polcyclo(cyclo,w))
        L=group[0]
        rank,proof=rank_bound(L)
        evidence.append(dict(cyclotomic_order=cyclo,compositum=str(L),rank_upper_bound=rank,**proof))
        if rank<=r_max:
            return dict(field=str(f),status='height certified',strict_multiplicative_height_bound=box,
                        certificate=evidence)
    return dict(field=str(f),status='height criterion unresolved',certificate=evidence)


def exhaustive_points(co,result):
    # Absolute Mahler measure: each coefficient of a primitive minimal
    # polynomial of degree e is bounded by binomial(e,j) H(alpha)^e.
    bound=result['strict_multiplicative_height_bound']
    f=p.Polrev(co,w)
    nf=p.nfinit(f)
    n=len(co)-1
    points=set()
    for e in range(1,n+1):
        if n%e:continue
        bounds=[comb(e,j)*bound**e for j in range(e+1)]
        def visit(j,coeff):
            if j<0:
                if gcd(*coeff)!=1:return
                g=p.Polrev(coeff,x)
                if not p.polisirreducible(g):return
                for a in p.nfroots(nf,g):
                    for b in p.nfroots(nf,x*x-a**6+1):
                        if b*b!=a**6-1:raise ArithmeticError('Invalid point')
                        points.add((str(p.lift(a)),str(p.lift(b))))
                return
            lo=1 if j==e else -bounds[j]
            for c in range(lo,bounds[j]+1):
                coeff[j]=c
                visit(j-1,coeff)
        visit(e,[0]*(e+1))
    return dict(result,status='complete',verified_affine_points=sorted(points),rational_infinities=2)


if __name__=='__main__':
    app=argparse.ArgumentParser(description=__doc__)
    app.add_argument('ascending_polynomial')
    app.add_argument('--enumerate',action='store_true',help='Exhaust the finite absolute coefficient boxes. The work grows rapidly with degree.')
    args=app.parse_args()
    verify_identities()
    co=json.loads(args.ascending_polynomial)
    result=certificate(co)
    if args.enumerate and result['status']=='height certified':
        result=exhaustive_points(co,result)
    print(json.dumps(result,indent=2))
