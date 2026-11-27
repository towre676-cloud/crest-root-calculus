"""Check identities, all 61 individual parameters, and a sample of the infinite theorem."""
import json
from itertools import combinations
from math import gcd, isqrt
from pathlib import Path

import sympy as sp
from descent import check, direct_rank, pari, source

HERE=Path(__file__).resolve().parent


def check_maps():
    x,d=sp.symbols('x d')
    f=source(x,d)
    maps=[(-108,-32*d*x/(x*x+1)**2,8*d*(3*x*x-1)/(x*x+1)**3,4),
          (-12,-32*d*x**3/(x*x+9)**2,-8*d*x*(x*x-27)/(x*x+9)**3,-4*x)]
    for a,X,H,omega in maps:
        check(sp.cancel(X**3+a*d*d*X-f*H**2)==0,'twist cover identity')
        check(sp.cancel(sp.diff(X,x)-omega*H)==0,'independent differentials')
    check(sp.cancel(x**6*source(3/x,3*d)-81*f)==0,'isomorphism to triple twist')


def full_panel():
    high={r['d']:r for r in json.loads((HERE/'panel.json').read_text())}
    branches=json.loads((HERE/'descent-receipt.json').read_text())
    rows=[]
    for d in range(1,101):
        if any(d%(k*k)==0 for k in range(2,isqrt(d)+1)):
            continue
        row=high[d] if d in high else direct_rank(d)
        check(row['torsion']==[2,2] and all(l==u for l,u in row['ranks']),
              'elliptic groups')
        affine=[[0,0]]
        if min(u for _,u in row['ranks'])==0:
            method='rational quartic quotient'
        else:
            branches_d=branches['squareclasses'][str(d)]
            check(all(not b.get('rational_x') for b in branches_d),'unexpected rank-zero lift')
            positive=sorted(b['a'] for b in branches_d if b.get('ranks',[0,0])[1])
            check(positive=={22:[1],55:[10],66:[3]}.get(d,[]),'unresolved squareclass')
            method='rational squareclass descent'
            if d in (22,55):
                survivors=[b['delta'] for b in branches['quadratic_sieve'][str(d)]
                           if b['obstruction'] is None]
                check(survivors==([[-5,-2]] if d==22 else []),'quadratic sieve')
                method='elliptic Chabauty' if d==22 else 'quadratic local obstructions'
            if d==66:
                method='isomorphism to d=22'
            if d in (22,66):
                x,y=(1,44) if d==22 else (3,396)
                check(source(x,d)==y*y,'exceptional points')
                affine.extend([[x,-y],[x,y]])
        rows.append(dict(**row,method=method,affine=affine,rational_infinities=1))
    check(len(rows)==61 and sum(len(r['affine'])==1 for r in rows)==59,'panel count')
    return rows


def uniform_examples():
    # The proof of the infinite result is in the accompanying proof. This is a regression.
    admissible=[p for p in sp.primerange(5,160) if p%12==5]
    special=[None]+[p for p in sp.primerange(7,110) if p%24==7]
    examples=0
    for k in range(3):
        for ps in combinations(admissible,k):
            u=int(sp.prod(ps))
            for q in special:
                n=u*(q or 1)
                if n>40000:
                    continue
                m=n if (u%8==(1 if q else 5)) else 3*n
                check(m%8 in ((1,7) if q else (3,5)),'choice of twist')
                check(gcd(11*19,m)==1,'torsion reduction primes')
                E=pari.ellinit([-12*m*m,0])
                check(int(E.ellrank()[1])==0 and int(E.elltors()[0])==2,
                      'uniform theorem disagrees with independent PARI rank computation')
                examples+=1
    return examples


if __name__=='__main__':
    check_maps()
    rows=full_panel()
    out=HERE/'all-points.json'
    if out.exists():
        check(rows==json.loads(out.read_text()),'final point receipt')
    else:
        out.write_text(json.dumps(rows,indent=2)+'\n')
    print('Complete individual parameters:',len(rows))
    print('Twists with rational non-Weierstrass points:',[r['d'] for r in rows if len(r['affine'])>1])
    print('Independent finite tests of the uniform theorem:',uniform_examples())
