"""Classify C6(K) for supported fields, with K given by an irreducible monic equation."""
import argparse
import json
from fractions import Fraction
from math import gcd, isqrt
from functools import lru_cache
from pathlib import Path
from cypari2 import Pari
from quadratic_descent import classify_twist

p=Pari()
w=p('w')
xvar=p('x')

def sqrt_rational(value):
    r=Fraction(str(value))
    if r<0:return None
    a,b=isqrt(r.numerator),isqrt(r.denominator)
    return p(a)/b if a*a==r.numerator and b*b==r.denominator else None

def signed_squarefree(n):
    n=int(n);f=p.factorint(abs(n));v=-1 if n<0 else 1
    for q,e in zip(f[0],f[1]):
        if not p.isprime(q):raise ArithmeticError('unproved factor in field radicand')
        if int(e)%2:v*=int(q)
    return v

def allowed(d):
    return all((int(q)==3 or int(q)%12==5) and p.isprime(q) for q in p.factorint(abs(d))[0])

def torsion_points(E):
    size,orders,gens=E.elltors();pts=[p([0])]
    for n,T in zip(orders,gens):
        pts=[E.elladd(P,E.ellmul(T,j)) for P in pts for j in range(int(n))]
    if len(set(str(P) for P in pts))!=int(size):raise ArithmeticError('torsion enumeration')
    return pts

@lru_cache(maxsize=1)
def bounded_quadratic_abscissas():
    """Small-height witnesses for fields where the descent leaves an open component."""
    by_field={}
    for a in range(1,49):
        for b in range(1,49):
            if gcd(a,b)!=1 or a==b:continue
            radicand=signed_squarefree(a**6-b**6)
            by_field.setdefault(radicand,[]).append(Fraction(a,b))
    return by_field

def classify(co, local_descent=False):
    if not isinstance(co,list) or len(co)<2 or any(type(c)!=int for c in co) or co[-1]!=1:
        raise ValueError('Use an ascending JSON list of integer coefficients with leading coefficient 1.')
    f=p.Polrev(co,w)
    if not p.polisirreducible(f):raise ValueError('The defining polynomial must be irreducible.')
    degree=len(co)-1
    subfields=[] if degree==1 else sorted({signed_squarefree(p.poldisc(row[0]))
                                      for row in p.nfsubfields(f,2)})
    ans=dict(field=str(f),degree=degree,quadratic_subfields=subfields)
    root_unit=p.Mod(1,f)
    pts={}
    def add(a,b):
        a,b=a*root_unit,b*root_unit
        if b*b!=a**6-1:raise ArithmeticError('point failed the source equation')
        key=(str(p.lift(a)),str(p.lift(b)))
        pts[key]=list(key)
    def root(d):
        rr=p.nfroots(p.nfinit(f),xvar**2-d)
        if not len(rr):raise ArithmeticError('element missing from declared subfield')
        return rr[0]
    for a in (1,-1):add(a,0)
    if -1 in subfields:
        i=root(-1);add(0,i);add(0,-i)
    if -3 in subfields:
        omega=root(-3)
        for h in (-1,1):
            for ordinate in (-omega,omega):add((h+ordinate)/2,0)
    complete=False
    if degree==len(subfields)+1 and all(allowed(d) for d in subfields):
        method='multiquadratic tower';complete=True
    elif degree==2:
        d=subfields[0]
        s=root(d)
        if d==-2:
            for a in [s/2,-s/2]:
                for y in [3*a/2,-3*a/2]:add(a,y)
        if d==-3:
            for h in (-1,1):
                for y in (-s,s):add((h+y)/2,0)
        if abs(d)==7:
            a=p(2) if d==7 else p(1)/2
            v=3*s if d==7 else 3*s/8
            for u in [-a,a]:
                for y in [-v,v]:add(u,y)
            complete=True;method='sqrt(7) descent and reciprocal isomorphism'
        else:
            method='quadratic divisors and rational elliptic quotient'
            for b in [-1,1]:
                E=p.ellinit([0,b*d**3])
                if int(E.ellrank()[1])!=0:continue
                for P in torsion_points(E):
                    if len(P)==1 or P[0]==0:continue
                    xx=P[0]/d if b==-1 else -p(d)/P[0]
                    a=sqrt_rational(xx)
                    if a is None:continue
                    v=sqrt_rational(d*(a**6-1))
                    if v is not None:
                        for X in [-a,a]:
                            for Y in [-v/s,v/s]:add(X,Y)
                complete=True;break
            if not complete:
                abscissas, open_pairs=classify_twist(d)
                for value in sorted(abscissas):
                    xx=p(str(value))
                    v=sqrt_rational(d*(xx**6-1))
                    if v is None:raise ArithmeticError('invalid twist lift')
                    add(xx,v/s)
                    add(xx,-v/s)
                complete=not open_pairs
                ans['remaining_squareclass_pairs']=open_pairs
                method='quadratic divisors and squareclass descent'
        if not complete:
            for value in bounded_quadratic_abscissas().get(d,[]):
                xx=p(str(value))
                v=sqrt_rational(d*(xx**6-1))
                if v is None:raise ArithmeticError('quadratic search lift')
                for X in [-xx,xx]:
                    for Y in [-v/s,v/s]:add(X,Y)
    elif degree==4 and subfields==[-7,-1,7]:
        method='exact rational and non-Gaussian abscissas; nonreal Gaussian abscissas unresolved'
        i,s=root(-1),root(7)
        for a,y in ((p(2),3*s),(p(1)/2,3*i*s/8)):
            for xx in (-a,a):
                for yy in (-y,y):add(xx,yy)
        nf=p.nfinit(f)
        rows=json.loads(Path(__file__).with_name('quadratic-pairs-over-gaussian.json').read_text())
        def gaussian(pair):return p(pair[0])+i*p(pair[1])
        for row in rows:
            if row['extension']!=7:continue
            U=p.Polrev([gaussian(c) for c in row['u']],xvar)
            V=p.Polrev([gaussian(c) for c in row['v']],xvar) if row['v'] else p(0)
            roots=p.nfroots(nf,U)
            if len(roots)!=2:raise ArithmeticError('quadratic divisor failed to split')
            for z in roots:add(1+1/z,p.subst(V,xvar,z)/z**3)
        if len(pts)!=36:raise ArithmeticError('biquadratic point count')
    else:
        method='outside the field classes proved by this program'
        if degree==len(subfields)+1:
            from isogeny_over_field import lift_torsion
            for sign in (-1,1):
                if any(int(p.ellinit([0,sign*d**3]).ellrank()[1]) for d in [1]+subfields):continue
                for X,Y in lift_torsion(p.nfinit(f),sign):add(p(X),p(Y))
                method='rank-zero rational character spaces and exact torsion'
                complete=True
                break
        if not complete and (degree<=5 or local_descent):
            from isogeny_over_field import classify as local_bound
            result=local_bound(co,max_residues=512)
            ans['local_descent']=result
            if result['status']=='complete':
                for X,Y in result['verified_affine_points']:add(p(X),p(Y))
                complete=True
                method='certified number-field 2-isogeny descent'
    ans.update(status='complete' if complete else 'partial',method=method,
               verified_affine_points=list(pts.values()),rational_infinities=2)
    return ans

if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('ascending_polynomial',help='For example "[-7,0,1]" for Q(sqrt(7)).')
    cli.add_argument('--local-descent',action='store_true',help='Also try class-number-one local descent above degree five.')
    a=cli.parse_args()
    try: print(json.dumps(classify(json.loads(a.ascending_polynomial),a.local_descent),indent=2))
    except ValueError as e:cli.error(str(e))
