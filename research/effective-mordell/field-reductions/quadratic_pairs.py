"""Enumerate all exceptional quadratic abscissas of C6 over Q and Q(i).

The odd chart is z=1/(x-1), w=z^3 y. Each output pair (u,v) describes the
points u(z)=0, w=v(z). Gaussian coefficients are [real, imaginary] as rational
strings. The two points of each irreducible u are conjugate over the base.
"""
import json
from collections import Counter,deque
from pathlib import Path
from cypari2 import Pari

p=Pari()
z,b=p('x'),p('w')
i=p.Mod(b,b*b+1)
f=((z+1)**6-z**6)*i**0
nf=p.nfinit('w^2+1')

def check(c,msg):
    if not c: raise ArithmeticError(msg)

def monic(g):return g/p.pollead(g,z)
def degree(g):return int(p.poldegree(g,z)) if g!=0 else -1

def gaussian(c):
    c=p.lift(c)
    return [str(p.polcoef(c,j,b)) for j in (0,1)]

def coeffs(g):
    return [gaussian(p.polcoef(g,j,z)) for j in range(degree(g)+1)]

def add(D,H):
    u,v=D
    s,t=H
    a,b,d=p.gcdext(u,s)
    c,h,e=p.gcdext(d,v+t)
    new_u=monic(u*s/e**2)
    new_v=((c*(a*u*t+b*s*v)+h*(v*t+f))/e)%new_u
    while degree(new_u)>2:
        new_u=monic((f-new_v**2)/new_u)
        new_v=(-new_v)%new_u
    check((f-new_v**2)%new_u==0,'divisor off curve')
    return new_u,new_v

def group(gaussian_base):
    H=[(z+p('1/2'),0),(z,1),(3*z*z+3*z+1,0)]
    if gaussian_base:H.append((z+1,i))
    G=[(monic(u*i**0),v*i**0) for u,v in H]
    identity=(z**0*i**0,0*i**0)
    def key(D):return json.dumps([coeffs(g) for g in D])
    seen={key(identity):identity};q=deque([identity])
    bound=36 if gaussian_base else 12
    while q:
        D=q.popleft()
        for T in G:
            v=add(D,T);label=key(v)
            if label not in seen:
                seen[label]=v;q.append(v)
                check(len(seen)<=bound,'torsion bound exceeded')
    check(len(seen)==bound,'group incomplete')
    return list(seen.values())

def split(disc):
    return len(p.nfroots(nf,z**2-p.lift(disc)))>0

def verify():
    for gaussian_base,primes,upper in [(False,(5,7),12),(True,(5,13),36)]:
        from math import gcd
        sizes=[int(p.ellinit([0,-1],q).ellcard()*p.ellinit([0,1],q).ellcard()) for q in primes]
        check(gcd(*sizes)==upper,'bad reduction bound')
        ds=group(gaussian_base)
        rows=[]
        for u,v in ds:
            if degree(u)!=2:continue
            disc=p.poldisc(u,z)
            if gaussian_base:
                if split(disc):continue
                fields=[n for n in (2,3,7) if split(disc/n)]
            else:
                disc=p.polcoef(p.lift(disc),0,b)
                if p.issquare(disc):continue
                fields=[n for n in (-2,-3) if p.issquare(disc/n)]
            check(len(fields)==1,'unexpected exceptional quadratic')
            rows.append(dict(extension=fields[0],u=coeffs(u),v=coeffs(v)))
        expected={2:4,3:2,7:12} if gaussian_base else {-2:2,-3:2}
        check(Counter(row['extension'] for row in rows)==expected,'quadratic orbit count')
        print('Base',('Q(i)' if gaussian_base else 'Q'),'group',len(ds),'exceptional divisor pairs',expected)
        if gaussian_base:
            path=Path(__file__).with_name('quadratic-pairs-over-gaussian.json')
            rows.sort(key=lambda r:(r['extension'],r['u'],r['v']))
            if path.exists():check(json.loads(path.read_text())==rows,'archived quadratic divisors')
            else:path.write_text(json.dumps(rows,indent=2)+'\n')

if __name__=='__main__':verify()
