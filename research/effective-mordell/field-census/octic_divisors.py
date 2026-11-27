"""Classify all exceptional quadratic abscissas over Q(zeta12) and rational squareclasses."""
from collections import deque,Counter
from fractions import Fraction
import json
from pathlib import Path

from cypari2 import Pari

p=Pari()
z,w=p('x'),p('w')
base=w**4-w*w+1
nf=p.nfinit(base)
unit=p.Mod(1,base)
zz=p.Mod(w,base)
f=((z+1)**6-z**6)*unit

def degree(g):
    return -1 if g==0 else int(p.poldegree(g,z))

def monic(g):return g/p.pollead(g,z)

def add(D,H):
    u,v=D
    s,t=H
    a,b,d=p.gcdext(u,s)
    c,h,e=p.gcdext(d,v+t)
    U=monic(u*s/e**2)
    V=((c*(a*u*t+b*s*v)+h*(v*t+f))/e)%U
    while degree(U)>2:
        U=monic((f-V**2)/U)
        V=(-V)%U
    if (f-V**2)%U:raise ArithmeticError("Divisor is off the odd source chart")
    return U,V

def group(G):
    identity=(z**0*unit,0*unit)
    def key(D):
        return tuple(tuple(str(p.lift(p.polcoef(g,j,z))) for j in range(degree(g)+1)) for g in D)
    seen={key(identity):identity}
    queue=deque([identity])
    while queue:
        D=queue.popleft()
        for T in G:
            P=add(D,T)
            k=key(P)
            if k not in seen:
                seen[k]=P
                queue.append(P)
    return list(seen.values())

def rational_discriminants(dc):
    """Outside ramified rational primes 2,3, valuations determine d uniquely."""
    factors=p.idealfactor(nf,dc)
    primes=set()
    for pr,e in zip(factors[0],factors[1]):
        q=int(pr[0])
        if q>3 and int(e)%2:primes.add(q)
    d=1
    for q in sorted(primes):d*=q
    return [n for n in (d,2*d) if p.nfeltissquare(nf,dc/n)]

def coefficients(g):
    return [str(p.lift(p.polcoef(g,j,z))) for j in range(degree(g)+1)]

def certify():
    G=[(z,unit),(z+1,zz**3)]
    for j in range(1,6):
        X=zz**(2*j)
        G.append((z-1/(X-1),0*unit))
    H=group(G)
    if len(H)!=144:raise ArithmeticError("Cyclotomic Jacobian group order")
    # Prime-to-q torsion injects at good reduction, both 13 and 37 split completely.
    from math import gcd
    counts=[int(p.ellinit([0,-1],q).ellcard())**2 for q in (13,37)]
    if gcd(*counts)!=144:raise ArithmeticError("Jacobian torsion upper bound")
    records=[]
    for u,v in H:
        if degree(u)!=2:continue
        dc=p.poldisc(u,z)
        if dc==0 or p.nfeltissquare(nf,dc):continue
        rads=rational_discriminants(dc)
        if len(rads)>1:raise ArithmeticError("Nonindependent radicands")
        records.append(dict(u=coefficients(u),v=coefficients(v),rational_class=rads[0] if rads else None))
    if len(records)!=96:raise ArithmeticError("Quadratic divisor count")
    return records
if __name__=='__main__':
    result=certify()
    print(dict(Counter(r['rational_class'] for r in result)))
    Path('octic-quadratic-divisors.json').write_text(json.dumps(result,indent=2)+'\n')
