"""Exact points on C6 over Q(i,sqrt(d)) for squarefree d with inert odd support."""
from functools import lru_cache
from itertools import product
from math import prod
from pathlib import Path
import argparse
import json

from cypari2 import Pari

pari = Pari()
x = pari('x')
w = pari('w')
gaussian_field = pari.nfinit(w*w+1)
i = pari.Mod(w, w*w+1)


def support(d):
    if d <= 1 or int(pari.core(d)) != d:
        raise ValueError('Use a positive squarefree radicand greater than one.')
    primes = tuple(map(int, pari.factorint(d)[0]))
    if any(q != 2 and q%4 != 3 for q in primes):
        raise ValueError('The Gaussian squareclass table here requires all odd factors to be inert.')
    return primes


def divisors(primes):
    values = [1]
    for q in primes:
        values += [q*v for v in values]
    return sorted(values)


def gaussian_class(n):
    return (0,n//2) if n%2 == 0 else (n,0)


@lru_cache(None)
def residue_data(q,modulus):
    def mul(u,v):
        return ((u[0]*v[0]-u[1]*v[1])%modulus, (u[0]*v[1]+u[1]*v[0])%modulus)
    def add(u,v):
        return ((u[0]+v[0])%modulus, (u[1]+v[1])%modulus)
    def sub(u,v):
        return ((u[0]-v[0])%modulus, (u[1]-v[1])%modulus)
    def is_nonunit(u):
        if q == 2:
            return (u[0]+u[1])%2 == 0
        return u[0]%q == u[1]%q == 0
    elements = list(product(range(modulus), repeat=2))
    one = (1,0)
    proj = [(a,one) for a in elements]+[(one,b) for b in elements if is_nonunit(b)]
    forms = []
    for u,v in proj:
        uu,vv,uv = mul(u,u),mul(v,v),mul(u,v)
        forms.append((add(add(uu,vv),uv),sub(add(uu,vv),uv),sub(uu,vv)))
    return mul, {mul(z,z) for z in elements}, forms


def solve_local(d, classes):
    viable = set(product(classes,repeat=2))
    for q,modulus in ((2,16),(3,9),(7,49),(11,11),(19,19)):
        mul,squares,forms = residue_data(q,modulus)
        reps = {a:gaussian_class(a) for a in classes}
        tested = {}
        def ok(scale,value):
            key=(scale,value)
            if key not in tested:
                tested[key] = mul(gaussian_class(scale),value) in squares
            return tested[key]
        survivors=set()
        for a,b in viable:
            for A,B,D in forms:
                if ok(a,A) and ok(b,B):
                    # Do not remove square factors from d*a*b at residue primes.
                    # d, a and b individually choose the integral representative of the squareclass.
                    if mul(mul(gaussian_class(d),mul(reps[a],reps[b])),D) in squares:
                        survivors.add((a,b))
                        break
        viable=survivors
    return sorted(viable)


@lru_cache(None)
def elliptic_projection(kind,h):
    coeffs = [0,h,0,h*h,0] if kind=='even' else [0,9*h,0,30*h*h,36*h**3]
    E = pari.ellinit(coeffs)
    if any(int(T.ellrank()[1])>0 for T in (E,E.elltwist(-4))):
        return None
    EK=pari.ellinit(coeffs,gaussian_field)
    order,orders,gens=EK.elltors()
    torsion=[pari([0])]
    for n,T in zip(orders,gens):
        torsion=[EK.elladd(P,EK.ellmul(T,j)) for P in torsion for j in range(int(n))]
    if len(set(map(str,torsion)))!=int(order):
        raise ArithmeticError('Torsion enumeration failed.')
    candidates=set()
    for P in torsion:
        if len(P)==1:
            continue
        if kind=='even':
            candidates.update(str(pari.lift(z)) for z in pari.nfroots(gaussian_field,x*x-P[0]/h))
        elif P[0]!=0:
            candidates.add(str(pari.lift(1+6*h/P[0])))
    return sorted(candidates)


def solve(d):
    primes = support(d)
    classes=divisors(sorted({2,3,*primes}))
    pairs=solve_local(d,classes)
    exceptional=json.loads(Path(__file__).with_name('quadratic-pairs-over-gaussian.json').read_text())
    f=w**4-2*(d-1)*w*w+(d+1)**2
    nf=pari.nfinit(f)
    unit=pari.Mod(1,f)
    root=pari.Mod(w,f)
    # w = sqrt(d)+i. This choice fixes both inclusions in the output basis.
    imag=(root-(d+1)/root)/2
    rad=root-imag
    if imag**2!=-1 or rad**2!=d:
        raise ArithmeticError('Field embedding failed.')
    def include(z):
        return pari.subst(pari.lift(z),w,imag)*unit
    points=set()
    def add(X,Y):
        if Y**2!=X**6-1:
            raise ArithmeticError('Point off the source curve.')
        points.add((str(pari.lift(X*unit)),str(pari.lift(Y*unit))))
    for u in (1,-1):add(u,0)
    for y in (-imag,imag):add(0,y)
    for row in exceptional:
        if row['extension']!=d:continue
        def poly(coeffs):
            return sum((pari(a)+imag*pari(b))*x**j for j,(a,b) in enumerate(coeffs))
        U,V=poly(row['u']),poly(row['v'])
        roots=pari.nfroots(nf,U)
        if len(roots)!=2:raise ArithmeticError('Exceptional divisor did not split.')
        for z in roots:add(1+1/z,pari.subst(V,x,z)/z**3)
    open_pairs=[]
    for a,b in pairs:
        for kind,h,sign in [('even',a*b,1),('mixed',d*b,1),('mixed',d*a,-1)]:
            values=elliptic_projection(kind,int(pari.core(h)))
            if values is None:continue
            for value in values:
                X=sign*pari.Mod(pari(value),w*w+1)
                for v in pari.nfroots(gaussian_field,x*x-d*(X**6-1)):
                    add(include(X),include(v)/rad)
            break
        else:
            open_pairs.append((a,b))
    return dict(d=d,field=str(f),status='partial' if open_pairs else 'complete',
                local_pairs=pairs,remaining_pairs=open_pairs,affine_points=sorted(points),infinities=2)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('d',type=int)
    print(json.dumps(solve(cli.parse_args().d),indent=2))
