"""Height-certified rational PGL2 isomorphisms from split-prime branch triples."""
from itertools import permutations
from math import gcd, lcm

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_pow_mod

from frobenius import X
from partitions import sextic_chart
from richelot import primitive, transform


def integral_model(f):
    return f.clear_denoms(convert=True)[1].primitive()[1]


def root_bound(f):
    a=int(f.LC())
    # The polynomial a^5 f(x/a) is monic and integral. Its roots are a times the original roots.
    return 1+max(abs(int(f.nth(i))*a**(5-i)) for i in range(6))


def split_roots(f, g):
    p=5
    while True:
        p=int(sp.nextprime(p))
        if any(int(poly.LC())%p==0 for poly in (f,g)):
            continue
        if any(gf_pow_mod([1,0],p,[int(c)%p for c in poly.all_coeffs()],p,ZZ)!=[1,0] for poly in (f,g)):
            continue
        points=[[j for j in range(p) if int(poly.eval(j))%p==0] for poly in (f,g)]
        if all(len(s)==6 for s in points):
            return p,*points


def lift(f,roots,p,mod):
    d=f.diff()
    m=p
    while m<mod:
        m*=m
        roots=[(r-int(f.eval(r))*pow(int(d.eval(r)),-1,m))%m for r in roots]
    return roots


def determinant(a,b,c):
    return (a[0]*(b[1]*c[2]-b[2]*c[1])-a[1]*(b[0]*c[2]-b[2]*c[0])
            +a[2]*(b[0]*c[1]-b[1]*c[0]))


def cofactors(w,z,m):
    rows=[[s,1,-s*t,-t] for s,t in zip(w,z)]
    return [(-1)**i*determinant(*(row[:i]+row[i+1:] for row in rows))%m for i in range(4)]


def reconstruct(r,m,bound):
    if r==0:
        return sp.Integer(0)
    u,v,a,b=m,r,0,1
    while abs(v)>bound:
        k=u//v
        u,v,a,b=v,u-k*v,b,a-k*b
    if not b or abs(b)>bound or gcd(v,b)!=1 or (v-r*b)%m:
        return None
    return sp.Rational(v,b)


def align(f, g):
    f,cf=sextic_chart(f);g,cg=sextic_chart(g)
    f,g=integral_model(f),integral_model(g)
    A,B=abs(int(f.LC())),abs(int(g.LC()))
    R,S=root_bound(f),root_bound(g)
    H=6*max(A*S,A*B,R*S,B*R)**3
    p,z,w=split_roots(f,g)
    m=p
    while m<=2*H**2:
        m*=m
    z,w=lift(f,z,p,m),lift(g,w,p,m)
    out=[]
    for triple in permutations(z,3):
        v=cofactors(w[:3],triple,m)
        pivot=next(i for i,u in enumerate(v) if u%p)
        a,b,c,d=[u%p for u in v]
        if any((a*r+b)*pow((c*r+d)%p,-1,p)%p not in {s%p for s in z}
               if (c*r+d)%p else True for r in w):
            continue
        ratios=[reconstruct(u*pow(v[pivot],-1,m)%m,m,H) for u in v]
        if None in ratios:
            continue
        matrix=primitive(ratios)
        if max(map(abs,matrix))>H:
            continue
        h=sp.Poly(transform(f,matrix,6),X,domain=sp.QQ)
        if h.monic()==g.monic():
            original=sp.Matrix(2,2,cf)*sp.Matrix(2,2,matrix)*sp.Matrix(2,2,cg).inv()
            out.append(primitive(list(original)))
    return sorted(set(out))
