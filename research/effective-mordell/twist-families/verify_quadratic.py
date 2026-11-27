"""Certify the quartic cover, torsion fibres, and simple reduction over F_23.
Requires cypari2, SymPy and python-flint. See arithmetic.m for the K-rank.
"""
import json
from pathlib import Path
from itertools import product

import flint
import sympy as sp
from cypari2 import Pari

pari = Pari()
x,t,w = [pari(name) for name in ("x","t","w")]
K = pari.Mod(1,w**2+2)


def check(value, label):
    if not value:
        raise ArithmeticError(label)


def point_count(co,p,degree):
    field=flint.fq_default_ctx(p,degree)
    poly=flint.fq_default_poly_ctx(field)([field(i) for i in co])
    total=1  # unique point at infinity on the odd-degree model
    for vals in product(range(p),repeat=degree):
        r=poly(field(list(vals)))
        total+=1 if r.is_zero() else 2*int(r.is_square())
    return total


def rational_gcd(g):
    # x is rational, so both coefficients on the Q-basis (1,w) must vanish.
    v=pari.lift(g)
    parts=[sum(pari.polcoef(pari.polcoef(v,i,x),j,w)*x**i
               for i in range(int(pari.poldegree(g,x))+1)) for j in range(2)]
    result=pari.gcd(*parts)
    return result/pari.pollead(result,x)


def verify():
    rec=json.loads(Path(__file__).with_name('quadratic-cover.json').read_text())
    co=rec['coefficients']
    f=sum(c*x**k for k,c in enumerate(co))
    center=rec['reciprocal_chart_center']
    fc=sum(c*(center*x+1)**i*x**(6-i) for i,c in enumerate(co))
    item=rec['maps'][0]
    n,s,h,g,omega=[pari(item[key])*K for key in
                   ('numerator','half_poles','lift','cubic','differential')]
    A,B=pari.polcoef(g,1,t),pari.polcoef(g,0,t)
    check(g==t**3+A*t+B and 4*A**3+27*B**2!=0,'target smoothness')
    check(n**3+A*n*s**4+B*s**6==fc*h**2,'cover equation')
    check(pari.poldegree(pari.gcd(n,s),x)==0 and
          max(pari.poldegree(n,x),2*pari.poldegree(s,x))==4,'degree four')
    check(pari.deriv(n,x)*s-2*n*pari.deriv(s,x)==h*omega and
          0<=pari.poldegree(omega,x)<=1,'holomorphic differential')
    check(1728*4*A**3/(4*A**3+27*B**2)==pari(item['j'])*K,'j-invariant')
    n0=(x-center)**4*pari.subst(n,x,1/(x-center))
    s0=(x-center)**2*pari.subst(s,x,1/(x-center))
    check(pari.poldegree(pari.gcd(n0,s0),x)==0,'original-coordinate base locus')
    ww=w*K
    torsion_x=[18*ww-33,-36*ww+30,18*ww+3]
    check(g==(t-torsion_x[0])*(t-torsion_x[1])*(t-torsion_x[2]),
          'three nonzero 2-torsion points')
    fibres=[rational_gcd(s0)]+[rational_gcd(n0-u*s0**2) for u in torsion_x]
    print('Affine rational-fibre gcds (pole, then the three listed torsion points):',fibres)
    all_x=set()
    for pol in fibres:
        sy=sp.sympify(str(pol).replace('^','**'))
        all_x.update(sp.ground_roots(sy,sp.Symbol('x')))
    check(all_x=={-1,0,8},'finite rational abscissas')
    check(all(pari.subst(f,x,int(u))==0 for u in all_x),'finite source points')
    p=23
    mod=sp.Poly.from_list(co[::-1],sp.Symbol('x'),modulus=p)
    check(mod.degree()==5 and mod.is_sqf,'good reduction')
    N1,N2=[point_count(co,p,k) for k in (1,2)]
    a1=p+1-N1
    a2=(a1*a1-(p*p+1-N2))//2
    poly=sp.Poly.from_list([1,-a1,a2,-p*a1,p*p],sp.Symbol('T'),domain=sp.QQ)
    check(poly.all_coeffs()==[1,0,14,0,529] and poly.is_irreducible,'Q-simple certificate')
    print('Point counts in F_23 and F_(23^2):',N1,N2,'; Frobenius polynomial:',poly.as_expr())
    print('All four rational points verified using the K-rank and torsion in arithmetic.m.')


if __name__=='__main__':
    verify()
