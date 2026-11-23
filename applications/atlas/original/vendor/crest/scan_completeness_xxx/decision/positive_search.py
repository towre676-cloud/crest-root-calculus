"""Enumerate exact marked pencils over successively larger number fields.

The transverse 3-division polynomial supplies c and determines b at a given
point and fibre line. Only roots over the chosen number field are needed,
since all marked chart data are eventually enumerated in their own field.
"""
from functools import lru_cache
from itertools import product
from math import gcd
from pathlib import Path
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
OLD=next((p/'exhaustion_xxviii' for p in HERE.parents if (p/'exhaustion_xxviii').is_dir()),None)
if OLD is None:
    raise ImportError('Extract the integrated release, including exhaustion_xxviii.')
sys.path.insert(0,str(OLD))
import free_atlas as fa
import recover_maps as rm

T=rm.T
t,z=fa.t,fa.z
cvar=sp.Symbol('cvar')


def elements(d,height):
    for nums in product(range(-height,height+1),repeat=d):
        for denominator in range(1,height+1):
            if gcd(denominator,*nums)==1:
                yield sum(sp.Rational(a,denominator)*T**i for i,a in enumerate(nums) if a)


@lru_cache(maxsize=16)
def data_field(f):
    return rm.field_for(sp.Poly(f,T,domain=sp.QQ))


def scan_flag(F,frame,f,u0,v0,s0):
    """Test one transverse direction, or the tangent when s0 is None."""
    K,alpha=data_field(f)
    U,V=sp.symbols('U V')
    def canon(value):
        return K.to_sympy(K.from_sympy(sp.sympify(value)))
    def conj_value(expr):
        return canon(sp.sympify(expr).subs(T,alpha))
    u,v=conj_value(u0),conj_value(v0)
    x,y,n=[fa.COORDS[i] for i in frame]
    aff=F.as_expr().subs(n,1)
    at=lambda expr:canon(expr.subs({x:u,y:v},simultaneous=True))
    if at(aff)!=0:return []
    A=at(aff.diff(y))
    if A==0:return []
    d=canon(-at(aff.diff(x))/A)
    adapted=sp.Poly(aff.subs({x:U+u,y:V+d*U+v},simultaneous=True),U,V,domain=K)
    coeff=adapted.coeff_monomial
    coeffs=adapted.terms()
    def G(f1,f2,f3,divide):
        raw=sum(q*f1**i*f2**j*f3**(4-i-j) for (i,j),q in coeffs)
        poly=sp.Poly(raw,t,z,domain=K)
        return poly.exquo(sp.Poly(t**divide,t,z,domain=K)) if divide else poly
    rows=[]
    if s0 is None:
        B=coeff((2,0))
        if B==0:return []
        a=canon(4*B/(3*A))
        b=canon(4*coeff((1,1))/(3*A)-3*coeff((3,0))/(2*B))
        rows.append(('II',{'a':a,'b':b},G(t,1,-a*t*t-b*t-z,0)))
    else:
        s=conj_value(s0)
        h=sp.Poly(sum(q*s**i*(-1)**j*cvar**(4-i-j) for (i,j),q in coeffs),cvar,domain=K)
        jrow=sp.Poly(sum(-i*q*s**(i-1)*(-1)**j*cvar**(4-i-j) for (i,j),q in coeffs if i),cvar,domain=K)
        hp=h.diff()
        hpp=hp.diff()
        D=2*h*hpp-hp**2
        for c in D.ground_roots():
            hh=h.eval(c)
            if hh==0:continue
            j1=jrow.diff()
            b=canon(-(2*hh**2*j1.diff().eval(c)-2*hh*hp.eval(c)*j1.eval(c)+hp.eval(c)**2*jrow.eval(c))/
                    (2*hh**2*hpp.diff().eval(c)))
            rows.append(('I',{'s':s,'b':b,'c':c},G(-t+s*t*t,-t*t,z+b*t+c*t*t,2)))
    result=[]
    for kind,params,poly in rows:
        m=poly.coeff_monomial
        g=m((6,0))
        if g==0:continue
        e2=canon(m((5,0))/(2*g))
        f1=canon(m((4,1))/(2*g))
        e1=canon((m((4,0))/g-e2**2)/2)
        e0=canon(m((3,0))/(2*g)-e2*e1)
        f0=canon(m((3,1))/(2*g)-e2*f1)
        shape=t**3+e2*t*t+(e1+f1*z)*t+e0+f0*z
        residual=poly-sp.Poly(g*shape**2,t,z,domain=K)
        if any(i for (i,j),q in residual.terms()):continue
        values={'u':u,'v':v,'d':d,'hinv':canon(1/(A*g)),'g':g,'e2':e2,'e1':e1,
                'f1':f1,'e0':e0,'f0':f0,**params}
        field_poly=sp.Poly(f,T,domain=sp.QQ)
        coords={k:str(rm.coefficient_as_polynomial(K.from_sympy(value),K)) for k,value in values.items()}
        length=field_poly.degree()
        # This is a parameterization of a verified point (possibly over a
        # nonminimal field), not a length assertion about the whole chart.
        record={'dimension':0,'points':length,'length':length,
                'minimal_polynomial':str(field_poly.as_expr()),'coordinates':coords}
        result.extend(rm.extract_chart(F,frame,kind,record))
    return result


def exhaustive_ticks(F):
    """Every turn does finite work and yields, including failed trials."""
    import itertools
    for height in itertools.count(1):
        for d in range(1,min(1680,height)+1):
            polynomials=([T] if d==1 else
                (T**d+sum(a*T**i for i,a in enumerate(c)) for c in product(range(-height,height+1),repeat=d)))
            for f in polynomials:
                if not sp.Poly(f,T,domain=sp.QQ).is_irreducible:
                    yield []
                    continue
                for frame in fa.FRAMES:
                    x,y,n=[fa.COORDS[i] for i in frame]
                    aff=F.as_expr().subs(n,1)
                    for u in elements(d,height):
                        for v in elements(d,height):
                            at=lambda expr:sp.Poly(expr.subs({x:u,y:v},simultaneous=True),T,domain=sp.QQ).rem(sp.Poly(f,T,domain=sp.QQ))
                            if not at(aff).is_zero or at(aff.diff(y)).is_zero:
                                yield []
                                continue
                            yield scan_flag(F,frame,f,u,v,None)
                            for s in elements(d,height):
                                yield scan_flag(F,frame,f,u,v,s)
