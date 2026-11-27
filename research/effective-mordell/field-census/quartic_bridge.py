"""Verify two reconstructed quartic maps and obtain an absolute height bound from them."""
from fractions import Fraction
import json
from math import lcm
from pathlib import Path
from cypari2 import Pari

pari = Pari()
x, w = pari('x'), pari('w')
root = pari.Mod(w, w*w-w+1)


def terms(poly):
    return [pari.polcoef(poly,i,x) for i in range(int(pari.poldegree(poly,x))+1)] if poly else []


def coefficients(value):
    v = pari.lift(value)
    return [Fraction(str(pari.polcoef(v,i,w))) for i in (0,1)]


def clear(polynomials):
    scale = lcm(*(a.denominator for f in polynomials for c in terms(f) for a in coefficients(c)))
    return [scale*f for f in polynomials], scale


def size(poly):
    # Both conjugates of w have complex absolute value one.
    return sum(abs(a) for c in terms(poly) for a in coefficients(c))


def bounds(r):
    (N,D),_ = clear([pari.numerator(r),pari.denominator(r)])
    degree = max(len(terms(N)),len(terms(D)))-1
    def reciprocal(f):
        return sum(c*x**(degree-i) for i,c in enumerate(terms(f)))
    bezout = []
    for f,g in ((N,D),(reciprocal(N),reciprocal(D))):
        a,b,h = pari.gcdext(f,g)
        if not h or int(pari.poldegree(h,x))>0:
            raise ArithmeticError('Common projective zero')
        a,b = a/h,b/h
        if a*f+b*g!=1:
            raise ArithmeticError('Incorrect Bezout coefficients')
        bezout.extend([a,b])
    exact,multiplier = clear(bezout)
    for offset,(f,g) in enumerate(((N,D),(reciprocal(N),reciprocal(D)))):
        if exact[2*offset]*f+exact[2*offset+1]*g!=multiplier:
            raise ArithmeticError('Integral Bezout certificate failed')
    upper=max(size(N),size(D))
    lower=max(size(exact[0])+size(exact[1]),size(exact[2])+size(exact[3]))
    return dict(degree=degree,forward=int(upper),inverse=int(lower))


def certificate():
    receipt=json.loads(Path(__file__).with_name('quartic-covers.json').read_text())
    source=(x**6+1)*root**0
    maps=[]
    forms=[]
    for row in receipt['maps'][:2]:
        def parse(key):
            return pari.subst(pari(row[key]),w,root)
        S=parse('half_poles')
        X=parse('numerator')/S**2
        V=parse('lift')/S**3
        if X**3-15*X+22!=source*V**2:
            raise ArithmeticError('Quartic map identity')
        forms.append(pari.deriv(X,x)/V)
        u=X-2
        X0=(u+2-3/u)/4
        V0=V*(1+3/u**2)/8
        if X0**3+1!=source*V0**2:
            raise ArithmeticError('Two-isogeny identity')
        maps.append((X0,V0))
    if not pari.deriv(forms[0]/forms[1],x):
        raise ArithmeticError('Dependent differential forms')
    singles=[bounds(a) for a,_ in maps]
    signed=[]
    for j in range(3):
        for sign in (-1,1):
            a,b=maps[0]
            c,e=maps[1]
            c,e=root**(2*j)*c,sign*e
            abscissa=source*((b-e)/(a-c))**2-a-c
            signed.append(dict(rotation=j,sign=sign,**bounds(abscissa)))
    dup=bounds((x**4-8*x)/(4*x**3+4)*root**0)
    degrees=[s['degree'] for s in signed]
    if [s['degree'] for s in singles]!=[8,8] or degrees!=[16,16,4,28,28,4] or dup['degree']!=4:
        raise ArithmeticError('Unexpected degrees')
    M=max(r['forward'] for r in singles)
    m=max(r['inverse'] for r in singles)
    # 265/153 is strictly smaller than sqrt(3). Each signed map keeps its own degree.
    h,k=265,153
    if h*h>=3*k*k:
        raise ArithmeticError('Invalid rational angular estimate')
    for row in signed:
        exponent=3*(k*row['degree']-(2*k-h)*8)
        ceiling=M**(3*k)*dup['forward']**k*m**(3*(h-k))*dup['inverse']**h*row['inverse']**(3*k)
        if exponent<=0 or 48659**exponent<=ceiling:
            raise ArithmeticError('Claimed field-independent height ceiling failed')
    direct=[]
    omega=root**2
    for sign in (-1,1):
        def addition(a,c,b,e):
            return source*((b-e)/(a-c))**2-a-c
        base=addition(x*x,x**-2,1,sign*x**-3)
        simple=bounds(base)
        if simple['degree']!=4 or simple['inverse']!=36:
            raise ArithmeticError('Reciprocal map coefficient bound')
        for j in range(3):
            scale=omega**j
            rotated=addition(x*x,scale**(-4)*x**-2,1,sign*x**-3)
            if pari.subst(base,x,scale*x)!=scale**2*rotated:
                raise ArithmeticError('Failed unit-rotation identity')
        direct.append(simple)
    return dict(target='Y^2=X^3+1',maps=singles,signed_sums=signed,duplication=dup,
                quartic_reconstruction_strict_height_bound=48659,reciprocal_sums=direct)

if __name__=='__main__':
    print(json.dumps(certificate(),indent=2))
