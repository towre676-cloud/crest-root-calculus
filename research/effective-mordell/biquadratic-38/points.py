"""Return the complete projective point list on y^2=x^6-1 over Q(i,sqrt(7))."""
import argparse
import json
from pathlib import Path
from cypari2 import Pari

p=Pari()
x,w=p('x'),p('w')
DEFAULT=[64,0,-12,0,1]


def check(test,detail):
    if not test:
        raise ArithmeticError(detail)


def classify(coefficients):
    if not isinstance(coefficients,list) or len(coefficients)!=5 or any(type(c)!=int for c in coefficients) or coefficients[-1]!=1:
        raise ValueError('Supply the ascending coefficients of a monic integral quartic.')
    f=p.Polrev(coefficients,w)
    if not p.polisirreducible(f):
        raise ValueError('The polynomial must be irreducible.')
    nf=p.nfinit(f)
    imag=p.nfroots(nf,x*x+1)
    radical=p.nfroots(nf,x*x-7)
    if len(imag)!=2 or len(radical)!=2:
        raise ValueError('The field must contain i and the square root of 7.')
    i,r=imag[0],radical[0]
    pts=set()
    def add(a,b):
        check(b*b==a**6-1,'Point does not satisfy the equation')
        pts.add((str(p.lift(a)),str(p.lift(b))))
    for sign in (-1,1):
        add(sign,0)
        add(0,sign*i)
        for v in (-1,1):
            add(2*sign,3*v*r)
            add(sign*p(1)/2,3*v*i*r/8)
    check(len(pts)==12,'Rational abscissas')
    src=Path(__file__).parent/'prior'/'quadratic-pairs-over-gaussian.json'
    records=[row for row in json.loads(src.read_text()) if row['extension']==7]
    check(len(records)==12,'Conjugate divisor count')
    for row in records:
        def poly(co):
            return p.Polrev([p(a)+i*p(b) for a,b in co],x) if co else p(0)
        U,V=poly(row['u']),poly(row['v'])
        roots=p.nfroots(nf,U)
        check(len(roots)==2,'Quadratic support failed to split')
        for t in roots:
            add(1+1/t,p.subst(V,x,t)/t**3)
    check(len(pts)==36,'Wrong number of affine points')
    return dict(field=str(f),degree=4,quadratic_subfields=[-7,-1,7],status='complete',
                method='Gaussian descent and relative elliptic Chabauty at 37',
                verified_affine_points=sorted(pts),rational_infinities=2,
                projective_points=38,rational_abscissa_including_infinities=14,
                outside_Gaussian_abscissa=24)


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('ascending_polynomial',nargs='?',default=json.dumps(DEFAULT))
    print(json.dumps(classify(json.loads(cli.parse_args().ascending_polynomial)),indent=2))
