"""Verify an elliptic generator and exact quadratic-discriminant witnesses for C6."""
import argparse
from fractions import Fraction
from math import isqrt
from cypari2 import Pari
pari=Pari()
pari.allocatemem(64000000,silent=True)


def witness(m):
    if not isinstance(m,int) or m<=0 or m%2==0:
        raise ValueError('Use a positive odd multiple.')
    E=pari.ellinit([0,21,0,441,0])
    P=pari([84,882])
    if E.ellmul(P,2)!=pari(['225/16','7335/64']):
        raise ArithmeticError('Incorrect doubling identity')
    X,Y=E.ellmul(P,m)
    u=X/21
    nr,dr=int(u.numerator()),int(u.denominator())
    a,b=isqrt(nr),isqrt(dr)
    if a*a!=nr or b*b!=dr or a==b:
        raise ArithmeticError('Expected a nonbranch square abscissa')
    z=Y/(21*pari(a)/b)
    if z*z!=21*((pari(a)/b)**4+(pari(a)/b)**2+1):
        raise ArithmeticError('Quartic lift failed')
    delta_bound=84*abs(a*a-b*b)
    height=max(a,b)
    return height,delta_bound,a>b


if __name__=='__main__':
    app=argparse.ArgumentParser(description=__doc__)
    app.add_argument('odd_multiple',nargs='?',type=int,default=213)
    args=app.parse_args()
    height,disc_bound,real=witness(args.odd_multiple)
    ratio_squared=Fraction(height*height,disc_bound)
    floor_hundredths=isqrt(10000*ratio_squared.numerator//ratio_squared.denominator)
    if args.odd_multiple==213 and (not real or floor_hundredths!=153):
        raise ArithmeticError('Example ratio check failed')
    print('Odd multiple',args.odd_multiple,'real quadratic field:',real)
    print('Certified ratio lower bound:',str(Fraction(floor_hundredths,100)))
    print('Abscissa height in bits:',height.bit_length())
