"""Exact geometric and rational elliptic-map counts in the cyclic quartic family."""
import argparse
import json
from pathlib import Path
import sys
from math import gcd, isqrt

import sympy as sp

from riemann_roch import X, Y, Z

E, x = sp.symbols('e x')
A, B = sp.symbols('a b')
Q = A*A+2*B*B
G = A**3-6*A*B*B
PARAMETER_MODULUS = 24871
PRIMES = (7,11,17,19)
ISOGENY_PRIMES = (2,3,5,7,11,13,17,19,37,43,67,163)


def norm_count(m):
    """Primitive Eisenstein vectors of norm m, up to simultaneous negation."""
    factors = sp.factorint(m)
    if factors.get(3,0)>1 or any(p%3==2 for p in factors):
        return 0
    return 3*2**sum(p%3==1 for p in factors)


def classes(degree):
    if type(degree) is not int or degree < 1:
        raise ValueError('Choose a positive map degree')
    geometric = int(sp.divisor_sigma(degree//3)) if degree%3==0 else 0
    if degree%2==0:
        for m in sp.divisors(degree//2):
            geometric += norm_count(m)*int(sp.divisor_sigma(degree//(2*m)))
    rational = [d for d in (2,6,12) if degree%d==0 and isqrt(degree//d)**2==degree//d]
    return {'degree':degree,'geometric_primitive':int(degree==3)+(norm_count(degree//2) if degree%2==0 else 0),
            'geometric_all':geometric,'rational_all':len(rational),'rational_base_degrees':rational}


def models(e):
    I=4*(16*e**2-54*e+81)
    J=-8*(128*e**3-675*e**2+1458*e-1458)
    return ((-27*I,-27*J),(261,1494-1024*e))


def trace(p,Acoef,Bcoef):
    if (4*Acoef**3+27*Bcoef**2)%p==0:
        raise ArithmeticError('Bad reduction in elliptic isogeny certificate')
    return -sum(int(sp.legendre_symbol((t**3+Acoef*t+Bcoef)%p,p)) for t in range(p))


def direct_source_counts(p, e, degree):
    """Count projective plane points by solving two successive quadratics over F_(p^degree)."""
    import numpy as np
    sys.path.insert(0,str(Path(__file__).resolve().parents[1]/'cubic-tower-factors'))
    from field import Field
    field = Field(p,degree)
    ag = field.log
    def times(a,b):
        return np.where((a==-1)|(b==-1),-1,(a+b)%field.order)
    def evaluate(coefficients):
        value = ag[coefficients[-1]%p]
        for c in reversed(coefficients[:-1]):
            value = field.add(times(value,ag),ag[c%p])
        return value
    # Let v=b^2. The affine equation is 4v^2 + linear(a)*v + constant(a)=0.
    linear = evaluate([6,-6,4])
    constant = evaluate([e,0,3,1,1])
    delta = field.add(times(linear,linear),times(ag[-16%p],constant))
    sqrt = np.where(delta==-1,-1,delta//2)
    minus_linear = times(ag[-1%p],linear)
    def roots(v):
        return np.where(v==-1,1,np.where(v%2==0,2,0))
    def branch(sign):
        numerator = field.add(minus_linear,times(ag[sign%p],sqrt))
        return roots(times(numerator,ag[pow(8,-1,p)]))
    counts = np.where(delta==-1,branch(1),np.where(delta%2==0,branch(1)+branch(-1),0))
    return int(counts.sum()+roots(ag[-2%p]))


def verify_family():
    e0=10
    require(PARAMETER_MODULUS==sp.prod(PRIMES),'Parameter modulus differs')
    # A smooth special fibre proves smoothness for the entire integral progression.
    q=X*X+2*Y*Y;g=X**3-6*X*Y*Y;source=q*q+Z*g+3*Z*Z*q+e0*Z**4
    for idx in (X,Y,Z):
        rest=[a for a in (X,Y,Z) if a!=idx]
        basis=sp.groebner([source.diff(v).subs(idx,1) for v in (X,Y,Z)],*rest,modulus=11)
        require(list(basis)==[1],'Singular source modulo 11')
    # Eliminate e from the generic source. This verifies the maps on every smooth fibre.
    parameter=-Q**2-G-3*Q
    t=Q;y=2*B*(3*A*A-2*B*B)
    f=-2*(x*x+3*x+E)**2+2*x**3
    require(sp.expand(y*y-f.subs({x:t,E:parameter}))==0,'Cyclic quotient equation failed')
    sqrt6=sp.sqrt(6)
    rot={A:(-A+sqrt6*B)/2,B:-B/2-sqrt6*A/4}
    inv={A:(-A-sqrt6*B)/2,B:-B/2+sqrt6*A/4}
    for u in (A,B):
        require(sp.expand(rot[u].xreplace(rot)-inv[u])==0 and sp.expand(rot[u].xreplace(inv)-u)==0,'Deck action failed')
    require(sp.expand(Q.xreplace(rot)-Q)==0 and sp.expand(G.xreplace(rot)-G)==0,'Cyclic invariance failed')
    EA,EB=models(E)
    ax,bx,cx,dx,ex=sp.Poly(f,x).all_coeffs()
    i=12*ax*ex-3*bx*dx+cx*cx;j=72*ax*cx*ex+9*bx*cx*dx-27*ax*dx*dx-27*bx*bx*ex-2*cx**3
    require(sp.expand(EA[0]+27*i)==0 and sp.expand(EA[1]+27*j)==0,'Incorrect Jacobian of the quartic target')
    # The rational reflection and the difference of its conjugates, on the quadratic twist.
    R0=3-16*A;R1=16*(2*Q-3*A+3)
    aR,bR=EB
    require(sp.expand(R1**2-R0**3-aR*R0-bR.subs(E,parameter))==0,'Reflection quotient failed')
    x0=3+8*A;u=(4*Q+3*A+6)/B
    x6=u**2-12*x0;y6=864*B+u*(x6-6*x0)
    require(sp.cancel(y6**2-x6**3-36*aR*x6-216*bR.subs(E,parameter))==0,'Twisted degree-six map failed')
    xa,ya=R0.xreplace(rot),R1.xreplace(rot)
    xb,yb=R0.xreplace(inv),R1.xreplace(inv)
    slope=sp.cancel((ya+yb)/(xa-xb))
    delta_x=sp.cancel(slope**2-xa-xb)
    delta_y=sp.cancel(slope*(xa-delta_x)-ya)
    require(sp.cancel(6*delta_x-x6)==0 and sp.cancel(6*sqrt6*delta_y-y6)==0,'Difference of conjugate maps differs')
    reductions=[]
    for coeffs in models(e0):
        row=[trace(p,*coeffs) for p in PRIMES]
        for ell in ISOGENY_PRIMES:
            require(any(p!=ell and all((u*u-a*u+p)%ell for u in range(ell)) for a,p in zip(row,PRIMES)),
                    'A prime rational isogeny survived the sieve')
        reductions.append(row)
    require(reductions==[[5,4,-2,5],[-2,5,5,-4]],'Elliptic traces changed')
    curve_counts = []
    for i,p in enumerate(PRIMES):
        row = [direct_source_counts(p,e0,n) for n in (1,2,3)]
        traces = [reductions[0][i],reductions[1][i],int(sp.legendre_symbol(6,p))*reductions[1][i]]
        predicted = [p+1-sum(traces), p**2+1-sum(a*a-2*p for a in traces),
                     p**3+1-sum(a**3-3*p*a for a in traces)]
        require(row==predicted,'Direct projective counts disagree with the three elliptic factors')
        curve_counts.append(row)
    require(sp.legendre_symbol(6,7)==-1,'Twist comparison prime failed')
    for index in (0,1):
        a,b=reductions[index][:2]
        require(a%7 and b%11 and not sp.sqrt(sp.Rational(a*a-28,b*b-44)).is_rational,
                'The two ordinary reductions failed to rule out CM')
    require(not sp.sqrt(sp.Rational(reductions[0][0]**2-28,reductions[1][0]**2-28)).is_rational,
            'The simple directions were not geometrically separated')
    for a,b in ((a,b) for a in range(-9,10) for b in range(-9,10) if gcd(a,b)==1):
        bound=gcd(2*a-b,2*b-a)
        require(bound in (1,3) and (bound!=3 or (a%3 and (b-2*a)%3==0)), 'Primitivity congruence differs')
    for m in range(1,101):
        lim=isqrt(2*m)+2
        represented=sum(a*a-a*b+b*b==m and gcd(a,b)==1 for a in range(-lim,lim+1) for b in range(-lim,lim+1))//2
        require(represented==norm_count(m),'Primitive norm representation formula differs')
    return {'parameter_residue':e0,'parameter_modulus':PARAMETER_MODULUS,'good_primes':list(PRIMES),
            'cyclic_jacobian_traces':reductions[0],'reflection_target_traces':reductions[1],
            'quadratic_twist':6,'source_counts_degrees_one_two_three':curve_counts,'rational_base_degrees':[2,6,12],
            'geometric_classes_through_12':[classes(n) for n in range(1,13)]}


def require(test, message):
    if not test:
        raise ArithmeticError(message)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('degree',nargs='?',type=int)
    parser.add_argument('--verify',action='store_true')
    args=parser.parse_args()
    if args.verify:
        answer=verify_family()
    elif args.degree:
        answer=classes(args.degree)
    else:
        parser.error('Supply a positive degree or --verify')
    print(json.dumps(answer,indent=2))
