"""Exact rational and Eisenstein squareclass sieves in the sqrt(7) proof."""
from functools import lru_cache
from itertools import product
from cypari2 import Pari
import sympy as sp

p=Pari()
def require(c,msg):
    if not c:raise ArithmeticError(msg)

def rational_filters():
    for k,size in ((1,2),(3,4),(7,2)):
        E=p.ellinit([0,k,0,k*k,0])
        require(list(E.ellrank()[:2])==[0,0] and int(E.elltors()[0])==size,'even quartic')
        _,invariants,gens=E.elltors()
        parts=[p([0])]
        for n,g in zip(invariants,gens):
            parts=[E.elladd(t,E.ellmul(g,j)) for t in parts for j in range(int(n))]
        for t in parts:
            if len(t)==1:continue
            xx=t[0]/k
            if p.issquare(xx):require(xx in [0,1],'unexpected rational abscissa')
    k=7
    E=p.ellinit([0,9*k,0,30*k*k,36*k**3])
    require(list(E.ellrank()[:2])==[0,0] and int(E.elltors()[0])==2 and
            E.ellisoncurve([-3*k,0]),'mixed quartic')
    x,r=sp.symbols('x r')
    t=(3+r*r)/(3-r*r)
    require(sp.cancel((t*t+t+1)*(3-r*r)**2-(r**4+27))==0,'first parameterization')
    require(sp.cancel((t*t-t+1)*(3-r*r)**2-3*(r**4+3))==0,'second parameterization')

def mul(x,y):
    a,b=x;c,d=y
    return a*c-b*d,a*d+b*c-b*d

def scale(v,p):return tuple(p*c for c in v)
def mod(v,m):return tuple(c%m for c in v)
@lru_cache(None)
def squares(m):return {mod(mul(a,a),m) for a in product(range(m),repeat=2)}
def local(delta,eta,p,m):
    square={a*a%m for a in range(m)}
    for u,v in [(j,1) for j in range(m)]+[(1,j) for j in range(0,m,p)]:
        if (u**4+3*v**4)%m not in square or (7*(u**4+27*v**4))%m not in square:continue
        beta=(u*u+v*v,2*v*v);gamma=(u*u+3*v*v,6*v*v)
        if mod(mul(delta,beta),m) in squares(m) and mod(mul(eta,gamma),m) in squares(m):
            return True
    return False

def eisenstein():
    # Basis (1,omega), omega^2+omega+1=0. The two displayed primes have norm 7.
    candidates=[((q,0),scale(t,j)) for q in (1,-1,2,-2) for t in [(3,2),(1,-2)]
                for j in (1,-1,2,-2)]
    for prime,m in ((2,16),(3,9)):
        candidates=[(d,e) for d,e in candidates if local(d,e,prime,m)]
        print('After modulus',m,candidates)
    require(candidates==[((-2,0),(2,-4))],'squareclass sieve')
    # y=3*sqrt(7); all four sign choices at x=+/-2.
    require(2**6-1==7*3**2,'quadratic point')
if __name__=='__main__':
    rational_filters()
    eisenstein()
