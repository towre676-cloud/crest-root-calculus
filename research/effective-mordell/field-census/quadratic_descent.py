"""Split the rational-abscissa twist into rational squareclasses and test elliptic fibres."""
from fractions import Fraction
from functools import lru_cache
from math import isqrt
from cypari2 import Pari
from eisenstein_quadratic import inert_supported, solve_residual_pair
from prime_density import certify_prime

p=Pari()

def sqroot(x):
    x=Fraction(x)
    if x<0:return None
    a,b=isqrt(x.numerator),isqrt(x.denominator)
    return Fraction(a,b) if a*a==x.numerator and b*b==x.denominator else None

def sf(n):
    acc=-1 if n<0 else 1
    f=p.factorint(abs(n))
    for q,e in zip(f[0],f[1]):
        if not p.isprime(q):raise ArithmeticError('unproved squareclass factor')
        if int(e)%2:acc*=int(q)
    return acc

@lru_cache(None)
def square_mod(m):return frozenset(t*t%m for t in range(m))

@lru_cache(None)
def forms(q,m):
    pairs=[(t,1) for t in range(m)]+[(1,t) for t in range(0,m,q)]
    return [((u*u+u*v+v*v)%m,(u*u-u*v+v*v)%m,(u*u-v*v)%m) for u,v in pairs]

def obstruction(d,a,b):
    for q,m in [(2,64),(3,81),(5,25),(7,49),(11,121),(13,169),(17,17),(19,19),
                (23,23),(29,29),(31,31),(37,37),(41,41),(43,43)]:
        good=square_mod(m)
        if not any(a*A%m in good and b*B%m in good and d*a*b*D%m in good for A,B,D in forms(q,m)):
            return m
    return None

@lru_cache(None)
def elliptic_fibres(kind,k):
    if kind=='even':
        E=p.ellinit([0,k,0,k*k,0])
    else:E=p.ellinit([0,9*k,0,30*k*k,36*k**3])
    if int(E.ellrank()[1])>0:return None
    count,orders,gens=E.elltors();pts=[p([0])]
    for n,g in zip(orders,gens):
        pts=[E.elladd(T,E.ellmul(g,j)) for T in pts for j in range(int(n))]
    if len({str(v) for v in pts})!=int(count):raise ArithmeticError('torsion enumeration')
    xs=set()
    for P in pts:
        if len(P)==1:continue
        X=Fraction(str(P[0]))
        if kind=='even':
            u=sqroot(X/k)
            if u is not None:xs.update((u,-u))
        elif X:
            xs.add(1+6*k/X)
    return tuple(sorted(xs))

def classify_twist(d):
    if inert_supported(d) or certify_prime(abs(d)):
        return set(),[]
    divisors=[1]
    for q in p.factorint(abs(3*d))[0]:
        if not p.isprime(q):raise ArithmeticError('unproved local support prime')
        j=int(q)
        if j!=2:divisors += [j*b for b in divisors]
    xs=set();open_pairs=[]
    for a in divisors:
        for b in divisors:
            if obstruction(d,a,b):continue
            for kind,k,sign in [('even',sf(a*b),1),('mixed',sf(d*b),1),('mixed',sf(d*a),-1)]:
                case=elliptic_fibres(kind,k)
                if case is None:continue
                for r in case:
                    X=sign*r
                    if X not in (0,-1,1) and sqroot(d*(X**6-1)) is not None:xs.add(X)
                break
            else:
                exact=solve_residual_pair(d,a,b)
                if exact is None:
                    open_pairs.append((a,b))
                else:
                    xs.update(exact)
    return xs,open_pairs
