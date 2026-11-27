"""Complete C6(K) when rank E_-(K) <= 1 is certified and i belongs to K.

Degree two and four use short Gaussian height lists and exact quadratic divisors.
In higher degree the finite fallback enumerates primitive-polynomial coefficient
boxes over Q(i); it is a proof of effectivity and may be very expensive.
"""
import argparse
from functools import lru_cache
import json
from itertools import product
from math import comb, gcd, isqrt, lcm
from pathlib import Path
from classgroup_descent import classify as selmer
from isogeny_over_field import pari
from rank_one_height import check, sixth, gaussian_mul, verify_identities


def contains_gaussian_square(z):
    a,b=z
    n=a*a+b*b
    h=isqrt(n)
    if h*h!=n or (h+a)%2:
        return False
    u=isqrt((h+a)//2)
    v=isqrt((h-a)//2)
    return 2*u*u==h+a and 2*v*v==h-a and 2*u*v==abs(b)


def disk(radius,positive=False):
    for a in range(-radius,radius+1):
        width=isqrt(radius*radius-a*a)
        for b in range(-width,width+1):
            if positive and not (a>0 and b>=0):
                continue
            yield a,b


@lru_cache(None)
def small_gaussians():
    ball=[(a,b) for a,b in disk(7) if a*a+b*b<=44]
    seen=set()
    for (a,b),(c,d) in product(ball,ball):
        if c<=0 or d<0:
            continue
        p,q,r=a*c+b*d,b*c-a*d,c*c+d*d
        g=gcd(gcd(p,q),r)
        seen.add((p//g,q//g,r//g))
    check(len(seen)==3165,"Gaussian height box")
    return sorted(seen)


def relative_radicand(nf,i,X):
    f=nf.nf_get_pol()
    w=f.variable()
    alpha=pari.Mod(w,f)
    conjugates=pari.nfroots(nf,pari.subst(f,w,X))
    matches=[z for z in conjugates if z!=alpha and pari.subst(pari.lift(i),w,z)==i]
    check(len(matches)==1,"Relative quadratic conjugate")
    delta=(alpha-matches[0])**2
    real=pari.nfelttrace(nf,delta)/4
    imag=-pari.nfelttrace(nf,i*delta)/4
    check(delta==real+i*imag,"Relative quadratic radicand")
    den=lcm(int(real.denominator()),int(imag.denominator()))
    return int(den**2*real),int(den**2*imag)


def points_low_degree(nf,i,degree):
    X=pari('x')
    points=set()
    def put(x,y):
        check(y*y==x**6-1,"Incorrect K-fibre")
        points.add((str(pari.lift(x)),str(pari.lift(y))))
    delta=None if degree==2 else relative_radicand(nf,i,X)
    for a,b,c in small_gaussians():
        real,imag=sixth((a,b))
        z=real-c**6,imag
        if not (contains_gaussian_square(z) or
                (delta is not None and contains_gaussian_square(gaussian_mul(z,delta)))):
            continue
        x=(a+i*b)/c
        for y in pari.nfroots(nf,X**2-(x**6-1)):
            put(x,y)
    if degree==4:
        records=json.loads(Path(__file__).with_name('quadratic-pairs-over-gaussian.json').read_text())
        def gauss(c):
            return pari(c[0])+i*pari(c[1])
        for row in records:
            U=pari.Polrev([gauss(c) for c in row['u']],X)
            V=pari.Polrev([gauss(c) for c in row['v']],X) if row['v'] else 0
            for z in pari.nfroots(nf,U):
                put(1+1/z,pari.subst(V,X,z)/z**3)
    return sorted(points)


def gaussian_polynomials(n,i):
    """Polynomial boxes from the relative Mahler measure; streamed in bounded space."""
    X=pari('x')
    measure=7**n
    def fill(j,coeffs):
        if j<0:
            yield pari.Polrev(coeffs,X)
            return
        for a,b in disk(comb(n,j)*measure,positive=(j==n)):
            coeffs[j]=a+b*i
            yield from fill(j-1,coeffs)
    yield from fill(n,[pari(0)]*(n+1))


def points_general(nf,i,degree):
    X=pari('x')
    points=set()
    relative_degree=degree//2
    for n in range(1,relative_degree+1):
        if relative_degree%n:
            continue
        for polynomial in gaussian_polynomials(n,i):
            for x in pari.nfroots(nf,polynomial):
                for y in pari.nfroots(nf,X**2-(x**6-1)):
                    check(y*y==x**6-1,"Incorrect K-fibre")
                    points.add((str(pari.lift(x)),str(pari.lift(y))))
    return sorted(points)


def classify(polynomial, bounded_only=False):
    if (not isinstance(polynomial,list) or not polynomial or polynomial[-1]!=1
            or any(type(c)!=int for c in polynomial) or len(polynomial)<2):
        raise ValueError("Supply an ascending monic integer polynomial")
    f=pari.Polrev(polynomial,pari('w'))
    if not pari.polisirreducible(f):
        raise ValueError("The defining polynomial must be irreducible")
    degree=len(polynomial)-1
    nf=pari.nfinit(f)
    units=pari.nfroots(nf,pari('x^2+1'))
    if not len(units):
        return dict(field=str(f),status='partial',reason='This entry point needs a field containing i.')
    if bounded_only:
        points=(points_low_degree if degree<=4 else points_general)(nf,units[0],degree)
        return dict(field=str(f),status='complete_if_elliptic_rank_at_most_one',
                    verified_affine_points=points,rational_infinities=2,
                    method='finite field enumeration without a rank or BNF computation')
    rank=selmer(polynomial)
    if rank['status']=='complete':
        return rank
    bound=rank['rank_upper_bounds']['-1']
    if bound>1:
        return dict(field=str(f),status='partial',rank_upper_bound=bound,
                    reason='The rank-one hypothesis has not been certified.')
    points=(points_low_degree if degree<=4 else points_general)(nf,units[0],degree)
    return dict(field=str(f),status='complete',rank_upper_bound=bound,
                verified_affine_points=points,rational_infinities=2,
                method='uniform rank-one abscissa height bound')


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('ascending_polynomial')
    cli.add_argument('--bounded-only',action='store_true',
                     help='Perform finite enumeration directly and report the rank hypothesis as conditional.')
    args=cli.parse_args()
    verify_identities()
    print(json.dumps(classify(json.loads(args.ascending_polynomial),args.bounded_only),indent=2))
