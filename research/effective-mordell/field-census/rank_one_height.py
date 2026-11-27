"""Uniform generator-free height bound and a practical biquadratic panel for C6."""
import argparse
from collections import Counter
from functools import lru_cache
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path
from cypari2 import Pari

pari=Pari()

def check(condition,message):
    if not condition:
        raise ArithmeticError(message)

def verify_identities():
    x=pari('x')
    j=pari.Mod(pari('w'),pari('w^2+1'))
    f=x**4+8*x
    g=4*x**3-4
    check(9*x**7==(9*x**3-8)*f-16*x*g,"First duplication identity")
    check(36==4*x**2*f-(x**3+9)*g,"Second duplication identity")
    for u in (j,-j):
        n=-2*u*x**3-3*x*x+2*u*x
        d=(x*x-1)*(x-u)**2
        # The ordinate of Q is i*y/x^3 and u records its sign.
        raw=(x**6-1)*(u-x**3)**2-(x**4+1)*(1-x**4)**2
        check(raw*d==n*x*x*(1-x**4)**2,"Reduced abscissa identity")
        a=4*u*x**3+5*x*x-4*u*x-6
        b=-8*x*x+6*u*x+3
        check(a*n+b*d==3,"Z^7 quartic identity")
        a=-6*x**3+4*u*x*x+5*x-4*u
        b=3*x**3-6*u*x*x-8*x
        check(a*n+b*d==3*x**7,"X^7 quartic identity")
    # Absolute projective H(x)^12 <= 36^3*156816 < 45^6.
    check(36**3*156816<45**6,"Gaussian component norm bound")


def gaussian_mul(x,y):
    a,b=x
    c,d=y
    return (a*c-b*d, a*d+b*c)


def cube(x):
    return gaussian_mul(gaussian_mul(x,x),x)


def sixth(x):
    c=cube(x)
    return gaussian_mul(c,c)


def squared_norm(z):
    return z[0]**2+z[1]**2


def rational_squareclass(z):
    a,b=z
    if not a and not b:
        return 0
    t=isqrt(a*a+b*b)
    if t*t!=a*a+b*b:
        return None
    representative=2*(a+t) if a+t else abs(a)
    return int(pari.core(representative))


@lru_cache(maxsize=1)
def gaussian_candidates():
    """All Gaussian abscissas under the proved height bound, with rational twist class."""
    ball=[(a,b) for a,b in product(range(-6,7),repeat=2) if a*a+b*b<=44]
    # Denominators in a single quadrant represent all associates.
    denoms=[(a,b) for a,b in ball if a>0 and b>=0]
    seen=set()
    values=[]
    for num,den in product(ball,denoms):
        a,b=num
        c,d=den
        A=a*c+b*d
        B=b*c-a*d
        C=c*c+d*d
        common=gcd(gcd(A,B),C)
        key=(A//common,B//common,C//common)
        if key in seen:
            continue
        seen.add(key)
        U,V=sixth(num),sixth(den)
        z=U[0]-V[0],U[1]-V[1]
        squareclass=rational_squareclass(z)
        if squareclass is not None:
            values.append((key,squareclass))
    return sorted(values),len(seen)


def verify_non_gaussian_table():
    from quadratic_pairs import verify
    verify()


def classify_biquadratic(d):
    if d<=1 or not pari.issquarefree(d):
        raise ValueError("Supply a positive squarefree integer greater than one.")
    ranks=[list(map(int,pari.ellinit([0,t*d**3]).ellrank()[:2])) for t in (-1,1)]
    bound=sum(b for _,b in ranks)
    w=pari('w')
    f=w**4-2*(d-1)*w**2+(d+1)**2
    nf=pari.nfinit(f)
    root=pari.Mod(w,f)
    i=(root**2-d-1)/(2*root)
    check(i*i==-1 and (root-i)**2==d,"Quartic field model")
    X=pari('x')
    points=set()
    def put(a,b):
        check(b*b==a**6-1,'Wrong fibre')
        points.add((str(pari.lift(a)),str(pari.lift(b))))
    candidates,n=gaussian_candidates()
    for (A,B,C),t in candidates:
        if t not in (0,1,d):
            continue
        x=(A+B*i)/C
        ys=pari.nfroots(nf,X**2-(x**6-1))
        check(len(ys)>0,'Gaussian squareclass lift')
        for y in ys:
            put(x,y)
    pairs=json.loads(Path(__file__).with_name('quadratic-pairs-over-gaussian.json').read_text())
    def gauss(pair):
        return pari(pair[0])+pari(pair[1])*i
    for row in pairs:
        if row['extension']!=d:
            continue
        U=pari.Polrev([gauss(c) for c in row['u']],X)
        V=pari.Polrev([gauss(c) for c in row['v']],X) if row['v'] else 0
        zs=pari.nfroots(nf,U)
        check(len(zs)==2,'Non-Gaussian divisor failed to split')
        for z in zs:
            put(1+1/z,pari.subst(V,X,z)/z**3)
    return dict(d=d,field=str(f),twist_rank_intervals=ranks,rank_upper_bound=bound,
                gaussian_abscissas_tested=n,verified_affine_points=sorted(points),
                rational_infinities=2,status='complete' if bound<=1 else 'partial')


def panel():
    return [classify_biquadratic(d) for d in range(2,201) if pari.issquarefree(d)]


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('d',type=int,nargs='?')
    cli.add_argument('--panel',action='store_true')
    cli.add_argument('--write',action='store_true')
    args=cli.parse_args()
    verify_identities()
    if args.panel:
        verify_non_gaussian_table()
        rows=json.loads(json.dumps(panel()))
        dest=Path(__file__).with_name('rank-one-biquadratic-receipts.json')
        if args.write:
            dest.write_text(json.dumps(rows,indent=2)+'\n')
        else:
            check(json.loads(dest.read_text())==rows,'Panel replay')
        print(Counter((r['rank_upper_bound'],r['status']) for r in rows))
        print('Point counts on complete cases:',Counter(2+len(r['verified_affine_points']) for r in rows if r['status']=='complete'))
    elif args.d is not None:
        print(json.dumps(classify_biquadratic(args.d),indent=2))
    else:
        cli.error('Supply a radicand or --panel.')
