"""Exact fibre checks for Stoll's positive squarefree twists through 100.

Requires Python 3.11+, SymPy and cypari2. Rebuilds the finite computations
reported in panel.json and prints the extra squareclass descent results.
"""
from fractions import Fraction
from functools import lru_cache
from itertools import product
from math import isqrt
import json
from pathlib import Path
import sympy as sp
from cypari2 import Pari

pari = Pari()
pari.allocatemem(128_000_000)


def check(value, message):
    if not value:
        raise ArithmeticError(message)


def sqrtq(value):
    f = Fraction(value)
    if f < 0:
        return None
    n,d = isqrt(f.numerator),isqrt(f.denominator)
    return Fraction(n,d) if n*n == f.numerator and d*d == f.denominator else None


def source(x,d):
    return d*(6*x**5+28*x**3+54*x)


def direct_rank(d):
    curves = [pari.ellinit([a*d*d,0]) for a in (-108,-12)]
    return dict(d=d,ranks=[[int(v) for v in E.ellrank()[:2]] for E in curves],
                torsion=[int(E.elltors()[0]) for E in curves])


@lru_cache(None)
def square_residues(m):
    return {i*i % m for i in range(m)}


def normalized_pairs(p,m):
    yield from ((u,1) for u in range(m))
    yield from ((1,v) for v in range(0,m,p))


def blocked(d,a,p,m):
    squares = square_residues(m)
    for u,v in normalized_pairs(p,m):
        rhs = a*d*(6*a**4*u**8+28*a*a*u**4*v**4+54*v**8)
        if rhs % m in squares:
            return False
    return True


def finite_torsion(E):
    size,orders,generators = E.elltors()
    points = [pari([0])]
    for order,G in zip(orders,generators):
        points = [E.elladd(P,E.ellmul(G,i)) for P in points for i in range(int(order))]
    check(len(set(map(str,points)))==int(size),"incomplete torsion subgroup")
    return points


def examine(d,a):
    for p,m in ((2,256),(3,243),(5,125),(7,49),(11,121),(13,169),
                (17,17),(19,19),(23,23),(29,29),(31,31),(37,37),(41,41),
                (43,43),(47,47),(53,53),(59,59),(61,61),(67,67),(71,71),
                (73,73),(79,79),(83,83),(89,89),(97,97)):
        if blocked(d,a,p,m):
            return dict(a=a,obstruction=m)
    A,B,C = 6*d*a,28*d*a,54*d*a
    E=pari.ellinit([0,B,0,A*C,0])
    ranks=[int(i) for i in E.ellrank()[:2]]
    if ranks[1]>0:
        return dict(a=a,ranks=ranks)
    found=set()
    for P in finite_torsion(E):
        if len(P)==1:
            continue
        xx=sqrtq(Fraction(str(P[0]))/A)
        if xx is None or xx<=0 or sqrtq(xx/a) is None:
            continue
        if sqrtq(source(xx,d)) is not None:
            found.add(xx)
    return dict(a=a,ranks=ranks,rational_x=sorted(map(str,found)))


def lifts(d):
    primes=list(sp.factorint(6*d))
    return [examine(d,int(sp.prod(q for q,k in zip(primes,bits) if k)))
            for bits in product((0,1),repeat=len(primes))]


def delta_sieve(x_class):
    candidates=[(sg*(3-2*h*j),sg*(3*h+j))
                for sg,h,j in product((1,-1),repeat=3)]
    records=[]
    for a,b in candidates:
        # Test delta*(3*u^4+(7+4*s)*v^4) for being a square in O_K.
        # For split primes, test both images in Z/(p^3). For 2, use O_K/(16).
        obstacle=None
        for p,m in ((3,27),(11,1331),(2,16)):
            if p==2:
                squares={((i*i-2*j*j)%m,2*i*j%m) for i,j in product(range(m),repeat=2)}
                def has_square(u,v):
                    t=3*x_class**2*u**4+7*v**4; q=4*v**4
                    return ((a*t-2*b*q)%m,(a*q+b*t)%m) in squares
            else:
                s=int(sp.sqrt_mod(-2,m))
                def has_square(u,v):
                    return all((a+b*r)*(3*x_class**2*u**4+(7+4*r)*v**4)%m in square_residues(m)
                               for r in (s,-s))
            if not any(has_square(u,v) for u,v in normalized_pairs(p,m)):
                obstacle=m
                break
        records.append(dict(delta=[a,b],obstruction=obstacle))
    expected={1:[[-5,-2]],10:[]}
    check([r['delta'] for r in records if r['obstruction'] is None]==expected[x_class],
          "quadratic descent did not isolate its class")
    return records


if __name__=="__main__":
    path=Path(__file__).resolve().parent
    records=[direct_rank(d) for d in range(31,101)
             if all(d%(i*i) for i in range(2,isqrt(d)+1))]
    check(records==json.loads((path/'panel.json').read_text()),"rank or torsion mismatch")
    residual=[14,22]+[r['d'] for r in records if all(pair[1] for pair in r['ranks'])]
    further={}
    for d in residual:
        further[str(d)]=lifts(d)
        print(d,"surviving positive-rank squareclasses",
              sorted(v['a'] for v in further[str(d)] if v.get('ranks',[0,0])[1]>0))
    evidence={'squareclasses':further,'quadratic_sieve':
              {str(d):delta_sieve(a) for d,a in ((22,1),(55,10))}}
    output=path/'descent-receipt.json'
    if output.exists():
        check(json.loads(output.read_text())==evidence,"descent receipt mismatch")
    else:
        output.write_text(json.dumps(evidence,indent=2)+'\n')
