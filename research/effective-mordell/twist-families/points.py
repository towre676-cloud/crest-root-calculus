"""Query the proved twist family and the two complete panels for positive d."""
import argparse
import json
from fractions import Fraction
from math import prod
from pathlib import Path

from cypari2 import Pari

HERE=Path(__file__).resolve().parent


def uniform(primes):
    q=[p for p in primes if p%12!=5]
    return len(q)<=1 and all(p%24==7 for p in q)


def point_list(d):
    if d<=0:
        raise ValueError('d must be a positive integer')
    pari=Pari()
    fac=pari.factorint(d)
    factors=[(int(p),int(e)) for p,e in zip(fac[0],fac[1])]
    if prod(p**e for p,e in factors)!=d or not all(pari.isprime(p) for p,e in factors):
        raise ArithmeticError('input factorization failed its primality or product check')
    k=prod(p**(e//2) for p,e in factors)
    sf=d//(k*k)
    atoms=[p for p,e in factors if e%2 and p!=3]
    base=prod(atoms)
    answer={'d':str(d),'squarefree_part':str(sf),'ordinate_scale':str(k)}
    if uniform(atoms):
        method='uniform theorem'
        points=[(Fraction(0),Fraction(0))]
    else:
        if base>200:
            return dict(**answer,status='outside the proved subclasses in this release')
        low=json.loads((HERE/'all-points.json').read_text())
        high=json.loads((HERE/'extension-receipt.json').read_text())['parameters']
        rows={r['d']:r for r in low}
        rows.update((int(i),r) for i,r in high.items())
        entry=rows[sf if sf<=200 else base]
        points=[(Fraction(x),Fraction(y)) for x,y in entry['affine']]
        method='certified parameter table'
        if sf>200:
            method+=' and the isomorphism to its triple twist'
            # The isomorphism exchanges the source origin and infinity. Both
            # are present, so keep the unique affine origin in the returned set.
            points=[(3/x,9*y/x**3) if x else (Fraction(0),Fraction(0)) for x,y in points]
    for x,y in points:
        if (y*k)**2!=d*(6*x**5+28*x**3+54*x):
            raise ArithmeticError('point does not lie on the input curve')
    return dict(**answer,status='complete',method=method,
                affine=[[str(x),str(y*k)] for x,y in sorted(points)],rational_infinities=1)


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('d',type=int,nargs='+')
    for d in parser.parse_args().d:
        if d<=0:
            parser.error('each d must be positive')
        print(json.dumps(point_list(d)))
