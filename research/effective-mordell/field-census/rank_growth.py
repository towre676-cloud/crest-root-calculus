"""Construct independent quadratic fields with high abscissas and certified rank contributions."""
import argparse
from math import isqrt
import json
from cypari2 import Pari

pari=Pari()

def witnesses(count):
    if count<1:
        raise ValueError('The number of independent radicals must be positive.')
    support={2,3}
    previous=3
    characters=set()
    rows=[]
    for _ in range(count):
        q=int(pari.nextprime(previous+1))
        while q in support:
            q=int(pari.nextprime(q+1))
        if not pari.isprime(q):
            raise ArithmeticError('Expected a proven prime')
        previous=q
        n=1+q
        N=n**6-1
        factors,exponents=pari.factorint(N)
        d=1
        for prime,val in zip(factors,exponents):
            if not pari.isprime(prime):
                raise ArithmeticError('The factorization has an unproved prime')
            if int(val)%2:
                d*=int(prime)
                support.add(int(prime))
        if N%(q*q)==0 or d%q!=0:
            raise ArithmeticError('A new ramified prime must have odd valuation')
        m=isqrt(N//d)
        if d*m*m!=N:
            raise ArithmeticError('Incorrect squarefree part')
        lifts=[(d,pari(d)*n*n,pari(d*d)*m),(-d,-pari(d)/(n*n),pari(d*d)*m/n**3)]
        for twist,X,Y in lifts:
            if Y*Y!=X**3-twist**3:
                raise ArithmeticError('Rational elliptic twist witness failed')
        if n**6<=132:  # 2 log(n) > log(132)/3, hence both witnesses have positive height
            raise ArithmeticError('Height fails to certify infinite order')
        four={int(pari.core(s*d)) for s in (1,-1,3,-3)}
        if characters&four or len(four)!=4:
            raise ArithmeticError('Character directions overlap')
        characters.update(four)
        rows.append(dict(prime=q,abscissa=n,radicand=d,square_factor=m,
                         characters=sorted(four),twist_witnesses=
                         [[t,str(X),str(Y)] for t,X,Y in lifts]))
    return dict(independent_radicals=count,compositum_degree=4*2**count,
                elliptic_rank_lower_bound=len(characters),witnesses=rows)

if __name__=='__main__':
    app=argparse.ArgumentParser(description=__doc__)
    app.add_argument('count',nargs='?',type=int,default=6)
    print(json.dumps(witnesses(app.parse_args().count),indent=2))
