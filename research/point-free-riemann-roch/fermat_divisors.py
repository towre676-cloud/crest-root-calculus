"""Parametrize rational effective divisors on X^4+Y^4+Z^4=0 from its rational Picard group."""
import argparse
import json

import sympy as sp

from riemann_roch import PlaneCurve, U, V, X, Y, Z, coefficients, kernel
from symmetric_power import chords


def bitangent_divisors(curve):
    found = list(chords(curve,1))
    if len(found) != 4 or any(quad.has(X) for _,_,quad in found):
        raise ArithmeticError('Unexpected rational bitangent conventions on the Fermat curve')
    return found


def parameterize(degree, include_bases=True):
    if type(degree) is not int or degree < 0:
        raise ValueError('Choose a nonnegative degree')
    if degree == 0:
        return [{'projective_dimension':0,'divisor':'empty divisor'}]
    if degree%2:
        return []
    n = degree//2
    curve = PlaneCurve(X**4+Y**4+Z**4)
    data = bitangent_divisors(curve)
    base = data[0][1]
    ideals = [curve.ideal([(1,[L.xreplace({X:U,Y:V,Z:sp.S.One}),
                                 q.xreplace({Y:V,Z:sp.S.One})])]) for _,L,q in data]
    answers = []
    for i,(_,line,_) in enumerate(data):
        if i == 0:
            power = (n+1)//2
            denominator = base**power
            which = [0] if n%2 else []
        else:
            k = n//2
            power = k+1
            denominator = base**k*line
            which = [0,i] if n%2 == 0 else [i]
        expected = 1 if n==1 else (3 if i==0 else 2) if n==2 else 2*n-2
        answer = {'class':i,'base_divisor':f'{n}D0' if i==0 else f'{n-1}D0+D{i}',
                  'projective_dimension':expected-1,'denominator':str(denominator)}
        if include_bases:
            mon = curve.monomials(power)
            cols = []
            for m in mon:
                col = {}
                for j in which:
                    rem = ideals[j].reduce(m)[1]
                    col.update({(j,*key):c for key,c in coefficients(rem).items()})
                cols.append(col)
            nums = [curve.homogenize(sum(v*m for v,m in zip(vec,mon)), power) for vec in kernel(cols)]
            if len(nums) != expected:
                raise ArithmeticError('The dimension differs from the Picard/Riemann--Roch formula')
            answer['numerators'] = list(map(str,nums))
        answers.append(answer)
    return answers


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('degree',type=int)
    parser.add_argument('--dimensions-only',action='store_true')
    args=parser.parse_args()
    print(json.dumps(parameterize(args.degree,not args.dimensions_only),indent=2))
