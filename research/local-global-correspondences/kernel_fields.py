"""The exact fields of definition of all fifteen polarized genus-two Richelot kernels."""
import argparse
from itertools import combinations
import json
from pathlib import Path

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_pow_mod

from frobenius import X, polynomial
from partitions import centered, mul_list, pairings, sextic_chart

Y=sp.Symbol('u')


def separating_prime(f):
    p=13
    while True:
        p=int(sp.nextprime(p))
        h=[int(c)%p for c in f.all_coeffs()]
        if gf_pow_mod([1,0],p,h,p,ZZ)!=[1,0]:
            continue
        roots=[j for j in range(p) if int(f.eval(j))%p==0]
        if len(roots)!=6:
            raise ArithmeticError('Split prime did not have six distinct roots')
        for k in range(min(p,106)):
            if len({(a+b+k*a*b)%p for a,b in combinations(roots,2)})!=15:
                continue
            cubic_roots=[[a+b+k*a*b for a,b in m] for m in pairings(roots)]
            for b in range(min(p,211)):
                if len({int(sp.prod(b-t for t in h))%p for h in cubic_roots})==15:
                    return p, k, b, roots


def coefficients_for_pairing(match,k,b,leading,modulus):
    # Ascending powers for each branch quadratic and its derivatives' Wronskians.
    triples=[[(u*v)%modulus,(-u-v)%modulus,1] for u,v in match]
    c0,c1,c2=triples
    delta=(c1[1]*c2[0]-c1[0]*c2[1]-c0[1]*(c2[0]-c1[0])+c0[0]*(c2[1]-c1[1]))%modulus
    brackets=[]
    for i in range(3):
        c,d=triples[(i+1)%3],triples[(i+2)%3]
        brackets.append([(c[1]*d[0]-c[0]*d[1])%modulus,2*(d[0]-c[0])%modulus,(d[1]-c[1])%modulus])
    dual=[(-leading*delta*v)%modulus for v in mul_list(mul_list(brackets[0],brackets[1]),brackets[2])]
    values=[(u+v+k*u*v)%modulus for u,v in match]
    m=[1]
    for v in values:
        m=mul_list(m,[-v,1])
    numerator=[0,0,0]
    for i,(u,v) in enumerate(match):
        base=mul_list([-values[(i+1)%3],1],[-values[(i+2)%3],1])
        for j,w in enumerate(base):
            numerator[j]+=u*v*w
    label=int(sp.prod(b-t for t in values))%modulus
    return label,[v%modulus for v in m[:3]+numerator+dual+[delta*delta]]


def field_vector(poly):
    # Constant coefficient first; trailing zero coefficients are omitted.
    return [str(poly.nth(i)) for i in range(max(0,int(poly.degree())+1))] if poly else []


def fields(f,models=False):
    f,chart=sextic_chart(f)
    a=int(f.LC())
    integral=sp.Poly.from_list([f.nth(i)*sp.Integer(a)**(5-i) for i in range(6,-1,-1)],X,domain=sp.ZZ)
    p,k,b,roots=separating_prime(integral)
    R=1+max(abs(int(c)) for c in integral.all_coeffs()[1:])
    theta_bound=2*R+k*R**2
    lambda_bound=(abs(b)+theta_bound)**3
    height=(1+lambda_bound)**15
    if models:
        C=(1+R)**2
        coefficient_bound=max((1+theta_bound)**3,3*R**2*(1+theta_bound)**2,10368*abs(a)*C**9)
        height=max(height,15*coefficient_bound*(1+lambda_bound)**14)
    modulus=p
    deriv=integral.diff()
    while modulus<=2*height:
        modulus*=modulus
        roots=[(r-int(integral.eval(r))*pow(int(deriv.eval(r)),-1,modulus))%modulus for r in roots]
    records=[coefficients_for_pairing(match,k,b,a,modulus) for match in pairings(roots)]
    local=[label for label,rec in records]
    if len({v%p for v in local})!=15:
        raise ArithmeticError('Kernel labels failed to separate')
    coeff=[1]
    for v in local:
        coeff=[c%modulus for c in mul_list(coeff,[-v,1])]
    result=sp.Poly.from_list([centered(c,modulus) for c in coeff[::-1]],Y,domain=sp.QQ)
    factors=result.factor_list()[1]
    if any(m!=1 for g,m in factors) or sum(g.degree() for g,m in factors)!=15:
        raise ArithmeticError('The resolvent is not etale of degree fifteen')
    factors=sorted((g for g,m in factors),key=lambda g:(g.degree(),g.all_coeffs()))
    out={'source_chart':list(chart),'monic_branch_multiplier':a,'split_prime':p,'pair_parameter':k,
         'matching_evaluation':b,'coefficient_height_bound':height,'hensel_modulus':modulus,
         'orbit_degrees':[g.degree() for g in factors],
         'kernel_field_polynomials_descending':[list(map(int,g.all_coeffs())) for g in factors]}
    if not models:
        return out
    interpolation=[[0]*15 for _ in range(14)]
    for i,(label,values) in enumerate(records):
        basis=[1]
        for j,v in enumerate(local):
            if i!=j:
                basis=[c%modulus for c in mul_list(basis,[-v,1])]
        for dest,val in zip(interpolation,values):
            for n,q in enumerate(basis):
                dest[n]=(dest[n]+val*q)%modulus
    numerators=[sp.Poly.from_list([centered(c,modulus) for c in row[::-1]],Y,domain=sp.QQ)
                for row in interpolation]
    out['orbit_models']=[]
    for g in factors:
        inv=result.diff().invert(g)
        values=[(n*inv).rem(g) for n in numerators]
        delta2=values[-1]
        sextic=([(q*delta2.invert(g)).rem(g) for q in values[6:13]] if delta2 else [])
        out['orbit_models'].append({'pair_cubic_lower_coefficients':list(map(field_vector,values[:3])),
                                    'pair_product_numerator':list(map(field_vector,values[3:6])),
                                    'target_sextic_ascending':list(map(field_vector,sextic)),
                                    'product_quotient':not bool(delta2)})
    return out


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--models',action='store_true',help='give every quotient over its kernel field')
    args=parser.parse_args()
    src=json.loads(args.input.read_text())
    print(json.dumps(fields(polynomial(src['source']),args.models),indent=2))
