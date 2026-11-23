"""Finite all-prime Frobenius obstruction for rational primitive polarized isogeny kernels."""
from math import gcd

import sympy as sp

from frobenius import X, frobenius, good
from rational_spectrum import irreducible_quartic


def lift(h,p,n=120):
    sums=[4]
    for k in range(1,2*n+1):
        value=sum(h[i]*sums[k-i] for i in range(1,min(4,k-1)+1))
        if k<=4:
            value+=k*h[k]
        sums.append(-value)
    q=p**n
    a=-sums[n]; b=(sums[n]**2-sums[2*n])//2
    return a,b,q


def integers(h,p):
    a,b,q=lift(h,p)
    lines=1+a+b+q*a+q*q
    # The six eigenvalues on wedge^2 contain q twice and four cross-products.
    t=2*q-b
    cross=1+t+q*(a*a-2*b+2*q)+q*q*t+q**4
    det_extreme=(1-q)**2*cross
    det_middle=a*a-4*(b-2*q)
    return p*lines,p*det_extreme,p*det_middle


def analyze(f,prime_bound=200):
    # Any prime not dividing this integer is good for the displayed hyperelliptic model.
    bad=abs(int(2*f.LC()*f.discriminant()))
    divisors={2,3,*map(int,sp.factorint(bad))}
    reduced=[0,0,0]
    records=[]
    primes=[int(p) for p in sp.primerange(3,prime_bound+1) if good(f,int(p))]
    for p in primes:
        counts,h=frobenius(f,p)
        obs=integers(h,p)
        next_gcd=[gcd(a,abs(b)) for a,b in zip(reduced,obs)]
        if next_gcd!=reduced:
            records.append({'p':p,'counts':counts,'weil_descending':h})
            reduced=next_gcd
        # Wait until all three nonzero gcds are small enough for direct factorization.
        if all(r and r.bit_length()<110 for r in reduced):
            break
    if any(r==0 for r in reduced):
        raise ValueError('Additional ordinary auxiliary primes are required')
    for g in reduced:
        divisors.update(map(int,sp.factorint(g)))
    todo=sorted(v for v in divisors if v!=2)
    excluded={};survivors=[]
    for ell in todo:
        for p in primes:
            if p==ell:
                continue
            counts,h=frobenius(f,p)
            if irreducible_quartic(h,ell):
                excluded[ell]=p
                break
        else:
            survivors.append(ell)
    return {'powered_frobenius_exponent':120,'auxiliary_primes':records,'three_gcds':reduced,
            'model_bad_prime_candidates':sorted(map(int,sp.factorint(bad))),
            'residual_odd_candidate_primes':todo,'irreducible_witnesses':excluded,
            'unresolved_odd_primes':survivors}

def rational_isogeny_class(source,prime_bound=200):
    """Sufficient certificate for the two-target, all-degree rational spectrum."""
    from frobenius import compare
    from partitions import stable_pairings,dual_algebra,Z
    working,_,parts=stable_pairings(source)
    factors=working.factor_list()[1]
    rat_two_torsion=any(q.degree()==2 for q,e in factors) or sum(q.degree()==1 for q,e in factors)>=2
    cert={'source_stable_pairings':len(parts),'source_has_rational_nonzero_2_torsion':rat_two_torsion,
          'all_odd_primes':analyze(source,prime_bound),'classification':'undetermined'}
    if len(parts)!=1 or rat_two_torsion or cert['all_odd_primes']['unresolved_odd_primes']:
        return cert
    dual=dual_algebra(working,parts[0])
    if dual is None:
        return cert
    h=sp.Poly(dual['dual'].as_expr().subs(Z,X),X,domain=sp.QQ)
    den,_=h.clear_denoms(convert=True)
    h=h.mul_ground(den*den)  # retain the rational twist
    contents=h.primitive()[0]
    square=sp.prod(p**(n//2) for p,n in sp.factorint(abs(int(contents))).items())
    h=h.mul_ground(sp.Rational(1,square**2))
    _,_,next_parts=stable_pairings(h)
    end=compare(source,h,list(sp.primerange(3,75)))['source_endomorphisms_Z_witnesses']
    cert.update({'endomorphism_witnesses':end,'dual_stable_pairings':len(next_parts),
                 'dual_sextic_ascending':[int(h.nth(i)) for i in range(7)]})
    if len(next_parts)==1 and end:
        cert['classification']='exact rational degrees n^4 and 4*n^4, one class per degree'
    return cert


if __name__=='__main__':
    import argparse
    import json
    from pathlib import Path
    from frobenius import polynomial
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--witness-limit',type=int,default=200)
    args=parser.parse_args()
    source=polynomial(json.loads(args.input.read_text())['source'])
    print(json.dumps(rational_isogeny_class(source,args.witness_limit),indent=2))
