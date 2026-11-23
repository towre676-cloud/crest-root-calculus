"""Equation-only cubic-quotient decision, using one good residue characteristic."""
from pathlib import Path
import argparse
from fractions import Fraction
from math import gcd,isqrt
import json
import sys
import sympy as sp

HERE=Path(__file__).resolve().parent
for directory in ('scan','adic','torsion'):
    sys.path.insert(0,str(HERE.parent/directory))
from positive_search import fa,T,exhaustive_ticks,scan_flag,OLD
from adic_lifting import forest
from field_scan import scan
from frobenius_fields import certificate
from frobenius_obstruction import analyze


def rational_candidate(a,modulus):
    bound=isqrt((modulus-1)//2)
    r0,r=modulus,int(a)%modulus
    d0,d=0,1
    while r>bound:
        q=r0//r
        r0,r=r,r0-q*r
        d0,d=d,d0-q*d
    if not d or abs(d)>bound or gcd(d,modulus)!=1:return None
    cand=Fraction(r,d)
    if (cand.numerator-a*cand.denominator)%modulus:return None
    return sp.Rational(cand.numerator,cand.denominator)


def proposal(flag,modulus_poly):
    modulus=flag['modulus']
    q=[tuple(rational_candidate(c,modulus) for c in flag[s]) for s in ('u','v','s')]
    if any(c is None for row in q for c in row):return
    f=sum(c*T**i for i,c in enumerate(modulus_poly))
    if all(all(not v for v in row[1:]) for row in q):
        f=T
    u,v,s=(sum(c*T**i for i,c in enumerate(row)) for row in q)
    if flag['kind']=='IIwide':
        s=sp.invert(s,f,T) if s!=0 else None
    return tuple(flag['frame']),f,u,v,s


def choose_prime(F,at,scope):
    _,integral=F.clear_denoms(convert=True)
    _,primitive=integral.primitive()
    data={'terms':[[*m,int(c)] for m,c in primitive.terms()]}
    candidates=[]
    prime=at or 5
    while len(candidates)<(1 if at else 4):
        loc=analyze(data['terms'],prime)
        if loc['smooth']:
            bound=certificate(prime,loc['frobenius_coefficients'])
            degree=bound['one_sufficient_scan_degree'] if scope=='geometric' else 2
            if not bound['elliptic_dimension']:degree=0
            candidates.append((prime**degree,prime,degree,bound,loc))
            if degree==0:break
        elif at:
            raise ValueError('The specified prime gives singular reduction of this model.')
        prime=int(sp.nextprime(prime))
    _,p,d,b,loc=min(candidates,key=lambda v:(v[0],v[1]))
    return data,p,d,b,loc


def decide(original,at=None,scope='geometric',max_precision=None):
    if scope not in ('geometric','rational'):raise ValueError('Unknown scope.')
    if at is not None and (not sp.isprime(at) or at<=3):
        raise ValueError('Choose a prime at least five.')
    if max_precision is not None and max_precision<1:raise ValueError('Precision must be positive.')
    if scope=='rational' and max_precision is None:
        raise ValueError('The short rational-class run requires a finite precision budget. For guaranteed termination use geometric scope.')
    F=fa.input_quartic(original)
    fa.smooth_check(F)
    data,p,d,bound,loc=choose_prime(F,at,scope)
    F=fa.input_quartic(data)
    status={'scope':scope,'prime':p,'scan_degree':d,'field_bound':bound,
            'point_counts':loc['point_counts'],'levels':[]}
    orbits={}
    def finish(lower,upper):
        if lower>upper:raise ArithmeticError('A verified lower bound exceeded the local upper bound.')
        status.update(complete=lower==upper,lower=lower,upper=upper,verified_orbits=list(orbits.values()))
        if lower==upper:
            status['rational_classes']=sum(v['orbit_size']==1 for v in orbits.values())
        return status
    if d==0:return finish(0,0)
    report=scan(data,p,d)
    status['scan']=report
    ids=[k for k,j in enumerate(report['p_frobenius_permutation']) if scope=='geometric' or j==k]
    search=exhaustive_ticks(F)
    if not ids:return finish(0,0)
    known_flags=set()
    def add(items):
        for row in items:
            if scope=='geometric' or row['orbit_size']==1:
                orbits.setdefault(row['chow_polynomial'],row)
    for level in forest(data,report,ramification='adaptive',indices=ids,include_seeds=True):
        upper=len(level['survivors'])
        for flag in level['flags']:
            candidate=proposal(flag,report['modulus_constant_first'])
            if candidate is not None and candidate not in known_flags:
                known_flags.add(candidate)
                add(scan_flag(F,*candidate))
        # The independent height search gets finite time even if the local
        # coordinates do not reconstruct in this particular number field.
        for _ in range(16):add(next(search))
        low=sum(row['orbit_size'] for row in orbits.values())
        status['levels'].append({'precision':level['precision'],'branches':sum(level['branches'].values()),
                                 'lower':low,'upper':upper})
        if low>=upper or level['precision']==max_precision:return finish(low,upper)


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('input',type=Path)
    ap.add_argument('--prime',type=int)
    ap.add_argument('--scope',choices=('geometric','rational'),default='geometric')
    ap.add_argument('--max-precision',type=int)
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    result=decide(json.loads(args.input.read_text()),args.prime,args.scope,args.max_precision)
    payload=json.dumps(result,indent=2)+'\n'
    if args.output: args.output.write_text(payload)
    else: print(payload,end='')
