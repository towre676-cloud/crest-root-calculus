"""Exhaust a rational Richelot component and certify primitive isogeny multipliers by distance."""
import argparse
from collections import deque
import json
from pathlib import Path

import sympy as sp

from all_primes import analyze
from frames import identify
from frobenius import X, compare, polynomial
from partitions import Z, dual_algebra, sextic_chart, stable_pairings
from richelot import rational_sqrt


def integral_twist(h):
    den,_=h.clear_denoms(convert=True)
    h=h.mul_ground(den**2)
    content=abs(int(h.primitive()[0]))
    square=sp.prod(p**(e//2) for p,e in sp.factorint(content).items())
    return h.mul_ground(sp.Rational(1,square**2))


def has_rational_two_torsion(f):
    g,_=sextic_chart(f)
    factors=g.factor_list()[1]
    return any(q.degree()==2 for q,e in factors) or sum(q.degree()==1 for q,e in factors)>=2


def same_model(f,g):
    if f==g:
        return True
    frames,_=identify(f,g,[5,7,11,13])
    return any(rational_sqrt(s) is not None for _,s in frames)


def component(f, limit=None):
    end=compare(f,f,list(sp.primerange(3,79)))['source_endomorphisms_Z_witnesses']
    if not end:
        return {'complete':False,'reason':'source End=Z not certified'}
    vertices=[f]; edges=[];i=0
    while i<len(vertices):
        g,_,parts=stable_pairings(vertices[i])
        links=[]
        for part in parts:
            dual=dual_algebra(g,part)
            if dual is None:
                raise ArithmeticError('Product quotient on an absolutely simple Jacobian')
            h=integral_twist(sp.Poly(dual['dual'].as_expr().subs(Z,X),X,domain=sp.QQ))
            index=next((j for j,v in enumerate(vertices) if same_model(h,v)),None)
            if index is None:
                if limit is not None and len(vertices)>=limit:
                    return {'complete':False,'reason':'vertex budget reached','vertices_reached':len(vertices)}
                index=len(vertices);vertices.append(h)
            links.append(index)
        if len(set(links))!=len(links) or i in links:
            raise ArithmeticError('Degree-four kernels collided on a source with End=Z')
        edges.append(sorted(links));i+=1
    if any(i not in edges[j] for i,adj in enumerate(edges) for j in adj):
        raise ArithmeticError('A rational dual edge is missing')
    odd=analyze(f)
    spectra=[]
    for i,g in enumerate(vertices):
        if has_rational_two_torsion(g) or odd['unresolved_odd_primes']:
            spectra.append(None);continue
        distances={i:0};queue=deque([i])
        while queue:
            v=queue.popleft()
            for w in edges[v]:
                if w not in distances:
                    distances[w]=distances[v]+1;queue.append(w)
        spectra.append([4**distances[j] for j in range(len(vertices))])
    return {'complete':True,'source_endomorphism_witnesses':end,'odd_prime_certificate':odd,
            'vertices_ascending':[[int(g.nth(j)) for j in range(7)] for g in vertices],
            'rational_richelot_edges':edges,'primitive_rational_isogeny_degrees':spectra}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    parser.add_argument('--limit',type=int)
    args=parser.parse_args()
    source=polynomial(json.loads(args.input.read_text())['source'])
    print(json.dumps(component(source,args.limit),indent=2))
