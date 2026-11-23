"""Replay the three all-degree examples, global correspondences, and independent controls."""
import argparse
import json
from itertools import combinations
from pathlib import Path

import sympy as sp

from all_primes import integers, rational_isogeny_class
from check_orbits import check_all
from discover import recover
from family import curves, packet
from frames import identify
from frobenius import X, compare, frobenius, polynomial
from kernel_fields import fields
from lattice_check import count_lattices
from rational_component import component
from rational_spectrum import exclude, irreducible_quartic, prove_irreducible, quadratic_tests
from richelot import correspondences, transform
from spectrum import counts

HERE=Path(__file__).resolve().parent


def require(condition,message):
    if not condition:
        raise ArithmeticError(message)


def independent_counts(f,p):
    """Count fibres using a table of squares in F_p[u], independently of the norm character."""
    c=[int(f.nth(j))%p for j in range(7)]
    squares=[0]*p
    for y in range(p):
        squares[y*y%p]+=1
    n1=squares[c[6]]+sum(squares[int(f.eval(x))%p] for x in range(p))
    nonsquare=next(i for i in range(1,p) if not squares[i])
    affine=[(a,b) for a in range(p) for b in range(p)]
    fibre={(a,b):0 for a,b in affine}
    def times(v,w):
        a,b=v;r,s=w
        return (a*r+nonsquare*b*s)%p,(a*s+b*r)%p
    for y in affine:
        fibre[times(y,y)]+=1
    n2=fibre[c[6],0]
    for x in affine:
        value=(0,0)
        for v in reversed(c):
            r,s=times(x,value);value=((r+v)%p,s)
        n2+=fibre[value]
    return n1,n2


def exterior_check(h,p):
    M=sp.zeros(4)
    for i in range(1,4):
        M[i,i-1]=1
    for i in range(4):
        M[i,3]=-h[4-i]
    A=M**120
    idx=list(combinations(range(4),2))
    wedge=sp.Matrix(6,6,lambda i,j:A.extract(idx[i],idx[j]).det())
    a,b,c=A.charpoly().all_coeffs()[1:4]
    expected=(p*(sp.eye(4)-A).det(),p*(sp.eye(6)-wedge).det(),p*(a*a-4*(b-2*p**120)))
    require(integers(h,p)==expected,'Exterior-power integer obstruction differs')


def make():
    summaries={}
    for stem,ps in [('example',[11,13,17,29]),('public',[7,17,19,29])]:
        inputs=json.loads((HERE/(stem+'.json')).read_text())
        f=polynomial(inputs['source'])
        g=polynomial(inputs['target'])
        r=recover(inputs['source'],inputs['target'],ps)
        require(r['geometric_hom'].startswith('Z generated') and r['rational_polarized_22_isogenies_up_to_sign']==1,
                'The full geometric Hom lattice was not found')
        arith=rational_isogeny_class(f)
        require(arith['classification'].startswith('exact rational degrees'),'All-prime rational classification failed')
        orbits=fields(f,True)
        require(check_all(f,orbits)==[1,6,8],'Incorrect geometric quotient field degrees')
        for p in ps[:2]:
            for curve in (f,g):
                require(independent_counts(curve,p)==frobenius(curve,p)[0],'Extension-field counts differ')
        summaries[stem]={'correspondence':r,'rational_spectrum':arith,'geometric_quotient_atlas':orbits}
    src=json.loads((HERE/'example.json').read_text())
    twist=recover(src['source'],[2*a for a in src['target']],[11,13,17,29])
    require(twist['rational_hom']=='zero' and twist['geometric_hom'].startswith('Z generated')
            and twist['graphs'][0]['field_squareclass']==2,'Twist descent failed')
    neg=src['source'].copy();neg[2]+=11
    collision=compare(polynomial(src['source']),polynomial(neg),[11,13])
    require(collision['geometric_hom_zero_prime']==13 and
            collision['reductions'][0]['source_weil_descending']==collision['reductions'][0]['target_weil_descending'],
            'Coincident local spectrum failed to separate globally')
    other=polynomial([0,-616,1222,-799,217,-25,1])
    dest=polynomial([14611320,66154320,92955186,20968290,-46101852,-29378016,-3466080])
    simple=correspondences(other,dest,[5,7,11,13])
    alternate=recover([int(other.nth(j)) for j in range(7)],[int(dest.nth(j)) for j in range(7)],[23,37])
    require(simple['source_pointwise_rational_kernels']==15 and len(simple['graphs_up_to_sign'])==1
            and alternate['graphs'][0]['target_mobius_matrix']==simple['graphs_up_to_sign'][0]['target_mobius_matrix']
            and alternate['geometric_hom'].startswith('Z generated'),'Split-branch independent construction failed')
    f=polynomial([-1,0,0,0,0,0,1])
    g=sp.Poly(3*transform(f,(3,2,1,1),6),X,domain=sp.QQ)
    mats,method=identify(f,g,[5,7])
    require(len(mats)==4 and method.startswith('height-certified'),'Exceptional binary sextic frame failed')
    q=sp.Poly(X**5-X+1,X,domain=sp.QQ)
    require(check_all(q,fields(q,True))==[15],'Quintic branch at infinity was not preserved')
    F,G=curves(1)
    for p in (17,29):
        h=frobenius(F,p)[1];exterior_check(h,p)
        for ell in sp.primerange(3,300):
            if ell!=p:
                require(prove_irreducible(quadratic_tests(h,p),ell)==irreducible_quartic(h,ell),
                        'Quadratic discriminants differ from polynomial factorization')
    lattice=[]
    for n in (1,2,3,4,6,8,9):
        c=counts(n)
        require(count_lattices(n)==(c['primitive_targets'],c['all_targets_at_this_degree']),
                'Integral Hermite-lattice enumeration disagrees with the formula')
        lattice.append(c)
    fam=packet()
    cert=summaries['example']['rational_spectrum']['all_odd_primes']
    used={3,17,29,59}|{p['p'] for p in cert['auxiliary_primes']}|set(cert['irreducible_witnesses'].values())
    family_mod=int(sp.prod(used))
    require(family_mod==30557478886935,'Changed congruence modulus in the all-prime family')
    star_source=polynomial(json.loads((HERE/'four-targets.json').read_text())['source'])
    star=component(star_source)
    require(star['complete'] and star['rational_richelot_edges']==[[1,2,3],[0],[0],[0]]
            and star['primitive_rational_isogeny_degrees']==
            [[1,4,4,4],[4,1,16,16],[4,16,1,16],[4,16,16,1]],'Incorrect rational degree-16 component')
    vertices=star['vertices_ascending']
    star['edge_graphs_from_center']=[]
    for leaf in vertices[1:]:
        edge=recover(vertices[0],leaf,[13,29])
        require(edge['rational_polarized_22_isogenies_up_to_sign']==1 and
                edge['rational_hom'].startswith('Z generated'),'Rational edge correspondence failed')
        star['edge_graphs_from_center'].append(edge['graphs'][0])
    return {'benchmarks':summaries,'quadratic_twist_field':twist['graphs'][0]['field_squareclass'],
            'local_collision':collision,'hermite_lattice_checks':lattice,
            'four_target_rational_component':star,
            'family':{**fam,'all_prime_family_modulus':family_mod}}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--extended',action='store_true',help='independently check all 664,578 odd primes below ten million')
    args=parser.parse_args()
    actual=json.loads(json.dumps(make()))
    expected=json.loads((HERE/'receipt.json').read_text())
    require(actual==expected,'Published genus-two research receipt differs from fresh computation')
    if args.extended:
        f=curves(1)[0]
        control=exclude(f,10_000_000,sp.primerange(3,702))
        require(control['all_powers_excluded']==664578 and control['assignments_sha256']==
                '161f5c932d79ea46b777cb7d2bc4504f977ffcb8ca2e6db2189536ee4626b8e4','Large finite cross-check differs')
    print('Verified three rational spectra, four-target degree-16 maps, 30 geometric model instances, and the infinite family.')
