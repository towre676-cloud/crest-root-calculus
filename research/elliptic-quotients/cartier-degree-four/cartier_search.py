"""Find rational elliptic map directions and certify geometric exhaustion when possible."""
import argparse
from itertools import permutations, product
import json
from math import gcd
from pathlib import Path
import sys

import numpy as np
import sympy as sp

from recover import ELLIPTIC, S, T, lift_target

sys.path.insert(0,str(ELLIPTIC / 'polarized-reduction'))
from screen import algebra, local_profile, read_terms
sys.path.insert(0,str(ELLIPTIC.parent / 'cubic-tower-factors'))
from field import Field

XYZ = sp.symbols('U V W')


def quartic(form, coords):
    return sum(c*coords[0]**i*coords[1]**j*coords[2]**k for (i,j,k),c in form.items())


def projective_points(height):
    # Each primitive rational point is represented once, with first nonzero coordinate positive.
    for h in range(1,height+1):
        shell = [v for v in product(range(-h,h+1), repeat=3)
                 if max(map(abs,v))==h and gcd(*v)==1 and next(x for x in v if x)>0]
        yield from sorted(shell,key=lambda v:(sum(x<0 for x in v),v))


def cartier(form, p):
    """Columns are images of the three canonical differentials modulo p."""
    edge = 2*p-1
    power = np.zeros((edge+1,edge+1), dtype=np.int64)
    power[0,0] = 1
    for _ in range(p-1):
        nxt = np.zeros_like(power)
        for (i,j,k),c in form.items():
            nxt[i:,j:] += (c%p)*power[:edge+1-i,:edge+1-j]
        power = nxt % p
    basis = [(2,1),(1,2),(1,1)]
    return np.array([[int(power[p*a-i,p*b-j]) for i,j in basis] for a,b in basis])


def possible_directions(form, height, primes):
    matrices = {p:cartier(form,p) for p in primes if algebra.smoothness(form,p)}
    candidates = np.asarray(list(projective_points(height)), dtype=np.int64)
    for p,matrix in matrices.items():
        good = np.all(np.cross(candidates,candidates@matrix.T)%p==0,axis=1)
        candidates = candidates[good]
    return [tuple(map(int,v)) for v in candidates], {str(p):m.tolist() for p,m in matrices.items()}


def axis_test(F, line):
    """First polar on the proposed line, then the full reflection identity."""
    h = sp.Matrix([line])
    base = h.nullspace()
    a,b,c = sp.symbols('a b c')
    coords = a*base[0]+b*base[1]
    polys = [sp.Poly(F.diff(z).subs(dict(zip(XYZ,coords)),simultaneous=True),a,b) for z in XYZ]
    constraints = sp.Matrix([[f.coeff_monomial(a**i*b**(3-i)) for f in polys] for i in range(4)])
    centers = constraints.nullspace()
    if len(centers)>1:
        raise ValueError('Singular source: polar restriction has rank below two')
    if not centers or h.dot(centers[0])==0:
        return None
    frame = sp.Matrix.hstack(*base,centers[0])
    transformed = F.subs(dict(zip(XYZ,frame*sp.Matrix([a,b,c]))),simultaneous=True).expand()
    if transformed != transformed.subs(c,-c):
        return None
    den, even = sp.Poly(transformed,a,b,c).clear_denoms(convert=True)
    even = even.primitive()[1]
    cover = sp.Poly(even.as_expr(),c)
    branch = sp.Poly(cover.nth(2)**2-4*cover.nth(4)*cover.nth(0),a,b)
    aa,bb,cc,dd,ee = [branch.coeff_monomial(a**i*b**(4-i)) for i in range(4,-1,-1)]
    I = 12*aa*ee-3*bb*dd+cc**2
    J = 72*aa*cc*ee+9*bb*cc*dd-27*aa*dd**2-27*bb**2*ee-2*cc**3
    target_j = sp.cancel(6912*I**3/(4*I**3-J**2))
    return {'center':list(map(str,centers[0])), 'j':str(target_j), 'frame':[[str(t) for t in row] for row in frame.tolist()],
            'even_model':[[*m,int(v)] for m,v in even.terms()]}


def suitable_chart(form, lines, height):
    F = quartic(form,XYZ)
    for point in projective_points(height):
        if quartic(form,point) or any(sum(a*b for a,b in zip(l,point))==0 for l in lines):
            continue
        for fixed in (2,1,0):
            if not point[fixed]:
                continue
            for i,j in permutations([k for k in range(3) if k != fixed]):
                slope = F.diff(XYZ[j]).subs(dict(zip(XYZ,point)))
                if slope:
                    values = list(point)
                    values[i] += T
                    values[j] += S
                    return (point, i, j, slope, quartic(form,values).expand())
    return None


def even_count(model, p, degree):
    """Sum the square-root counts in the fibres of the even quartic over P^1."""
    field = Field(p,degree)
    def coefficient(w_degree):
        val = np.full(field.q,-1,dtype=np.int64)
        infinity = -1
        for i,j,k,c in model:
            if k != w_degree or not c%p:
                continue
            mon = field.log.copy()*i % field.order if i else np.zeros(field.q,dtype=np.int64)
            if i:
                mon[0] = -1
            val = field.add(val,field.shift(mon,field.log[c%p]))
            if j==0:
                infinity = field.add(infinity,field.log[c%p])
        return np.append(val,infinity)
    a,b,c = [coefficient(d) for d in (4,2,0)]
    def mul(a,b):
        return np.where((a<0)|(b<0),-1,(a+b)%field.order)
    def negate(a):
        return field.shift(a,field.log[p-1])
    disc = field.add(mul(b,b),mul(field.shift(a,field.log[-4%p]),c))
    valid = (disc==-1)|(disc%2==0)
    sqrt_disc = np.where(disc==-1,-1,disc//2)
    inv = -field.shift(a,field.log[2]) % field.order
    def solutions(num):
        value = mul(num,inv)
        return np.where(value<0,1,np.where(value%2==0,2,0))
    plus = solutions(field.add(negate(b),sqrt_disc))
    minus = solutions(field.add(negate(b),negate(sqrt_disc)))
    return int(np.sum(valid*(plus+(disc!=-1)*minus)))


def geometric_bound(form, p, axes):
    if not sp.isprime(p) or not 5<=p<=97 or not algebra.smoothness(form,p):
        raise ValueError('Select a smooth good prime between 5 and 97')
    model = None
    for axis in axes:
        frame = sp.Matrix(axis['frame'])
        trial = {tuple(m[:3]):m[3] for m in axis['even_model']}
        den = sp.ilcm(*(v.q for v in frame))
        if den%p and int((den*frame).det())%p and trial.get((0,0,4),0)%p:
            model = axis['even_model']
            break
    if model is None and p>11:
        raise ValueError('For the direct projective counter use a good prime at most 11')
    counts = [even_count(model,p,k) if model else algebra.point_count(form,p,k) for k in (1,2,3)]
    weil = algebra.characteristic(p,counts)[::-1]
    result = local_profile(p,weil)
    blocks = result['elliptic_blocks']
    cap = None if any(b['geometric_multiplicity']>1 for b in blocks) else len(blocks)
    return {'point_counts':counts,'local_profile':result,'geometric_optimal_upper_bound':cap}


def discover(path, height=4, point_height=16, bound_prime=17, prime_limit=139):
    form = read_terms(path)
    F = quartic(form,XYZ)
    lines, matrices = possible_directions(form,height,[5,7,11,13,17,19,23])
    axes, possible = [], []
    for line in lines:
        axis = axis_test(F,line)
        if axis:
            axes.append({'line':line, **axis})
        else:
            possible.append(line)
    bound = geometric_bound(form,bound_prime,axes)
    if bound['geometric_optimal_upper_bound'] == 0:
        possible = []
    result = {'line_height':height,'point_height':point_height,'cartier_matrices':matrices,
              'double_quotients':axes, 'local_bound':bound,'primitive_degree_four':[]}
    chart = suitable_chart(form,possible,point_height) if possible else None
    if chart:
        point,i,j,slope,polynomial = chart
        result['chart'] = {'point':point,'t_coordinate':i,'s_coordinate':j}
        for line in possible:
            scale = slope/sum(a*b for a,b in zip(point,line))
            u,v = scale*line[i],scale*line[j]
            found = lift_target(polynomial,u,v,sp.primerange(29,prime_limit))
            if found:
                result['primitive_degree_four'].append({'line':line,'u':str(u),'v':str(v),**found})
    count = len(axes)+len(result['primitive_degree_four'])
    result['all_optimal_directions_exhausted'] = count == bound['geometric_optimal_upper_bound']
    return result


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('source',type=Path)
    parser.add_argument('--line-height',type=int,default=4)
    parser.add_argument('--point-height',type=int,default=16)
    parser.add_argument('--bound-prime',type=int,default=17)
    parser.add_argument('--prime-limit',type=int,default=139)
    args = parser.parse_args()
    print(json.dumps(discover(args.source,args.line_height,args.point_height,args.bound_prime,args.prime_limit),indent=2))
