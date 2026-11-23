"""Point-free Riemann--Roch and symmetric-power charts for embedded smooth integral curves over Q."""
import argparse
from collections import defaultdict
from itertools import combinations, count, product
import json
from math import comb, gcd
from pathlib import Path

import sympy as sp

from riemann_roch import kernel, polynomial, unit


def gb(equations, variables):
    return sp.groebner(equations, *variables, order='grevlex', domain=sp.QQ)


def hilbert(ideal):
    """Return the numerator Q and pole dimension b for H(s)=Q(s)/(1-s)^b."""
    count = len(ideal.gens)
    shifts = {(0,)*count: 1}
    for p in ideal.polys:
        leading = p.LM(order=ideal.order).exponents
        changes = defaultdict(int)
        for powers, weight in shifts.items():
            multiple = tuple(max(i,j) for i,j in zip(powers, leading))
            changes[multiple] -= weight
        for multiple, weight in changes.items():
            shifts[multiple] = shifts.get(multiple,0)+weight
            if not shifts[multiple]:
                del shifts[multiple]
    s = sp.Symbol('s')
    q = sp.Poly(sum(w*s**sum(e) for e,w in shifts.items()), s, domain=sp.QQ)
    while count and not q.eval(1):
        q = q.exquo(sp.Poly(1-s, s)); count -= 1
    return q, count


def degree_exponents(n, slots):
    if slots == 1:
        yield (n,)
    else:
        for i in range(n+1):
            for sub in degree_exponents(n-i,slots-1):
                yield (i,*sub)


def independent(columns):
    keys = sorted(set().union(*columns)) if columns else []
    matrix = sp.zeros(len(keys),len(columns))
    positions = {key:i for i,key in enumerate(keys)}
    for j,col in enumerate(columns):
        for key,c in col.items():
            matrix[positions[key],j] = c
    return matrix.rref()[1]


class EmbeddedCurve:
    """Validate a smooth geometrically integral curve and work over its field of definition."""
    def __init__(self, variables, equations):
        self.variables = tuple(variables)
        self.equations = tuple(equations)
        if len(variables) < 2 or not all(sp.Poly(f,*variables,domain=sp.QQ).is_homogeneous for f in equations if f):
            raise ValueError('Expected homogeneous rational equations in at least two coordinates')
        self.ideal = gb(equations,variables)
        q, pole = hilbert(self.ideal)
        if pole != 2:
            raise ValueError('The homogeneous ideal must present a projective curve')
        self.degree = int(q.eval(1))
        self.genus = int(1-self.degree+q.diff().eval(1))
        if self.genus < 0:
            raise ValueError('The input has negative arithmetic genus')
        self.leading = [p.LM(order=self.ideal.order).exponents for p in self.ideal.polys]
        self.charts = []
        for idx, chart_var in enumerate(variables):
            v = tuple(x for x in variables if x != chart_var)
            loc = gb([f.subs(chart_var,1) for f in equations],v)
            if unit(loc):
                continue
            equations_here = list(loc)
            jac = sp.Matrix([[f.diff(x) for x in v] for f in equations_here])
            codim = len(v)-1
            minors = [jac.extract(rows,cols).det() for rows in combinations(range(len(equations_here)),codim)
                      for cols in combinations(range(len(v)),codim)] if codim else [1]
            if not unit(gb(equations_here+minors,v)):
                raise ValueError('Jacobian rank drops below the smooth curve codimension')
            self.charts.append((idx,chart_var,v,loc))
        if not self.charts:
            raise ValueError('Empty projective source')
        # A one-dimensional projective complete intersection is connected over an algebraic closure.
        if sum(bool(f) for f in equations) != len(variables)-2:
            self.check_connectedness()

    def check_connectedness(self):
        # Taylor resolution + Groebner degeneration bound the regularity of the coordinate ring.
        bound = sum(max(m[i] for m in self.leading) for i in range(len(self.variables))) if self.leading else 0
        for parameter in count():
            line = sum(parameter**i*x for i,x in enumerate(self.variables))
            q,pole = hilbert(gb([*self.equations,line],self.variables))
            if pole == 1 and q.eval(1) == self.degree:
                break
        # A length-degree hyperplane slice certifies injectivity of multiplication by line
        # as a sheaf, by the Hilbert polynomial exact sequence.
        mon = self.ambient_sections(bound)
        cols = [{} for _ in mon]
        for idx,z,variables,loc in self.charts:
            denominator = line.subs(z,1)**bound
            target = gb([*loc,denominator],variables)
            for j,f in enumerate(mon):
                rem = sp.Poly(target.reduce(f.subs(z,1))[1],*variables,domain=sp.QQ)
                cols[j].update({(idx,*a):c for a,c in rem.as_dict().items()})
        if len(kernel(cols)) != 1:
            raise ValueError('The projective source lacks geometric connectedness or has isolated scheme components')

    @classmethod
    def from_data(cls,data):
        variables = sp.symbols('x0:'+str(data['coordinates']))
        return cls(variables,[polynomial(f,variables) for f in data['equations']])

    def pieces(self,data):
        return [(c.get('multiplicity',1),[polynomial(g,self.variables) for g in c['generators']]) for c in data]

    def divisor_data(self,pieces):
        h = sp.S.One
        degree = 0
        for mult, generators in pieces:
            if type(mult) is not int or mult < 1 or not generators:
                raise ValueError('Expected a finite divisor and positive multiplicity')
            generators = [sp.sympify(g) for g in generators]
            if any(not sp.Poly(g,*self.variables,domain=sp.QQ).is_homogeneous for g in generators if g):
                raise ValueError('Projective ideal generators must be homogeneous')
            q,pole = hilbert(gb([*self.equations,*generators],self.variables))
            if pole > 1:
                raise ValueError('Projective divisor support is positive-dimensional')
            if not pole:
                continue
            degree += mult*int(q.eval(1))
            # A generator nonzero on a dense affine open cuts a finite divisor on an integral source.
            nz = [g for g in generators if any(loc.reduce(g.subs(z,1))[1] for _,z,_,loc in self.charts)]
            h *= min(nz,key=lambda q: sp.Poly(q,*self.variables).total_degree())**mult
        return h,degree

    def ambient_sections(self,n):
        mon = [sp.prod(x**j for x,j in zip(self.variables,e)) for e in degree_exponents(n,len(self.variables))
               if not any(all(i<=j for i,j in zip(m,e)) for m in self.leading)]
        cols = []
        for f in mon:
            col = {}
            for idx,z,v,loc in self.charts:
                nf = sp.Poly(loc.reduce(f.subs(z,1))[1],*v,domain=sp.QQ)
                col.update({(idx,*a):c for a,c in nf.as_dict().items()})
            cols.append(col)
        return [mon[i] for i in independent(cols)]

    def space(self,positive=(),negative=()):
        h,dp = self.divisor_data(positive)
        _,dm = self.divisor_data(negative)
        n0 = sp.Poly(h,*self.variables).total_degree()
        n = max(n0,(2*self.genus-2)//self.degree+1,0)
        while True:
            mon = self.ambient_sections(n)
            if len(mon) == n*self.degree+1-self.genus:
                break
            n += 1
        # All affine patches are checked, including zeros introduced by the extra factor.
        nonzero_coordinate = self.charts[0][1]
        h = sp.expand(h*nonzero_coordinate**(n-n0))
        cols = [{} for _ in mon]
        for idx,z,variables,loc in self.charts:
            def product(pieces):
                result = gb([1],variables)
                for mult, generators in pieces:
                    factor = gb([*loc,*(sp.sympify(g).subs(z,1) for g in generators)],variables)
                    if unit(factor):
                        continue
                    for _ in range(mult):
                        if unit(result):
                            result = factor
                        else:
                            pieces_g = [loc.reduce(a*b)[1] for a in result for b in factor]
                            result = gb([*loc,*pieces_g],variables)
                return result
            plus,minus = product(positive),product(negative)
            target = gb([*loc,*(h.subs(z,1)*g for g in minus)],variables)
            tests = [loc.reduce(g)[1] for g in plus]
            for j,f in enumerate(mon):
                for k,g in enumerate(tests):
                    if g:
                        rem = sp.Poly(target.reduce(f.subs(z,1)*g)[1],*variables,domain=sp.QQ)
                        cols[j].update({(idx,k,*m):c for m,c in rem.as_dict().items()})
        nums = [sp.expand(sum(c*f for c,f in zip(row,mon))) for row in kernel(cols)]
        return {'curve_degree':self.degree,'genus':self.genus,'degree':dp-dm,'positive_degree':dp,'negative_degree':dm,
                'ambient_dimension':len(mon),'dimension':len(nums),'denominator_degree':n,'denominator':h,'numerators':nums}

    def _symmetric_data(self,m,i,slope):
        n = len(self.variables)-1
        if m < 1 or not 0<=i<=n*m or not 0<=slope<=(n-1)*comb(m,2):
            raise ValueError('Degree, chart index or separating slope is outside the atlas')
        t = sp.Symbol('t')
        a = sp.symbols(f'a0:{m}')
        b = [sp.symbols(f'b{k}_0:{m}') for k in range(2,n+1)]
        ring = sp.QQ.poly_ring(*a,*(q for row in b for q in row))
        def poly(f):
            return sp.Poly(f,t,domain=ring)
        h = poly(t**m+sum(w*t**j for j,w in enumerate(a)))
        gs = [poly(sum(w*t**j for j,w in enumerate(row))) for row in b]
        coords = [poly(t)-sum((slope**(k+1)*g for k,g in enumerate(gs)),poly(0)),*gs]
        coords.insert(0,poly(1)-sum((i**(k+1)*f for k,f in enumerate(coords)),poly(0)))
        def power(f,e):
            return (f**e).rem(h)
        residues = []
        for f in self.equations:
            rem = poly(0)
            for exponents,c in sp.Poly(f,*self.variables,domain=sp.QQ).terms():
                term = poly(c)
                for g,e in zip(coords,exponents):
                    term = (term*power(g,e)).rem(h)
                rem += term
            residues.append([rem.nth(j) for j in range(m)])
        return t, (*a,*(v for row in b for v in row)), h.as_expr(),[c.as_expr() for c in coords],residues

    def symmetric_chart(self,m,i,slope):
        _,_,h,coords,residues = self._symmetric_data(m,i,slope)
        return {'degree':m,'chart_index':i,'separating_slope':slope,'h':str(h),
                'coordinates':list(map(str,coords)),'equations':[list(map(str,row)) for row in residues]}

    def search(self,m,height):
        """All effective rational divisors of the specified canonical atlas height and degree."""
        if m < 1 or height < 1:
            raise ValueError('Positive divisor degree and search height are required')
        n = len(self.variables)-1
        qden = sp.Symbol('den')
        for i in range(n*m+1):
            for s in range((n-1)*comb(m,2)+1):
                t,params,h,coords,equations = self._symmetric_data(m,i,s)
                tests = []
                for row in equations:
                    for f in row:
                        if f:
                            poly = sp.Poly(f,*params,domain=sp.QQ).clear_denoms()[1]
                            homogeneous = poly.homogenize(qden).as_expr()
                            tests.append(sp.lambdify((*params,qden),homogeneous,'math'))
                for denominator in range(1,height+1):
                    for nums in product(range(-height,height+1),repeat=len(params)):
                        if gcd(denominator,*nums)!=1 or any(test(*nums,denominator) for test in tests):
                            continue
                        sub = dict(zip(params,[sp.Rational(k,denominator) for k in nums]))
                        mod = sp.Poly(h.subs(sub),t,domain=sp.QQ)
                        position = [sp.Poly(v.subs(sub),t,domain=sp.QQ).rem(mod) for v in coords]
                        if first_chart(position,mod) == (i,s):
                            yield {'chart_index':i,'separating_slope':s,'denominator':denominator,
                                   'coefficient_numerators':list(nums),'h':str(mod.as_expr()),
                                   'coordinates':[str(v.as_expr()) for v in position]}


def first_chart(coords, modulus):
    """Find the first chart and slope with invertible denominator and a power basis."""
    n,m = len(coords)-1,modulus.degree()
    t = modulus.gens[0]
    zero = sp.Poly(0,t,domain=sp.QQ);one=sp.Poly(1,t,domain=sp.QQ)
    for i in range(n*m+1):
        line = sum((i**k*c for k,c in enumerate(coords)),zero).rem(modulus)
        if line.gcd(modulus).degree():
            continue
        inverse = line.invert(modulus)
        for s in range((n-1)*comb(m,2)+1):
            primitive = (inverse*sum((s**k*c for k,c in enumerate(coords[1:])),zero)).rem(modulus)
            powers=[one]
            for _ in range(m-1):
                powers.append((powers[-1]*primitive).rem(modulus))
            if sp.Matrix([[p.nth(k) for p in powers] for k in range(m)]).det():
                return i,s
    raise ArithmeticError('The finite atlas did not contain this effective degree-m divisor')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input',type=Path)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument('--atlas-degree',type=int)
    mode.add_argument('--search-degree',type=int)
    parser.add_argument('--height',type=int,default=1,help='coefficient height of the first admissible interpolation chart')
    parser.add_argument('--frame-index',type=int,default=0)
    parser.add_argument('--slope',type=int,default=0)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding='utf-8'))
    curve = EmbeddedCurve.from_data(data)
    if args.search_degree is not None:
        result = list(curve.search(args.search_degree,args.height))
    elif args.atlas_degree is None:
        result = curve.space(curve.pieces(data.get('projective_positive',[])),curve.pieces(data.get('projective_negative',[])))
        result['denominator'] = str(result['denominator'])
        result['numerators'] = list(map(str,result['numerators']))
    else:
        result = curve.symmetric_chart(args.atlas_degree,args.frame_index,args.slope)
    print(json.dumps(result,indent=2))
