"""Riemann--Roch spaces on smooth projective plane curves over Q, using finite divisor ideals."""
import argparse
import json
from pathlib import Path

import sympy as sp

X, Y, Z = sp.symbols('X Y Z')
U, V = sp.symbols('u v')
T = sp.Symbol('t')
QQ = sp.QQ


def polynomial(terms, variables):
    """Parse exact [exponents..., coefficient] rows without evaluating source text."""
    result = sp.S.Zero
    for *exponents, coefficient in terms:
        if len(exponents) != len(variables) or any(type(e) is not int or e < 0 for e in exponents):
            raise ValueError('Expected one nonnegative integer exponent per variable')
        result += sp.Rational(coefficient)*sp.prod(x**e for x, e in zip(variables, exponents))
    return sp.expand(result)


def terms(expr, variables):
    return [[*powers, str(c)] for powers, c in sp.Poly(expr, *variables, domain=QQ).terms() if c]


def groebner(generators):
    return sp.groebner(generators, V, T, order='lex', domain=QQ)


def unit(ideal):
    return list(ideal) == [sp.S.One]


def length(ideal):
    if unit(ideal):
        return 0
    if not ideal.is_zero_dimensional:
        raise ValueError('An effective divisor must be a finite scheme on the affine curve')
    leading = [p.LM(order=ideal.order).exponents for p in ideal.polys]
    nv = min(i for i, j in leading if j == 0)
    nt = min(j for i, j in leading if i == 0)
    return sum(not any(a <= i and b <= j for a, b in leading)
               for i in range(nv) for j in range(nt))


def kernel(columns):
    """Regular Q-basis of the kernel of sparse rational columns."""
    rows = sorted({key for col in columns for key in col})
    matrix = sp.zeros(len(rows), len(columns))
    row_index = {key: i for i, key in enumerate(rows)}
    for j, col in enumerate(columns):
        for key, coeff in col.items():
            matrix[row_index[key], j] = coeff
    return matrix.nullspace()


def coefficients(expr):
    return sp.Poly(expr, V, T, domain=QQ).as_dict()


class PlaneCurve:
    """The supplied frame (u:v:w) contains all divisor support in w != 0."""
    def __init__(self, form, frame=None):
        poly = sp.Poly(form, X, Y, Z, domain=QQ)
        self.degree = poly.total_degree()
        if self.degree < 2 or not poly.is_homogeneous:
            raise ValueError('Expected a nonzero homogeneous plane curve of degree at least two')
        self.form = poly.as_expr()
        for var in (X, Y, Z):
            rest = [s for s in (X, Y, Z) if s != var]
            jac = [self.form.diff(s).subs(var, 1) for s in (X, Y, Z)]
            if list(sp.groebner(jac, *rest, domain=QQ)) != [1]:
                raise ValueError('The projective source is singular')
        basis = sp.eye(3) if frame is None else sp.Matrix(frame).applyfunc(sp.Rational)
        if basis.shape != (3, 3) or not basis.det():
            raise ValueError('The projective frame must be an invertible rational 3 by 3 matrix')
        self.frame = basis
        self.frame_forms = tuple(basis*sp.Matrix([X, Y, Z]))
        normalized = self.form.xreplace(dict(zip((X, Y, Z), basis.inv()*sp.Matrix([X, Y, Z])))).expand()
        self.shear = next(s for s in range(self.degree+1) if normalized.subs({X: -s, Y: 1, Z: 0}))
        self.affine = normalized.xreplace({X: T-self.shear*V, Y: V, Z: sp.S.One}).expand()
        self.curve_ideal = groebner([self.affine])
        self.infinity = sp.Poly(normalized.subs(Z, 0).subs({X: 1-self.shear*V, Y: V}), V, domain=QQ)
        self.genus = (self.degree-1)*(self.degree-2)//2

    @classmethod
    def from_data(cls, data):
        return cls(polynomial(data['terms'], (X, Y, Z)), data.get('frame'))

    def monomials(self, degree):
        return [T**(n-j)*V**j for n in range(degree+1) for j in range(min(n, self.degree-1)+1)]

    def homogenize(self, expr, degree):
        a, b, c = self.frame_forms
        return sum(coefficient*(a+self.shear*b)**i*b**j*c**(degree-i-j)
                   for (i, j), coefficient in sp.Poly(expr, T, V, domain=QQ).terms()).expand()

    def ideal(self, pieces=()):
        """Add effective divisors by multiplying their finite ideals on the curve."""
        result = groebner([1])
        for multiplicity, generators in pieces:
            if type(multiplicity) is not int or multiplicity < 1 or not generators:
                raise ValueError('Each divisor piece needs positive multiplicity and nonempty generators')
            gens = [sp.sympify(g).subs(U, T-self.shear*V) for g in generators]
            factor = groebner([self.affine, *gens])
            length(factor)
            if unit(factor):
                continue
            for _ in range(multiplicity):
                if unit(result):
                    result = factor
                else:
                    products = [self.curve_ideal.reduce(a*b)[1] for a in result for b in factor]
                    result = groebner([self.affine, *products])
        return result

    def denominator(self, divisor):
        """Find a lowest-degree section which vanishes on D and misses C at infinity."""
        if unit(divisor):
            return sp.S.One, 0
        upper = min(sp.degree(p, T) for p in divisor if not p.has(V))
        for n in range(1, int(upper)+1):
            mon = self.monomials(n)
            cols = [coefficients(divisor.reduce(m)[1]) for m in mon]
            sections = [sp.expand(sum(a*m for a, m in zip(row, mon))) for row in kernel(cols)]
            leading = [sp.Poly(sum(c*V**j for (i, j), c in sp.Poly(h, T, V).terms() if i+j == n), V,
                               domain=QQ) for h in sections]
            common = self.infinity
            for h in leading:
                common = common.gcd(h)
            if common.degree() > 0:
                continue
            # Each geometric point at infinity rules out at most len(sections)-1 scalar choices.
            for value in range(self.degree*(len(sections)-1)+1):
                weights = [value**j for j in range(len(sections))]
                top = sum((w*h for w, h in zip(weights, leading)), sp.Poly(0, V, domain=QQ))
                if self.infinity.gcd(top).degree() == 0:
                    return sp.expand(sum(w*h for w, h in zip(weights, sections))), n
        raise ArithmeticError('The eliminant was expected to supply an admissible denominator')

    def space(self, positive=(), negative=()):
        """Return a common projective denominator and a full independent Q-basis of numerators."""
        plus, minus = self.ideal(positive), self.ideal(negative)
        h, n = self.denominator(plus)
        target = groebner([self.affine, *(h*g for g in minus)])
        tests = [self.curve_ideal.reduce(g)[1] for g in plus]
        tests = [g for g in tests if g]
        mon = self.monomials(n)
        cols = []
        for m in mon:
            col = {}
            for k, g in enumerate(tests):
                remainder = target.reduce(m*g)[1]
                col.update({(k, *powers): c for powers, c in coefficients(remainder).items()})
            cols.append(col)
        basis = [self.homogenize(sum(a*m for a, m in zip(row, mon)), n) for row in kernel(cols)]
        return {'positive_degree': length(plus), 'negative_degree': length(minus),
                'degree': length(plus)-length(minus), 'dimension': len(basis),
                'denominator_degree': n, 'denominator': self.homogenize(h, n), 'numerators': basis}

    def check_duality(self, positive, negative, result):
        # K=(d-3)H. The line t=0 avoids all points of the curve at infinity.
        h = U+self.shear*V
        m = self.degree-3
        surplus = [(m, [h])] if m > 0 else []
        deficit = [(-m, [h])] if m < 0 else []
        other = self.space(list(negative)+surplus, list(positive)+deficit)
        if result['dimension']-other['dimension'] != result['degree']+1-self.genus:
            raise ArithmeticError('Riemann--Roch duality failed')
        return other['dimension']


def from_json(pieces):
    return [(p.get('multiplicity', 1), [polynomial(g, (U, V)) for g in p['generators']]) for p in pieces]


def output(result):
    return {**result, 'denominator': str(result['denominator']), 'numerators': list(map(str, result['numerators']))}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--duality', action='store_true', help='also compute L(K-D) and verify the dimension identity')
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding='utf-8'))
    if 'projective_positive' in data or 'projective_negative' in data:
        if 'positive' in data or 'negative' in data:
            parser.error('Choose either projective or affine divisor data')
        from symmetric_power import from_projective
        curve, pos, neg = from_projective(data)
    else:
        curve = PlaneCurve.from_data(data)
        pos, neg = from_json(data.get('positive', [])), from_json(data.get('negative', []))
    result = curve.space(pos, neg)
    if args.duality:
        result['dual_dimension'] = curve.check_duality(pos, neg, result)
    print(json.dumps(output(result), indent=2))
