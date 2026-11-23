"""Charts for symmetric powers and exhaustive bounded chord searches over Q."""
import argparse
from functools import reduce
from itertools import product
import json
from math import comb, gcd
from pathlib import Path

import sympy as sp

from riemann_roch import PlaneCurve, T, U, V, X, Y, Z, polynomial, terms

P = (X, Y, Z)
A, B = sp.symbols('A B')


def moment_frame(i):
    return sp.Matrix([[0, 1, 0], [0, 0, 1], [1, i, i*i]])


def div_data(groups):
    return [(group.get('multiplicity', 1), [polynomial(g, P) for g in group['generators']]) for group in groups]


def affine_chart(form, groups, specified=None):
    """Choose a rational line at infinity disjoint from the combined projective support."""
    curve = sp.Poly(form, *P, domain=sp.QQ)
    bound = 0
    for _, generators in groups:
        if any(not sp.Poly(g, *P, domain=sp.QQ).is_homogeneous for g in generators if g):
            raise ValueError('Projective generators must be homogeneous')
        degrees = [sp.Poly(g, *P).total_degree() for g in generators if sp.reduced(g, [form], *P)[1]]
        if not degrees:
            raise ValueError('Positive-dimensional projective divisor input')
        bound += curve.total_degree()*min(degrees)

    def avoids(frame):
        if frame.shape != (3, 3) or not frame.det():
            raise ValueError('The specified projective frame is singular')
        line = dict(zip(P, frame.inv()*sp.Matrix([A, B, 0])))
        f = form.xreplace(line)
        for _, generators in groups:
            restrictions = [f, *(g.xreplace(line) for g in generators)]
            if all(poly.subs({A: 1, B: 0}) == 0 for poly in restrictions):
                return False
            gcd_poly = reduce(sp.gcd, (sp.Poly(poly.subs(B, 1), A, domain=sp.QQ) for poly in restrictions))
            if gcd_poly.degree() > 0:
                return False
        return True

    if specified is not None:
        shape = sp.Matrix(specified).applyfunc(sp.Rational)
        if not avoids(shape):
            raise ValueError('The specified line at infinity contains divisor support')
        return shape
    if avoids(sp.eye(3)):
        return sp.eye(3)
    for i in range(2*bound+1):
        shape = moment_frame(i)
        if avoids(shape):
            return shape
    raise ArithmeticError('No disjoint projective chart found within the finite cover')


def affine_pieces(frame, groups):
    affine = dict(zip(P, frame.inv()*sp.Matrix([U, V, 1])))
    return [(m, [g.xreplace(affine).expand() for g in gs]) for m, gs in groups]


def from_projective(data):
    form = polynomial(data['terms'], P)
    pos = div_data(data.get('projective_positive', []))
    neg = div_data(data.get('projective_negative', []))
    frame = affine_chart(form, pos+neg, data.get('frame'))
    return PlaneCurve(form, frame), affine_pieces(frame, pos), affine_pieces(frame, neg)


def encode(curve, pieces):
    """Write a supplied effective divisor as a monic h and an interpolant g on one chart."""
    from riemann_roch import length
    original = curve.ideal(pieces)
    m = length(original)
    if m < 1:
        raise ValueError('Choose a nonempty effective divisor')
    for s in range(comb(m, 2)+1):
        gens = [p.subs(T, T+(curve.shear-s)*V) for p in original]
        basis = sp.groebner(gens, V, T, order='lex', domain=sp.QQ)
        h = next((p for p in basis if not p.has(V) and sp.degree(p, T) == m), None)
        if h is not None:
            g = basis.reduce(V)[1]
            if g.has(V):
                raise ArithmeticError('A separating element did not generate the divisor algebra')
            return s, sp.Poly(h, T).monic().as_expr(), g
    raise ArithmeticError('A curvilinear divisor failed the separating-slope bound')


def chart_equations(curve, m, frame_index, slope):
    """Coefficient equations for h dividing F(T-s*g(T), g(T), 1) in the selected frame."""
    if m < 1 or not 0 <= frame_index <= 2*m or not 0 <= slope <= comb(m, 2):
        raise ValueError('Invalid degree, chart index or projection slope')
    frame = moment_frame(frame_index)
    a, b = sp.symbols(f'a0:{m}'), sp.symbols(f'b0:{m}')
    ring = sp.QQ.poly_ring(*a, *b)
    h = sp.Poly(T**m+sum(a[i]*T**i for i in range(m)), T, domain=ring)
    g = sp.Poly(sum(b[i]*T**i for i in range(m)), T, domain=ring)
    first = sp.Poly(T,T,domain=ring)-slope*g
    affine = dict(zip(P, frame.inv()*sp.Matrix([U, V, 1])))
    f = sp.Poly(curve.form.xreplace(affine), U, V, domain=sp.QQ)
    one = sp.Poly(1, T, domain=ring)
    up, vp = [one], [one]
    for _ in range(curve.degree):
        up.append((up[-1]*first).rem(h))
        vp.append((vp[-1]*g).rem(h))
    residue = sum((c*(up[i]*vp[j]).rem(h) for (i, j), c in f.terms()), sp.Poly(0, T, domain=ring))
    return {'frame':frame.tolist(), 'slope':slope, 'h':str(h.as_expr()), 'g':str(g.as_expr()),
            'equations': [str(residue.nth(i)) for i in range(m)]}


def binary_degree_two(form):
    """All rational degree-two subdivisors of a binary form, retaining multiplicities."""
    factors = sp.factor_list(form, A, B)[1]
    linear = []
    result = []
    for h, multiplicity in factors:
        degree = sp.Poly(h, A, B).total_degree()
        if degree == 1:
            linear.append(h)
            if multiplicity >= 2:
                result.append(h*h)
        elif degree == 2:
            result.append(h)
    result.extend(g*h for i, g in enumerate(linear) for h in linear[i+1:])
    return result


def chords(curve, height):
    """Enumerate all rational length-two divisors with primitive span coefficients <= height."""
    if height < 1:
        raise ValueError('Positive chord height required')
    for values in product(range(-height, height+1), repeat=3):
        if gcd(*values) != 1 or next(c for c in values if c) < 0:
            continue
        line = sum(c*x for c, x in zip(values, P))
        pivot = next(i for i, c in enumerate(values) if c)
        free = [j for j in range(3) if j != pivot]
        chart = {P[free[0]]: A, P[free[1]]: B, P[pivot]: -(values[free[0]]*A+values[free[1]]*B)/values[pivot]}
        for factor in binary_degree_two(curve.form.xreplace(chart)):
            equation = factor.xreplace({A: P[free[0]], B: P[free[1]]}).expand()
            yield values, line, equation


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    choice = parser.add_mutually_exclusive_group(required=True)
    choice.add_argument('--chord-height', type=int)
    choice.add_argument('--atlas-degree', type=int)
    choice.add_argument('--encode', action='store_true')
    parser.add_argument('--frame-index', type=int, default=0)
    parser.add_argument('--slope', type=int, default=0)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding='utf-8'))
    curve = PlaneCurve.from_data(data)
    if args.chord_height is not None:
        answer = [{'line': list(line), 'projective_generators': [terms(L, P), terms(g, P)]}
                  for line, L, g in chords(curve, args.chord_height)]
    elif args.encode:
        curve, pos, _ = from_projective(data)
        s, h, g = encode(curve, pos)
        answer = {'frame':curve.frame.tolist(), 'slope':s, 'h':str(h), 'g':str(g)}
    else:
        answer = chart_equations(curve, args.atlas_degree, args.frame_index, args.slope)
    print(json.dumps(answer, indent=2, default=str))
