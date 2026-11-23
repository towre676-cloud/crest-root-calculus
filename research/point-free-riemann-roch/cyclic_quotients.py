"""Exhaust rationally defined geometrically cyclic cubic genus-one quotients of a plane quartic."""
import argparse
from math import gcd
from pathlib import Path
import json

import sympy as sp

from riemann_roch import PlaneCurve, T, U, V, X, Y, Z, terms
from symmetric_power import A, B, P, affine_chart, affine_pieces

R, S = sp.symbols('r s')
H = sp.Symbol('h')


def square_tests(f):
    """Coefficients of obstructions for a binary quartic with nonzero first coefficient."""
    f = sp.Poly(f, A, B)
    a, b, c, d, e = [f.coeff_monomial(A**(4-i)*B**i) for i in range(5)]
    k = 4*a*c-b*b
    return a, sp.expand(8*a*a*d-b*k), sp.expand(64*a**3*e-k*k)


def rational_roots(f, variable):
    p = sp.Poly(f, variable, domain=sp.QQ)
    if p.is_zero:
        raise ArithmeticError('A zero elimination polynomial cannot bound the lines')
    return sorted(p.ground_roots())


def one_parameter(f):
    """Handle a pencil, including members that vanish at the default binary anchor."""
    candidate = None
    for k in range(5):
        _, g, j = square_tests(f.subs(B, k*A+B))
        for eq in (g, j):
            eq = sp.Poly(eq, S, domain=sp.QQ)
            if not eq.is_zero:
                candidate = eq if candidate is None else candidate.gcd(eq)
        if candidate is not None and candidate.degree() == 0:
            break
    if candidate is None:
        raise ArithmeticError('A pencil would consist entirely of bitangents')
    return rational_roots(candidate.as_expr(), S)


def normalize_line(line):
    vals = [sp.Rational(c) for c in line]
    common = sp.ilcm(*(c.q for c in vals))
    vals = [int(c*common) for c in vals]
    g = gcd(*vals)
    if next(c for c in vals if c) < 0:
        g = -g
    return tuple(c//g for c in vals)


def line_form(curve, line):
    pivot = next(i for i, c in enumerate(line) if c)
    free = [i for i in range(3) if i != pivot]
    subs = {P[free[0]]: A, P[free[1]]: B,
            P[pivot]: (H-line[free[0]]*A-line[free[1]]*B)/line[pivot]}
    return sp.Poly(curve.form.xreplace(subs), H, domain=sp.QQ.poly_ring(A, B)), free


def root_quadratic(binary):
    c, factors = sp.factor_list(binary, A, B)
    if any(m % 2 for _, m in factors):
        return None
    return c, sp.expand(sp.prod(g**(m//2) for g, m in factors))


def rational_bitangents(curve):
    if curve.degree != 4:
        raise ValueError('The cyclic genus-one routine requires a smooth plane quartic')
    f = curve.form.xreplace({X: -R*A-S*B, Y: A, Z: B})
    a, g, j = square_tests(f)
    if a.has(S) or not a:
        raise ArithmeticError('Unexpected binary anchor')
    nonzero = [p for p in (g, j) if p]
    constant = next((p for p in nonzero if not p.has(S)), None)
    if constant is not None:
        elim = constant
    elif len(nonzero) == 2:
        elim = sp.resultant(g, j, S)
    else:
        raise ArithmeticError('An infinite generic family of square binary sections was detected')
    candidates = set()
    for r in rational_roots(elim, R):
        if a.subs(R, r):
            polys = [sp.Poly(p.subs(R, r), S, domain=sp.QQ) for p in (g, j)]
            common = polys[0].gcd(polys[1])
            for s in rational_roots(common.as_expr(), S):
                candidates.add(normalize_line((1, r, s)))
    for r in rational_roots(a, R):
        for s in one_parameter(f.subs(R, r)):
            candidates.add(normalize_line((1, r, s)))
    for s in one_parameter(curve.form.xreplace({X: A, Y: -S*B, Z: B})):
        candidates.add(normalize_line((0, 1, s)))
    candidates.add((0, 0, 1))
    answers = []
    for line in sorted(candidates):
        binary = line_form(curve, line)[0].nth(0)
        if root_quadratic(binary) is not None:
            answers.append(line)
    return answers


def fixed_functions(curve, line, with_spaces=True):
    """Recover the unique cyclic quotient, if one exists for this bitangent."""
    section, free = line_form(curve, line)
    root = root_quadratic(section.nth(0))
    if root is None:
        return None
    c, q = root
    quad = sp.Poly(q, A, B)
    a, b, d = [quad.coeff_monomial(m) for m in (A*A, A*B, B*B)]
    disc = b*b-4*a*d
    if not disc:
        return None

    def laplace(p):
        return d*p.diff(A,2)-b*p.diff(A,B)+a*p.diff(B,2)

    z0, z1 = sp.symbols('z0 z1')
    normal = section.nth(1)+2*c*q*(z0*q.diff(A)+z1*q.diff(B))
    lap = sp.Poly(laplace(normal), A, B)
    linear, rhs = sp.linear_eq_to_matrix([lap.coeff_monomial(A), lap.coeff_monomial(B)], (z0,z1))
    if not linear.det():
        raise ArithmeticError('A nonsingular quadratic contact pair must determine a unique center')
    center = linear.inv()*rhs
    centered = sp.Poly(section.as_expr().xreplace({A:A+center[0]*H, B:B+center[1]*H}), H, domain=sp.QQ.poly_ring(A,B))
    quart, cubic, square, first, e = [centered.nth(i) for i in range(5)]
    d0 = sp.cancel(square/q)
    if first or d0.free_symbols:
        return None
    if laplace(cubic).expand() or (quart-c*q*q).expand() or not e:
        raise ArithmeticError('Invalid centered expansion')
    alt = (q.diff(B)*cubic.diff(A)-q.diff(A)*cubic.diff(B))/3
    b0 = sp.cancel((alt**2-disc*cubic**2)/q**3)
    if b0.free_symbols or not b0:
        raise ArithmeticError('The harmonic cubic must give a nonzero multiple of q^3')
    target = sp.Poly(disc*(c*T*T+d0*T+e)**2+b0*T**3, T, domain=sp.QQ)
    if target.degree() != 4 or target.gcd(target.diff()).degree() != 0:
        raise ArithmeticError('The genus-one target lacks a squarefree quartic model')
    L = sum(k*x for k,x in zip(line, P))
    affine_change = {A:P[free[0]]-center[0]*L, B:P[free[1]]-center[1]*L}
    t_num = q.xreplace(affine_change).expand()
    w_num = alt.xreplace(affine_change).expand()
    # Clear L^8 in w^2 = target(t) and certify the identity on the whole source.
    identity = L**2*w_num**2-sum(k*t_num**i*L**(8-2*i) for (i,),k in target.terms())
    if sp.reduced(identity, [curve.form], *P)[1]:
        raise ArithmeticError('The recovered map identity failed')
    result = {'line': list(line), 'contact_quadratic': str(q.xreplace({A:P[free[0]], B:P[free[1]]})),
              'center_in_line_frame':list(map(str,center)), 'deck_extension_squareclass':str(sp.cancel(-3/disc)),
              'rational_deck_transformations': sp.sqrt(sp.cancel(-3/disc)).is_rational is True,
              't_numerator':str(t_num), 'w_numerator':str(w_num), 'common_linear_form':str(L),
              'target_quartic':str(target.as_expr()), 'quotient_degree':3}
    if with_spaces:
        projective_q = q.xreplace({A:P[free[0]], B:P[free[1]]})
        groups = [(1,[L, projective_q])]
        frame = affine_chart(curve.form, groups)
        model = curve if frame == curve.frame else PlaneCurve(curve.form, frame)
        pieces = affine_pieces(frame, groups)[0][1]
        rr_dims = []
        for m, n, base in ((3,2,t_num),(6,3,w_num)):
            rr = model.space([(m,pieces)])
            if rr['denominator_degree'] != n:
                raise ArithmeticError('Unexpected polar degree')
            scale = sp.cancel(rr['denominator']/L**n)
            if scale.free_symbols:
                raise ArithmeticError('The minimal denominator must be proportional to the bitangent power')
            cols = [sp.Poly(p,*P).as_dict() for p in rr['numerators']+[scale*base]]
            all_mons = sorted(set().union(*cols))
            rank = sp.Matrix([[p.get(monomial,0) for p in cols] for monomial in all_mons]).rank()
            if rank != rr['dimension']:
                raise ArithmeticError('The computed function is outside its certified Riemann--Roch space')
            rr_dims.append(rank)
        result['dimensions_L3D_L6D'] = rr_dims
    return result


def recover(curve, with_spaces=True):
    lines = rational_bitangents(curve)
    maps = [fixed_functions(curve, l, with_spaces) for l in lines]
    return {'rational_bitangents': [list(l) for l in lines], 'cyclic_cubic_quotients': [m for m in maps if m is not None]}


if __name__ == '__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('input', type=Path)
    parser.add_argument('--skip-spaces', action='store_true', help='omit the independent Riemann--Roch membership checks')
    args = parser.parse_args()
    curve = PlaneCurve.from_data(json.loads(args.input.read_text(encoding='utf-8')))
    print(json.dumps(recover(curve, not args.skip_spaces), indent=2))
