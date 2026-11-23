"""Replay two coordinate presentations and check geometric exhaustion and an incomplete control."""
from fractions import Fraction
import json
from pathlib import Path

import sympy as sp

from cartier_search import XYZ, algebra, discover, even_count, quartic, read_terms


HERE = Path(__file__).resolve().parent


def verify_family():
    k,x,y,g = sp.symbols('k x y g')
    b = k**2-k+1
    d = x**3+(k-2)*x**2+b*x
    h = 16*g**2-(24*k**2-4*k+41)*g+9*k**4+21*k**3+4*k**2+5*k+25
    quartic_x = (x**2+(k-1)*x+b-g)**2-4*d
    if sp.factor(sp.discriminant(quartic_x,x)+256*g**2*(g-(k-1)**2)*h) != 0:
        raise ArithmeticError('The family branch formula failed')
    if sp.expand(quartic_x.subs(g,0)-(x**2+(k-3)*x+b)**2) != 0:
        raise ArithmeticError('The zero fibre is not even')
    units = [sp.discriminant(d,x), sp.discriminant(x**2+(k-3)*x+b,x),
             sp.discriminant(h,g), h.subs(g,0), h.subs(g,(k-1)**2), k-1]
    if any(v.subs(k,2)%17 == 0 for v in units):
        raise ArithmeticError('The branch sections collide in the progression')
    g2 = x**2+x+3-2*y
    hs = (g*h).subs(k,2).subs(g,g2)
    u = (-3*x**2+2*x*y+y-3)/2
    v = (4*x**3-3*x**2+6*x+5*y+3)/4
    form = read_terms(HERE/'quartic.json')
    relation = sum(c*u**i*v**j*hs**(e//2) for (i,j,e),c in form.items())
    if sp.rem(relation,y**2-x**3-3*x,y).expand() != 0:
        raise ArithmeticError('The family does not specialize to the plane equation')
    linear = 24*k**2-4*k+41
    const = 16*h.subs(g,0)
    invariant = 256*(linear**2-3*const)**3/(const**2*(linear**2-4*const))
    if invariant.subs(k,2) != sp.Rational(17779581,32166277) or invariant.subs(k,2)==invariant.subs(k,19):
        raise ArithmeticError('The varying elliptic target was not verified')


def run():
    verify_family()
    expected = json.loads((HERE/'receipt.json').read_text())
    results = {}
    for name in ('quartic','mixed'):
        path = HERE/f'{name}.json'
        answer = discover(path)
        if json.loads(json.dumps(answer)) != expected[name]:
            raise ArithmeticError(f'The reconstruction receipt differs for {name}')
        if not answer['all_optimal_directions_exhausted']:
            raise ArithmeticError('The geometric bound was not attained')
        invariants = []
        for m in answer['primitive_degree_four']:
            a,b = map(Fraction,m['target'])
            invariants.append(1728*4*a**3/(4*a**3+27*b*b))
        if sorted(invariants) != [Fraction(17779581,32166277),Fraction(11664,13)]:
            raise ArithmeticError('Wrong targets')
        blocks = answer['local_bound']['local_profile']['elliptic_blocks']
        if sorted(abs(b['projector_multiplier']) for b in blocks) != [16,32,32]:
            raise ArithmeticError('Wrong specialization denominator bound')
        if answer['local_bound']['point_counts'] != [30,312,4878]:
            raise ArithmeticError('Wrong point counts')
        if len(answer['double_quotients']) != 1 or len(invariants) != 2:
            raise ArithmeticError('The three directions did not separate')
        if answer['double_quotients'][0]['j'] != '1728':
            raise ArithmeticError('Wrong elliptic double quotient')
        form = read_terms(path)
        even = answer['double_quotients'][0]['even_model']
        for k in (1,2,3):
            if even_count(even,5,k) != algebra.point_count(form,5,k):
                raise ArithmeticError('Square-root and direct projective counts disagree')
        results[name] = answer
    f = quartic(read_terms(HERE/'quartic.json'),[XYZ[0]+XYZ[1], XYZ[1]+XYZ[2], sum(XYZ)])
    if sp.expand(f-quartic(read_terms(HERE/'mixed.json'),XYZ)) != 0:
        raise ArithmeticError('The coordinate change failed')
    limited = discover(HERE/'quartic.json',height=1,prime_limit=31)
    if limited['all_optimal_directions_exhausted'] or limited['primitive_degree_four']:
        raise ArithmeticError('Cutoff incorrectly asserted completeness')
    return results


if __name__=='__main__':
    run()
    print('Verified both positive quartic presentations, both new degree-four maps and all three optimal directions.')
