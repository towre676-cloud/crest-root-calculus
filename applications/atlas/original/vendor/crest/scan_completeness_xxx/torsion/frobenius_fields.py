"""Field degrees sufficient to find every geometric degree-three quotient in a good fibre.

python frobenius_fields.py --universal
python frobenius_fields.py --prime 7 --weil 1 0 0 -21 0 0 343

Coefficients are given in descending order. A supplied Weil polynomial must be
the characteristic polynomial of the Jacobian in the stated characteristic.
"""
import argparse
from collections import Counter
from itertools import product
import json
from math import gcd, lcm
import sympy as sp

T,U = sp.symbols('T U')


def scalar_order(coefficients):
    """First positive n for which T^n equals a nonzero scalar modulo (3,P)."""
    coeffs = [int(c) % 3 for c in coefficients]
    if coeffs[0] != 1 or not coeffs[-1] or len(coeffs) < 2:
        raise ValueError('Expected a positive-degree monic polynomial with invertible constant modulo three.')
    tail = coeffs[:0:-1]
    d = len(tail)
    power = [1] + [0] * (d-1)
    # The class of T is a unit of the finite algebra, so the loop terminates.
    for n in range(1,3**d):
        top = power[-1]
        power = [(old - top*c) % 3 for old,c in zip([0]+power[:-1],tail)]
        if not any(power[1:]) and power[0]:
            return n
    raise ArithmeticError('A unit failed to return to the scalar subgroup.')


def power_mod(poly, exponent):
    power = sp.Poly(T, T, domain=sp.QQ)
    result = sp.Poly(1, T, domain=sp.QQ)
    while exponent:
        if exponent & 1:
            result = (result*power).rem(poly)
        power = (power*power).rem(poly)
        exponent //= 2
    return result


def power_polynomial(poly, exponent):
    reduced = power_mod(poly,exponent)
    return sp.Poly(sp.resultant(poly.as_expr(),U-reduced.as_expr(),T),U,domain=sp.QQ)


def root_ratio_orders(poly):
    """All orders of roots of unity among ratios of eigenvalues."""
    squarefree = poly.sqf_part()
    expr = squarefree.as_expr()
    ratio = sp.Poly(sp.resultant(expr,expr.subs(T,U*T),T),U,domain=sp.QQ)
    orders = []
    for fac,_ in ratio.factor_list()[1]:
        if fac.is_cyclotomic:
            d = fac.degree()
            # n <= 2 phi(n)^2 by multiplicativity and the prime-power formula.
            n = next(n for n in range(1,2*d*d+1)
                     if sp.totient(n)==d and sp.cyclotomic_poly(n,U)==fac.monic().as_expr())
            orders.append(n)
    return sorted(orders)


def integral_coeffs(q):
    if any(not c.is_Integer for c in q.all_coeffs()):
        raise ArithmeticError('An expected integral factor has fractional coefficients.')
    return [int(c) for c in q.all_coeffs()]


def reciprocal_shape(a,b,c,mu):
    return [1,a,b,c,mu*b,mu*mu*a,mu**3]


def universal():
    result = {}
    for mu in (1,2):
        distribution = Counter()
        witness = {}
        for a,b,c in product(range(3),repeat=3):
            coeffs = reciprocal_shape(a,b,c,mu)
            n = scalar_order(coeffs)
            distribution[n] += 1
            witness[n] = [x%3 for x in coeffs]
        maximum = max(distribution)
        result[str(mu)] = {'distribution':dict(sorted(distribution.items())),
                           'maximum_scalar_exponent':maximum,
                           'maximizing_polynomial_mod3':witness[maximum],
                           'maximum_single_field_degree_for_marks':2*maximum}
    assert sum(result['1']['distribution'].values())==sum(result['2']['distribution'].values())==27
    assert (result['1']['maximum_scalar_exponent'],result['2']['maximum_scalar_exponent'])==(36,40)
    return result


def certificate(p,coefficients):
    if not sp.isprime(p) or p<=3:
        raise ValueError('Use a good prime greater than three.')
    if len(coefficients)!=7 or coefficients != reciprocal_shape(*coefficients[1:4],p):
        raise ValueError('Expected a degree-six reciprocal Weil characteristic polynomial.')
    P = sp.Poly.from_list(coefficients,T,domain=sp.QQ)
    orders = root_ratio_orders(P)
    power = lcm(*orders)
    e = scalar_order(coefficients)
    # An even extension exposes all supersingular elliptic powers as rational linear factors.
    split = lcm(2,power)
    Psplit = power_polynomial(P,split)
    factors = P.factor_list()[1]
    components = []
    for q,mult in Psplit.factor_list()[1]:
        ordinary = q.degree()==2 and q.nth(1)%p!=0 and q.nth(0)==p**split
        supersingular = q.degree()==1
        if not ordinary and not supersingular:
            continue
        component = sp.Poly(1,T,domain=sp.QQ)
        for factor,reps in factors:
            image = power_polynomial(factor,split).sqf_part().monic()
            if image==q:
                component *= factor**reps
        assert component.degree()==q.degree()*mult and component.degree()%2==0
        dim = component.degree()//2
        coeff = integral_coeffs(component)
        local_orders = root_ratio_orders(component)
        local_e = scalar_order(coeff)
        bound = gcd(local_e,lcm(*local_orders)) if dim>1 else 1
        assert gcd(power,e)%bound==0
        components.append({'powered_simple_polynomial':integral_coeffs(q),'geometric_multiplicity':dim,
                           'weil_polynomial':coeff,'ratio_orders':local_orders,
                           'mod3_scalar_exponent':local_e,'class_field_divides':bound})
    combined = lcm(*(r['class_field_divides'] for r in components))
    return {'p':p,'weil_polynomial':coefficients,'root_ratio_orders':orders,
            'endomorphism_field_divides':power,'mod3_scalar_exponent':e,
            'whole_jacobian_class_field_divides':gcd(power,e),
            'even_splitting_power':split,'elliptic_components':components,
            'elliptic_dimension':sum(r['geometric_multiplicity'] for r in components),
            'one_sufficient_scan_degree':2*combined if components else None,
            'separate_sufficient_scan_degrees':sorted(set(2*r['class_field_divides'] for r in components))}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('--universal',action='store_true')
    ap.add_argument('--prime',type=int)
    ap.add_argument('--weil',type=int,nargs=7)
    args=ap.parse_args()
    if args.universal:
        out=universal()
    elif args.prime and args.weil:
        out=certificate(args.prime,args.weil)
    else:
        ap.error('Choose --universal or supply both --prime and --weil.')
    print(json.dumps(out,indent=2))
