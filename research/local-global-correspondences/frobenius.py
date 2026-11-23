"""Good reduction, genus-two Frobenius and verifiable geometric Hom bounds."""
from functools import lru_cache
from math import isqrt

import numpy as np
import sympy as sp

X, T, U = sp.symbols('x T U')


def polynomial(coefficients):
    """Seven ascending integral coefficients of a squarefree binary sextic."""
    if len(coefficients) not in (6, 7) or any(type(c) != int for c in coefficients):
        raise ValueError('Use six or seven integer coefficients in ascending order')
    f = sp.Poly.from_list(list(reversed(coefficients)), X, domain=sp.QQ)
    if f.degree() not in (5, 6) or f.gcd(f.diff()).degree():
        raise ValueError('A smooth genus-two quintic or sextic is required')
    return f


def good(f, p):
    if not sp.isprime(p) or not 3 <= p <= 1009:
        raise ValueError('Use an odd prime at most 1009 for direct counting')
    c = [int(f.nth(i)) % p for i in range(7)]
    fp = sp.Poly.from_list(c[::-1], X, modulus=p)
    return fp.degree() in (5, 6) and not fp.gcd(fp.diff()).degree()


def character_table(p):
    chi = -np.ones(p, dtype=np.int64)
    chi[0] = 0
    for i in range(1, p):
        chi[i*i % p] = 1
    return chi


def frobenius(f, p):
    """Count both fields by quadratic characters, including the infinity fibre."""
    if not good(f, p):
        raise ValueError(f'Bad reduction at {p}')
    chi = character_table(p)
    c = [int(f.nth(i)) % p for i in range(7)]
    x = np.arange(p, dtype=np.int64)
    v = np.zeros_like(x)
    for a in reversed(c):
        v = (v*x+a) % p
    n1 = int(p+1+chi[v].sum()+chi[c[6]])
    delta = int(np.where(chi == -1)[0][0])
    # Over F_p[u]/(u^2-delta), the quadratic character equals chi of the norm.
    r, s = np.meshgrid(x, x)
    a, b = np.zeros_like(r), np.zeros_like(r)
    for k in reversed(c):
        a, b = (a*r+delta*b*s+k) % p, (a*s+b*r) % p
    norm = (a*a-delta*b*b) % p
    n2 = int(p*p+1+chi[norm].sum()+int(c[6] != 0))
    tr = p+1-n1
    tr2 = p*p+1-n2
    if (tr*tr-tr2) % 2:
        raise ArithmeticError('Inconsistent point counts')
    h = (1, -tr, (tr*tr-tr2)//2, -p*tr, p*p)
    return (n1, n2), h


@lru_cache(None)
def cyclotomic_ratios(h, k):
    """Cyclotomic factors of the cross-root resultant, including multiplicities."""
    f = sp.Poly.from_list(h, T, domain=sp.ZZ).as_expr()
    g = sp.Poly.from_list(k, T, domain=sp.ZZ).as_expr()
    ratio = sp.Poly(sp.resultant(f, g.subs(T, U*T), T), U).primitive()[1]
    return tuple((tuple(map(int, q.all_coeffs())), e)
                 for q,e in ratio.factor_list()[1] if q.is_cyclotomic)


def squareclass(n):
    return int(sp.prod(p for p,k in sp.factorint(abs(n)).items() if k % 2))*(1 if n>0 else -1)


def field_witness(p, h):
    """Sufficient test for a quartic geometric reduction field with one quadratic subfield."""
    f = sp.Poly.from_list(h, T, domain=sp.ZZ)
    if not f.is_irreducible or cyclotomic_ratios(h, h) != (((1, -1), 4),):
        return None
    disc = int(f.discriminant())
    if disc >= 0 and isqrt(disc)**2 == disc:
        return None
    real_disc = h[1]**2-4*(h[2]-2*p)
    if real_disc <= 0 or isqrt(real_disc)**2 == real_disc:
        return None
    return {'prime':p, 'quartic_discriminant':disc, 'real_quadratic_squareclass':squareclass(real_disc)}


def compare(source, target, primes):
    """Bounds valid over Q and its algebraic closure. Absence of a witness is inconclusive."""
    rows = []
    fields = {}
    geometric_zero = rational_zero = None
    for p in primes:
        if not good(source, p) or not good(target, p):
            continue
        c,h = frobenius(source, p)
        d,k = frobenius(target, p)
        roots = cyclotomic_ratios(h, k)
        if geometric_zero is None and not roots:
            geometric_zero = p
        if rational_zero is None and all(q != (1,-1) for q,e in roots):
            rational_zero = p
        witness = field_witness(p, h)
        if witness:
            fields.setdefault(witness['real_quadratic_squareclass'], witness)
        rows.append({'p':p, 'source_counts':c, 'target_counts':d,
                     'source_weil_descending':h, 'target_weil_descending':k,
                     'cyclotomic_root_ratio_factors':roots})
    return {'reductions':rows, 'geometric_hom_zero_prime':geometric_zero,
            'rational_hom_zero_prime':rational_zero,
            'source_endomorphisms_Z_witnesses':list(fields.values())[:2] if len(fields)>=2 else []}
