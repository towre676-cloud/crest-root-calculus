"""Certified-map records for bielliptic curves y^2 = prod (D x^2 - n_i)."""
from fractions import Fraction
from math import lcm

import sympy as sp

x, u = sp.symbols("x u")


def normalise(h, U, W):
    """Map (u, Y) = (U(x), W(x) y) onto Y^2 = h(u); return N, D, a, b, g2, g3 for X = h3 U + h2/3."""
    hp = sp.Poly(h, u)
    h3, h2 = hp.coeff_monomial(u ** 3), hp.coeff_monomial(u ** 2)
    w = sp.symbols("w")
    expr = sp.expand(4 * h3 ** 2 * h.subs(u, (w - h2 / 3) / h3))
    P = sp.Poly(expr, w)
    assert P.coeff_monomial(w ** 3) == 4 and P.coeff_monomial(w ** 2) == 0
    g2, g3 = -P.coeff_monomial(w), -P.coeff_monomial(1)
    X = sp.cancel(h3 * U + h2 / 3)
    N, D = sp.fraction(X)
    ab = sp.simplify(sp.diff(U, x) / (2 * W))
    abp = sp.Poly(sp.expand(ab), x)
    assert abp.degree() <= 1
    fr = lambda z: Fraction(int(sp.fraction(sp.nsimplify(z))[0]), int(sp.fraction(sp.nsimplify(z))[1]))
    cl = lambda p: [(fr(c), Fraction(0)) for c in sp.Poly(p, x).all_coeffs()[::-1]]
    return {"N": cl(N), "D": cl(D), "a": (fr(abp.coeff_monomial(1)), Fraction(0)),
            "b": (fr(abp.coeff_monomial(x)), Fraction(0)), "g2": (fr(g2), Fraction(0)), "g3": (fr(g3), Fraction(0)),
            "degree": 2, "gaussian": False}


def bielliptic(cs):
    cs = [Fraction(c) for c in cs]
    Dn = lcm(*[c.denominator for c in cs])
    ns = [int(c * Dn) for c in cs]
    f = sp.Poly(sp.expand(sp.prod([Dn * x ** 2 - n for n in ns])), x)
    fco = [int(c) for c in f.all_coeffs()[::-1]]
    m1 = normalise(sp.expand(sp.prod([u - n for n in ns])), Dn * x ** 2, sp.Integer(1))
    m2 = normalise(sp.expand(sp.prod([Dn - n * u for n in ns])), 1 / x ** 2, 1 / x ** 3)
    for m in (m1, m2):
        N, D = (sp.Poly(sum(sp.Rational(c[0])*x**i for i, c in enumerate(m[key])), x)
                for key in ("N", "D"))
        a, b, g2, g3 = (sp.Rational(m[key][0]) for key in ("a", "b", "g2", "g3"))
        W = N.diff()*D-N*D.diff()
        linear = sp.Poly(a+b*x, x)
        if (N.gcd(D).degree() != 0 or max(N.degree(), D.degree()) != 2
                or W.is_zero or not g2**3-27*g3**2
                or f*W**2 != linear**2*D*(4*N**3-g2*N*D**2-g3*D**3)):
            raise ArithmeticError("The bielliptic map does not verify over Q.")
    return fco, [m1, m2]
