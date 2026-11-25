"""Richelot (2,2)-isogenous genus-two curves and their correspondences, exactly over Q."""
from itertools import combinations

import sympy as sp

x, z = sp.symbols("x z")


def pairings(roots):
    roots = list(roots)
    if not roots:
        yield []
        return
    a = roots[0]
    for b in roots[1:]:
        rest = [r for r in roots[1:] if r != b]
        for tail in pairings(rest):
            yield [(a, b)] + tail


def richelot(G):
    """G = [G1, G2, G3] quadratics in x.  Returns (Delta, [H1, H2, H3])."""
    cs = [sp.Poly(g, x).all_coeffs()[::-1] + [0] * (3 - len(sp.Poly(g, x).all_coeffs())) for g in G]
    cs = [[sp.Poly(g, x).coeff_monomial(x ** k) for k in range(3)] for g in G]
    Delta = sp.Matrix(cs).det()
    d = [sp.diff(g, x) for g in G]
    H = [sp.expand(d[(i + 1) % 3] * G[(i + 2) % 3] - G[(i + 1) % 3] * d[(i + 2) % 3]) for i in range(3)]
    return sp.nsimplify(Delta), H


if __name__ == "__main__":
    f = sp.expand((x**2 - 1) * (x**2 - 4) * (x**2 - 9))
    target = sp.Poly(11 * (x**2 - 10*x + 13) * (2*x**2 + 2*x - 7) * (7*x**2 + 18*x + 3), x)
    for P in pairings([1, -1, 2, -2, 3, -3]):
        G = [sp.expand((x - a) * (x - b)) for a, b in P]
        D, H = richelot(G)
        if D == 0:
            print(P, "Delta=0 (product)")
            continue
        rhs = sp.Poly(sp.expand(H[0] * H[1] * H[2] / D), x)
        ratio = sp.simplify(target.as_expr() / rhs.as_expr())
        print(P, "Delta", D, "ratio", ratio if ratio.is_number else "-")
