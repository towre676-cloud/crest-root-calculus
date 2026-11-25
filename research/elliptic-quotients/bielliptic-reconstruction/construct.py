"""Ground-truth test curves: glue E1 to a k-isogenous E2 along each of the six 2-torsion bijections.

E1 : y^2 = x(x+32)(x+27) has a rational 3-isogeny with kernel polynomial x+36.  For each
bijection beta between the roots of E1 and the roots of E2 = E1/<ker>, the genus-two curve
y^2 = prod (x^2 - c_i), c_i = a_i + s, has quotients E1 (by x -> -x) and a twist of E2 with
gluing a_i <-> b_beta(i).  The ground-truth permutation pi = beta^{-1} rho is computed with
PARI's exact Velu isogeny.  Each curve is then scrambled by a rational Mobius change of x so
that its involution must be rediscovered.
"""
import itertools
import json
from fractions import Fraction

import cypari2

pari = cypari2.Pari()
pari.allocatemem(10 ** 9)


def cycle_type(perm):
    seen, cyc = set(), []
    for i in range(len(perm)):
        if i in seen:
            continue
        j, n = i, 0
        while j not in seen:
            seen.add(j)
            j = perm[j]
            n += 1
        cyc.append(n)
    return {(1, 1, 1): "identity", (1, 2): "transposition", (3,): "3-cycle"}[tuple(sorted(cyc))]


def _roots(poly):
    out = []
    F = pari.factor(poly)
    for i in range(len(F[0])):
        q = F[0][i]
        if pari.poldegree(q) == 1:
            r = -pari.polcoef(q, 0) / pari.polcoef(q, 1)
            out.append(Fraction(int(pari.numerator(r)), int(pari.denominator(r))))
    return out


def two_torsion_roots(E):
    """Roots of 4x^3 + b2 x^2 + 2 b4 x + b6 (x-coordinates of 2-torsion)."""
    b2, b4, b6 = E[5], E[6], E[7]
    return _roots(pari(f"4*x^3 + ({b2})*x^2 + 2*({b4})*x + ({b6})"))


def build(E1coeffs, kernel, k, scramble=(2, 1, 1, 3)):
    E1 = pari.ellinit(E1coeffs)
    iso = pari.ellisogeny(E1, pari(kernel))
    E2 = pari.ellinit(iso[0])
    xmap = iso[1][0] / iso[1][2] ** 2
    a = two_torsion_roots(E1)
    b = two_torsion_roots(E2)
    assert len(a) == 3 and len(b) == 3
    rho = []
    for ai in a:
        img = pari.subst(xmap, "x", pari(f"{ai.numerator}/{ai.denominator}"))
        img = Fraction(int(pari.numerator(img)), int(pari.denominator(img)))
        rho.append(b.index(img))
    curves = []
    for beta in itertools.permutations(range(3)):
        bb = [b[beta[i]] for i in range(3)]
        R = (bb[1] - bb[0]) / (bb[2] - bb[0])
        den = (a[0] - a[1]) - R * (a[0] - a[2])
        if den == 0:
            continue
        s = (R * (a[0] - a[2]) * a[1] - (a[0] - a[1]) * a[2]) / den
        c = [ai + s for ai in a]
        if 0 in c or len(set(c)) < 3:
            continue
        # pi = beta^{-1} rho on E1 labels
        binv = [beta.index(j) for j in range(3)]
        pi = [binv[rho[i]] for i in range(3)]
        # f(x) = prod (x^2 - c_i), scrambled by x -> (p x + q)/(r x + t)
        p, q, r, t = scramble
        X = pari(f"({p}*x+{q})")
        D = pari(f"({r}*x+{t})")
        f = pari(1)
        for ci in c:
            f *= X ** 2 - pari(f"{ci.numerator}/{ci.denominator}") * D ** 2
        f = f / pari.content(f)
        curves.append({"beta": list(beta), "pi": pi, "type": cycle_type(pi), "k": k, "shift": str(s),
                       "unscrambled_c": [str(ci) for ci in c], "f": [int(pari.polcoef(f, i)) for i in range(7)]})
    return curves


if __name__ == "__main__":
    curves = build([0, 59, 0, 864, 0], "x+36", 3)
    json.dump(curves, open("test_curves.json", "w"), indent=1)
    for C in curves:
        print(C["beta"], C["type"], C["f"])
