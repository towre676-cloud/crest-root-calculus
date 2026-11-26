"""Verify the explicit bielliptic equations for the opposite-secant family."""

from fractions import Fraction
from itertools import combinations, product
from pathlib import Path
import sys

from cypari2 import Pari
from flint import acb, acb_mat, acb_poly, arb, ctx

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "twin-curves"))

from certify_link import gram, humbert_from_gluing, isometry, require

PARI = Pari(128_000_000)


def verify_parity(shift, shear):
    """The gluing conditions depend only on the product form modulo four."""
    for p, q, a, b, c in product(range(4), repeat=5):
        T = ((a, shift*a-c*b), (b, 225*a-shift*b))
        E = [
            [0, p, -T[1][0], T[0][0]],
            [-p, 0, -T[1][1], T[0][1]],
            [T[1][0], T[1][1], 0, q],
            [-T[0][0], -T[0][1], -q, 0],
        ]
        W = ((1, 0, 1, 0), (0, 1, shear, 1))
        EW = [[sum(row[i]*v[i] for i in range(4)) for v in W] for row in E]
        mixed = all(entry % 2 == 0 for row in EW for entry in row)
        paired = sum(W[0][i]*EW[i][1] for i in range(4)) % 4 == 0
        expected = p == q and (p-a) % 2 == 0 and b % 2 == 0
        require((mixed and paired) == expected, "Gluing congruences differ from the uniform formula.")


def target(f, side):
    k = (4*f+1)//25
    shift = (114, 39)[side]
    extra = (231, 27)[side]
    c = (k+extra)//4
    require(225*c-shift*shift == 9*f, "The ideal has the wrong norm form.")
    return shift, [[225, 0, (3, 78)[side]], [0, 4, 0], [(3, 78)[side], 0, k+27*side]]


def certify_forms(f):
    field = -4*f
    order = ((Fraction(1), Fraction(0)), (Fraction(0), Fraction(3)))
    targets = []
    for side, perm in enumerate(((0, 1, 2), (0, 2, 1))):
        shift, expected = target(f, side)
        ideal = ((Fraction(225), Fraction(0)), (Fraction(6*f-shift), Fraction(3)))
        computed = humbert_from_gluing(order, ideal, perm, field)["humbert_form"]
        isometry(gram(computed), gram(expected))
        targets.append(gram(expected))
    require(not PARI.qfisom(*targets), "The two target lattices must be distinct.")
    return targets


def labelled_roots(tau):
    j = tau.modular_j()
    c = j/(j-1728)
    polynomial = acb_poly([-2*c, -3*c, 0, 1])
    g2, g3 = acb.elliptic_invariants(tau)
    scale = 2*g2/(3*g3)
    roots = [scale*acb.elliptic_p(w, tau) for w in (acb(1)/2, tau/2, (1+tau)/2)]
    require(all(polynomial(root).contains(0) for root in roots)
            and all(not a.overlaps(b) for a,b in combinations(roots, 2)),
            "The labelled cubic roots are not isolated.")
    return j, roots


def coefficient_enclosures(f, side):
    shift, _ = target(f, side)
    tau = acb(0, 3)*arb(f).sqrt()
    j_left, left = labelled_roots(tau)
    j_right, right = labelled_roots((tau-shift)/225)
    require(not j_left.overlaps(j_right), "The two elliptic curves may be isomorphic.")
    perm = ((0, 1, 2), (0, 2, 1))[side]
    def chart(roots):
        u, v, w = roots
        return acb_mat([[v-w, -u*(v-w)], [v-u, -w*(v-u)]])

    mobius = chart([right[i] for i in perm]).inv()*chart(left)
    a, b, c, d = mobius.entries()
    require(not c.contains(0) and not (a*d-b*c).contains(0),
            "The root bijection is not certified projective with a finite pole.")
    u = [v+d/c for v in left]
    require(all(not v.contains(0) for v in u) and all(not v.overlaps(w) for v,w in combinations(u, 2)),
            "The derived even sextic is not certified smooth.")
    for v,i in zip(u, perm, strict=True):
        require((a/c+(b*c-a*d)/(c*c*v)).overlaps(right[i]), "The reciprocal labels disagree.")
    sextic = acb_poly([1])
    for v in u:
        sextic *= acb_poly([-v, 0, 1])
    return (j_left, j_right), sextic


def main():
    for shift, shear in ((114, 0), (39, 1)):
        verify_parity(shift, shear)
    for m in (9, 33, 57, 81):
        f = (25*m*m-101)//4
        require(f % 300 == 181, "The parameter is outside the proved family.")
        forms = certify_forms(f)
        rows = [list(map(int, Q.qfrep(m*m))) for Q in forms]
        first = next(d for d in range(1, m+1) if rows[0][d*d-1] != rows[1][d*d-1])
        require(first == m and 2*(rows[0][m*m-1]-rows[1][m*m-1]) == 4,
                "The exact representation counts do not reproduce the first difference.")
        with ctx.workprec(150*m+900):
            for side in range(2):
                js, sextic = coefficient_enclosures(f, side)
                require(sextic.degree() == 6 and all(sextic[i].is_finite() for i in range(7)),
                        "The algebraic equation did not yield finite coefficient enclosures.")
                if m == 9:
                    class_polynomial = PARI.polclass(-36*f)
                    H = acb_poly([int(class_polynomial.polcoef(k)) for k in range(class_polynomial.poldegree()+1)])
                    require(all(H(j).contains(0) for j in js), "A ring-class cross-check failed.")
        print(f"Verified both sextic constructions at M={m}. First difference: +4 classes at degree {m}.",
              flush=True)


if __name__ == "__main__":
    main()
