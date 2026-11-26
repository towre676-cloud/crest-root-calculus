"""Igusa invariants of y^2 = x^6 - a x^4 + b x^2 - 1 as exact polynomials in the dihedral invariants
u = a^3 + b^3, v = a b (they are invariant under (a, b) -> (zeta a, zeta^2 b) and under a <-> b up to the
twist by -1, which fixes even-weight invariants up to sign).  Solved exactly by linear algebra on monomials."""
import json
import cypari2
pari = cypari2.Pari()


def to_uv(P):
    """Express P(a, b) as a polynomial in u = a^3 + b^3, v = a b (exact; raises if impossible)."""
    deg = int(pari.poldegree(P, "a")) + int(pari.poldegree(P, "b"))
    mons = [(i, j) for i in range(deg // 3 + 1) for j in range(deg // 2 + 1) if 3 * i + 2 * j <= deg]
    exps = [pari("(a^3 + b^3)^%d * (a*b)^%d" % (i, j)) for i, j in mons]
    # coefficient vectors over monomials a^p b^q, p, q <= deg
    keys = [(p, q) for p in range(deg + 1) for q in range(deg + 1)]
    def vec(Q):
        return [pari.polcoef(pari.polcoef(Q, p, "a"), q, "b") for p, q in keys]
    M = pari.matrix(len(keys), len(mons), [x for row in zip(*[vec(e) for e in exps]) for x in row])
    c = pari.matinverseimage(M, pari.matrix(len(keys), 1, vec(P)))
    assert int(pari.matsize(c)[1]) == 1, "not a polynomial in u, v"
    expr = sum((c[k, 0] * pari("u^%d * v^%d" % mons[k]) for k in range(len(mons))), pari(0))
    assert pari.substvec(expr, ["u", "v"], [pari("a^3 + b^3"), pari("a*b")]) == P
    return expr


GEN = pari("genus2igusa(x^6 - a*x^4 + b*x^2 - 1)")
J = [to_uv(g) for g in GEN]
json.dump({"J2,J4,J6,J8,J10 in (u,v)": [str(x) for x in J]}, open("dihedral_igusa.json", "w"), indent=1)
for name, x in zip(("J2", "J4", "J6", "J8", "J10"), J):
    print(name, "=", x)
