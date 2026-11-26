"""Principally polarised products C/a1 x C/a2 (a_i proper ideals of orders of K) and their refined Humbert
forms, exactly.  Shioda-Mitani: every abelian surface isogenous to E x E (E with CM by O_K) is isomorphic to
such a product, so every class of the genus closure arises this way."""
import itertools, math, sys
sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
import cypari2
from gluing import K, E_of, pf, q_from_ns, reduce_form, isometric

pari = cypari2.Pari()


def ideals(d, f):
    """Z-bases (pairs (p, q) = p + q w) of representatives of proper ideal classes of O_f = Z + f O_K,
    from reduced forms (a, b, c) of discriminant f^2 d:  a = Z a + Z (-b + f sqrt d)/2,  sqrt d = 2w - d."""
    D = f * f * d
    out = []
    for a in range(1, int(math.isqrt(-D // 3)) + 2):
        for b in range(-a + 1, a + 1):
            if (b * b - D) % (4 * a):
                continue
            c = (b * b - D) // (4 * a)
            if c < a or (c == a and b < 0) or math.gcd(math.gcd(a, b), c) != 1:
                continue
            out.append(((Fr(a), Fr(0)), (Fr(-b - f * d, 2), Fr(f))))
    return out


def ns_general(Kf, B):
    basisv = [(Fr(1), Fr(0), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0), Fr(0)),
              (Fr(0), Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(0), Fr(1))]
    rows = [[E_of(Kf, bv, B[i], B[j]) for bv in basisv] for i, j in itertools.combinations(range(4), 2)]
    den = 1
    for r in rows:
        for x in r:
            den = den * x.denominator // math.gcd(den, x.denominator)
    R = pari.matrix(4, 6, [int(rows[j][i] * den) for i in range(4) for j in range(6)])
    Hr = pari.mathnf(R)
    Dm = pari.matsolve(pari.mattranspose(Hr) / den, pari.matid(4))
    ns = [tuple(Fr(str(Dm[i, j])) for i in range(4)) for j in range(4)]
    s = 1 if pf(Kf, B, (Fr(1), Fr(1), Fr(0), Fr(0))) > 0 else -1      # orientation from a positive H
    P = lambda v: pf(Kf, B, v) * s
    G = [[(P(tuple(a + b for a, b in zip(u, v))) - P(u) - P(v)) if i != j else 2 * P(u)
          for j, v in enumerate(ns)] for i, u in enumerate(ns)]
    assert all(x.denominator == 1 for r in G for x in r)
    return ns, [[int(x) for x in r] for r in G]


def positive(Kf, v):
    a, c, bp, bq = v
    nb = Kf.mul((bp, bq), Kf.conj((bp, bq)))[0]
    return a > 0 and a * c - nb > 0


def principal_polarisations(Kf, B, box=3):
    ns, G = ns_general(Kf, B)
    out = []
    for co in itertools.product(range(-box, box + 1), repeat=4):
        if sum(co[i] * G[i][j] * co[j] for i in range(4) for j in range(4)) != 2:
            continue
        v = tuple(sum(co[k] * ns[k][t] for k in range(4)) for t in range(4))
        if not positive(Kf, v):
            continue
        q = q_from_ns(G, list(co))
        if int(pari.qfminim(pari.matrix(3, 3, [x for r in q for x in r]), 1, 1)[0]) > 0:
            continue                                                    # product polarisation
        out.append({"coords": co, "H": v, "q": reduce_form(q)})
    return out


def products(d, fs):
    Kf = K(d)
    lat = []
    for f in fs:
        for a in ideals(d, f):
            lat.append((f, a))
    for (f1, a1), (f2, a2) in itertools.combinations_with_replacement(lat, 2):
        B = [(a1[0][0], a1[0][1], Fr(0), Fr(0)), (a1[1][0], a1[1][1], Fr(0), Fr(0)),
             (Fr(0), Fr(0), a2[0][0], a2[0][1]), (Fr(0), Fr(0), a2[1][0], a2[1][1])]
        yield (f1, a1, f2, a2), Kf, B


def principal_polarisations_fast(Kf, B, box=11):
    import numpy as np
    ns, G = ns_general(Kf, B)
    Gm = np.array(G, dtype=np.int64)
    r = np.arange(-box, box + 1, dtype=np.int64)
    grid = np.stack(np.meshgrid(r, r, r, r, indexing="ij"), -1).reshape(-1, 4)
    vals = np.einsum("ni,ij,nj->n", grid, Gm, grid)
    out = []
    for co in grid[vals == 2]:
        co = tuple(int(x) for x in co)
        v = tuple(sum(co[k] * ns[k][t] for k in range(4)) for t in range(4))
        if not positive(Kf, v):
            continue
        q = q_from_ns(G, list(co))
        if int(pari.qfminim(pari.matrix(3, 3, [x for r_ in q for x in r_]), 1, 1)[0]) > 0:
            continue
        out.append({"coords": co, "H": v, "q": reduce_form(q)})
    return out
