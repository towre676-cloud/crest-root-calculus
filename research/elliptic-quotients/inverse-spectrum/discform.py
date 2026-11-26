"""Exact isomorphism test for discriminant quadratic forms of even lattices (Nikulin): genus test."""
import itertools
import math
from fractions import Fraction as Fr

import cypari2

pari = cypari2.Pari()


def data(q):
    """Smith invariants e_i > 1, Gram b(g_i, g_j) in Q/Z and q(g_i) in Q/2Z for the lattice L = (Z^3, 2q)."""
    n = len(q)
    A = pari.matrix(n, n, [2 * x for r in q for x in r])
    U, V, D = pari.matsnf(A, 1)
    # U A V = D.  L^# = A^{-1} Z^n; generators A^{-1} U^{-1} e_i of order D_ii; in L^# coordinates
    Ainv = pari.matsolve(A, pari.matid(n))
    Ui = pari.matsolve(U, pari.matid(n))
    gens, invs = [], []
    for i in range(n):
        e = int(D[i, i])
        if e > 1:
            g = Ainv * pari.matrix(n, 1, [Ui[k, i] for k in range(n)])
            gens.append([Fr(str(g[k, 0])) for k in range(n)])
            invs.append(e)
    Ai = [[2 * x for x in r] for r in q]
    bil = lambda x, y: sum(x[i] * Ai[i][j] * y[j] for i in range(n) for j in range(n))
    m = len(gens)
    Bm = [[bil(gens[i], gens[j]) % 1 for j in range(m)] for i in range(m)]
    Qm = [bil(gens[i], gens[i]) % 2 for i in range(m)]
    return invs, Bm, Qm


def _elements(invs):
    return list(itertools.product(*[range(e) for e in invs]))


def _q(co, B, Q):
    v = sum(c * c * Q[i] for i, c in enumerate(co))
    v += 2 * sum(co[i] * co[j] * B[i][j] for i in range(len(co)) for j in range(i + 1, len(co)))
    return v % 2


def _b(x, y, B):
    return sum(x[i] * y[j] * B[i][j] for i in range(len(x)) for j in range(len(y))) % 1


def _order(co, invs):
    o = 1
    for c, e in zip(co, invs):
        o = o * (e // math.gcd(e, c)) // math.gcd(o, e // math.gcd(e, c))
    return o


def isomorphic(q1, q2):
    i1, B1, Q1 = data(q1)
    i2, B2, Q2 = data(q2)
    if math.prod(i1) != math.prod(i2) or sorted(i1) != sorted(i2):
        return False
    E2 = _elements(i2)
    tgt = [[y for y in E2 if _order(y, i2) == e and _q(y, B2, Q2) == Q1[k]] for k, e in enumerate(i1)]
    m = len(i1)

    def surj(imgs):
        seen = set()
        for co in _elements(i1):
            seen.add(tuple(sum(c * im[t] for c, im in zip(co, imgs)) % i2[t] for t in range(len(i2))))
        return len(seen) == math.prod(i2)

    def bt(k, imgs):
        if k == m:
            return surj(imgs)
        for y in tgt[k]:
            if all(_b(y, imgs[j], B2) == B1[k][j] for j in range(k)):
                if bt(k + 1, imgs + [y]):
                    return True
        return False
    return bt(0, [])
