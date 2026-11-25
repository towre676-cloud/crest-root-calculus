"""Refined Humbert invariants of level-l gluings E x E -> J for CM curves E = C/O_K, exactly.

K = Q(sqrt d), d a fundamental discriminant of class number 1, O_K = Z[w], w = (d + sqrt d)/2.
Elements of K are pairs (p, q) = p + q w.  E_h(x, y) = w-coefficient of x * conj(y) is the
principal polarisation of C/O_K (unimodular on O_K).  A gluing is psi in GL_2(F_l) acting on
O_K / l = F_l^2 in the basis (1, w); Lambda_J = O_K^2 + Z (x, psi x)/l, theta = l (E_h + E_h) restricted to Lambda_J (pullback of theta_J is l theta_0).
psi is admissible iff theta is integral and unimodular on Lambda_J (i.e. det psi = -1).

NS(J) = { H = [[a, b], [conj b, c]] : a, c in Q, b in K, E_H integral on Lambda_J },
E_H(x, y) = w-coefficient of x^T H conj(y).  D^2 = 2 Pf(E_D) (orientation fixed by Pf(E_theta) = 1),
D.D' = Pf(E_D + E_D') - Pf(E_D) - Pf(E_D').  Kani's refined Humbert invariant is
q(D) = (D.theta)^2 - 2 D^2 on NS / Z theta: positive definite ternary (rho = 4), and degree-n
elliptic subcovers correspond to primitive representations of n^2 by q (Kani).
"""
from __future__ import annotations

import itertools
from math import gcd
from fractions import Fraction as Fr

import cypari2

pari = cypari2.Pari()


class K:
    def __init__(self, d):
        self.d = d
        self.c = Fr(d * d - d, 4)          # w^2 = d w - c

    def mul(self, x, y):
        p1, q1 = x
        p2, q2 = y
        qq = q1 * q2
        return (p1 * p2 - qq * self.c, p1 * q2 + q1 * p2 + qq * self.d)

    def conj(self, x):
        p, q = x
        return (p + q * self.d, -q)

    def add(self, *xs):
        return (sum(x[0] for x in xs), sum(x[1] for x in xs))


def lattice_J(Kf, l, psi):
    """Z-basis (4 vectors in Q^4, coordinates (p1, q1, p2, q2)) of Lambda_J."""
    gens = [[1, 0, 0, 0], [0, 1, 0, 0], [0, 0, 1, 0], [0, 0, 0, 1]]
    for col in range(2):
        x = [1, 0] if col == 0 else [0, 1]
        y = [psi[0][col], psi[1][col]]
        gens.append([x[0], x[1], y[0], y[1]])   # divided by l below
    M = pari.matrix(4, 6, [(l if j < 4 else 1) * gens[j][i] for i in range(4) for j in range(6)])
    H = pari.mathnf(M)
    return [[Fr(int(H[i, j]), l) for i in range(4)] for j in range(4)]


def E_of(Kf, v, x, y):
    a, c, bp, bq = v
    x1, x2 = (x[0], x[1]), (x[2], x[3])
    y1, y2 = (y[0], y[1]), (y[2], y[3])
    b = (bp, bq)
    t = Kf.add(Kf.mul(Kf.mul(x1, (a, 0)), Kf.conj(y1)), Kf.mul(Kf.mul(x1, b), Kf.conj(y2)),
               Kf.mul(Kf.mul(x2, Kf.conj(b)), Kf.conj(y1)), Kf.mul(Kf.mul(x2, (c, 0)), Kf.conj(y2)))
    return t[1]


def pf(Kf, B, v):
    e = lambda i, j: E_of(Kf, v, B[i], B[j])
    return e(0, 1) * e(2, 3) - e(0, 2) * e(1, 3) + e(0, 3) * e(1, 2)


def ns_and_q(d, l, psi):
    Kf = K(d)
    B = lattice_J(Kf, l, psi)
    theta = (Fr(l), Fr(l), Fr(0), Fr(0))
    Et = [[E_of(Kf, theta, B[i], B[j]) for j in range(4)] for i in range(4)]
    if any(x.denominator != 1 for r in Et for x in r):
        return None
    s = pf(Kf, B, theta)
    if abs(s) != 1:
        return None
    # functionals v -> E_v(B_i, B_j), i < j, rows of a 6 x 4 rational matrix (linear in v)
    basisv = [(Fr(1), Fr(0), Fr(0), Fr(0)), (Fr(0), Fr(1), Fr(0), Fr(0)),
              (Fr(0), Fr(0), Fr(1), Fr(0)), (Fr(0), Fr(0), Fr(0), Fr(1))]
    rows = [[E_of(Kf, bv, B[i], B[j]) for bv in basisv] for i, j in itertools.combinations(range(4), 2)]
    den = 1
    for r in rows:
        for x in r:
            den = den * x.denominator // gcd(den, x.denominator)
    R = pari.matrix(4, 6, [int(rows[j][i] * den) for i in range(4) for j in range(6)])
    Hr = pari.mathnf(R)                          # 4 x 4, columns span (den * row lattice)
    # dual basis: v such that r . v in Z  <=>  v in columns of (Hr^T / den)^{-1}
    Dm = pari.matsolve(pari.mattranspose(Hr) / den, pari.matid(4))
    ns = [[Fr(str(Dm[i, j])) for i in range(4)] for j in range(4)]
    P = lambda v: pf(Kf, B, v) * s
    def dot(u, v):
        w = tuple(a + b for a, b in zip(u, v))
        return P(w) - P(u) - P(v)
    Gns = [[dot(u, v) if i != j else 2 * P(u) for j, v in enumerate(ns)] for i, u in enumerate(ns)]
    assert all(x.denominator == 1 for r in Gns for x in r)
    thm = pari.matsolve(pari.matrix(4, 4, [str(ns[j][i]) for i in range(4) for j in range(4)]),
                        pari.matrix(4, 1, [str(x) for x in theta]))
    th = [Fr(str(thm[i, 0])) for i in range(4)]
    assert all(x.denominator == 1 for x in th)
    G = [[int(x) for x in r] for r in Gns]
    t = [int(x) for x in th]
    return {"NS_gram": G, "theta": t, "q": q_from_ns(G, t)}


def q_from_ns(G, t):
    """Kani's refined Humbert form on NS / Z theta from an NS Gram matrix and theta."""
    tdot = [sum(G[i][j] * t[j] for j in range(4)) for i in range(4)]
    assert sum(t[i] * tdot[i] for i in range(4)) == 2
    Bq = [[tdot[i] * tdot[j] - 2 * G[i][j] for j in range(4)] for i in range(4)]
    # Extend theta to an integral basis before taking the quotient.
    T = pari.matsnf(pari.matrix(4, 1, t), 1)
    Ux = pari.matsolve(T[0], pari.matid(4))       # columns: first is +-t / gcd
    cols = [[int(Ux[i, j]) for i in range(4)] for j in (3, 0, 1, 2)]
    assert cols[0] in (t, [-x for x in t])
    q3 = [[sum(cols[a][i] * Bq[i][j] * cols[b][j] for i in range(4) for j in range(4)) for b in range(1, 4)]
          for a in range(1, 4)]
    first = [sum(cols[0][i] * Bq[i][j] * cols[0][j] for i in range(4) for j in range(4))]
    assert first == [0]
    return q3
