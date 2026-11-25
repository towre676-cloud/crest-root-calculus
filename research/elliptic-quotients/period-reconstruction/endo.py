"""Endomorphisms, elliptic subvarieties and optimal degrees from the period matrix.

Let Pr be the real 4x4 matrix whose columns are the basis cycles in R^4 = C^2 and J0 the
standard complex structure.  Jl = Pr^{-1} J0 Pr is the complex structure in lattice
coordinates.  An integer matrix R is an endomorphism iff Jl R = R Jl.  The integer kernel of
this linear condition is found by LLL at 70 digits; every returned generator is re-checked
to 10^-50.  A rank-one complex endomorphism is exactly an integer R != 0 with det R = 0, its
saturated image L is the lattice of an elliptic subvariety, and the optimal degree of the
corresponding elliptic quotient of the curve is |E(u1, u2)| for a basis of L (polarized-
reduction identity: degree of the restricted principal polarization).
"""
from __future__ import annotations

import itertools
from fractions import Fraction

import cypari2
import mpmath as mp

from periods import period_matrix, riemann_check

pari = cypari2.Pari()
pari.allocatemem(10 ** 9)


def real_frame(Pi):
    Pr = mp.matrix(4, 4)
    for c in range(4):
        Pr[0, c], Pr[1, c] = mp.re(Pi[0, c]), mp.im(Pi[0, c])
        Pr[2, c], Pr[3, c] = mp.re(Pi[1, c]), mp.im(Pi[1, c])
    J0 = mp.matrix([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
    return Pr, Pr ** -1 * J0 * Pr


def integer_kernel(cols, scale=mp.mpf(10) ** 55, tol=mp.mpf(10) ** -45):
    """Integer vectors c with sum c_i cols[i] = 0 (cols are real vectors)."""
    n, m = len(cols), len(cols[0])
    rows = []
    for r in range(m):
        rows.append([int(mp.nint(scale * cols[i][r])) for i in range(n)])
    for r in range(n):
        rows.append([1 if i == r else 0 for i in range(n)])
    M = pari.matrix(len(rows), n, [x for row in rows for x in row])
    U = pari.qflll(M)
    B = M * U
    ker = []
    for j in range(n):
        top = [B[r, j] for r in range(m)]
        vec = [int(U[i, j]) for i in range(n)]
        if not any(vec):
            continue
        resid = max(abs(sum(vec[i] * cols[i][r] for i in range(n))) for r in range(m))
        if resid < tol:
            ker.append(vec)
    return ker


def endomorphisms(Jl):
    basis = []
    for a in range(4):
        for b in range(4):
            Eab = mp.zeros(4, 4)
            Eab[a, b] = 1
            C = Jl * Eab - Eab * Jl
            basis.append([C[i, j] for i in range(4) for j in range(4)])
    ker = integer_kernel(basis)
    mats = [[[v[4 * a + b] for b in range(4)] for a in range(4)] for v in ker]
    return saturate_basis(mats)


def saturate_basis(mats):
    """Z-basis of the saturation of the span (HNF over Z^16 then saturation by matkerint)."""
    if not mats:
        return []
    vecs = [[m[a][b] for a in range(4) for b in range(4)] for m in mats]
    M = pari.matrix(16, len(vecs), [vecs[j][i] for i in range(16) for j in range(len(vecs))])
    # saturation: kernel of the kernel of the transpose
    K = pari.matkerint(pari.mattranspose(M))
    S = pari.matkerint(pari.mattranspose(K)) if len(K) else pari.matid(16)
    out = []
    for j in range(len(S)):
        v = [int(S[i, j]) for i in range(16)]
        out.append([[v[4 * a + b] for b in range(4)] for a in range(4)])
    return out


def neron_severi(Jl):
    """Integral alternating forms E with Jl^T E Jl = E."""
    pairs = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]
    cols = []
    for (a, b) in pairs:
        Eab = mp.zeros(4, 4)
        Eab[a, b], Eab[b, a] = 1, -1
        C = Jl.T * Eab * Jl - Eab
        cols.append([C[i, j] for i in range(4) for j in range(4)])
    ker = integer_kernel(cols)
    return len(ker)


def matmul(A, B):
    return [[sum(A[i][k] * B[k][j] for k in range(4)) for j in range(4)] for i in range(4)]


def det4(A):
    return int(pari.matdet(pari.matrix(4, 4, [x for row in A for x in row])))


def saturated_image(R):
    M = pari.matrix(4, 4, [x for row in R for x in row])
    K = pari.matkerint(pari.mattranspose(M))       # orthogonal complement (rows)
    L = pari.matkerint(pari.mattranspose(K))       # saturated column space
    L = pari.mathnf(L)
    return tuple(tuple(int(L[i, j]) for j in range(len(L))) for i in range(4))


def form(E, u, v):
    return sum(u[i] * E[i][j] * v[j] for i in range(4) for j in range(4))


def elliptic_subvarieties(f, box=4):
    Pi, E, S = period_matrix(f)
    first, ev, det = riemann_check(Pi, E)
    if max(ev) < 0:
        E = [[-x for x in row] for row in E]
    Pr, Jl = real_frame(Pi)
    End = endomorphisms(Jl)
    rho = neron_severi(Jl)
    found = {}
    rng = range(-box, box + 1)
    for c in itertools.product(rng, repeat=len(End)):
        if not any(c):
            continue
        R = [[sum(c[t] * End[t][a][b] for t in range(len(End))) for b in range(4)] for a in range(4)]
        if det4(R) != 0:
            continue
        L = saturated_image(R)
        if len(L[0]) != 2 or L in found:
            continue
        u = [L[i][0] for i in range(4)]
        v = [L[i][1] for i in range(4)]
        found[L] = abs(form(E, u, v))
    return {"Pi": Pi, "E": E, "End": End, "rank_End": len(End), "rank_NS": rho, "riemann_residual": first,
            "relation_residual": S["relation_residual"], "subvarieties": found, "S": S}


# ------------------------------------------------------------------ complete enumeration
PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def hermitian_gram(E, Jl):
    """Real Gram matrix of H(u, u) = E(Jl u, u) (sign fixed to be positive definite)."""
    Em = mp.matrix(E)
    G = Jl.T * Em.T
    G = (G + G.T) / 2
    if mp.det(G) < 0 or G[0, 0] < 0:
        G = -G
    return G


def plucker_integral(u, Ju, bound=10 ** 9, tol=1e-30):
    p = [u[i] * Ju[j] - u[j] * Ju[i] for i, j in PAIRS]
    piv = max(range(6), key=lambda t: abs(p[t]))
    if abs(p[piv]) < tol:
        return None
    r = [x / p[piv] for x in p]
    fr = []
    for x in r:
        q = pari.bestappr(pari(mp.nstr(x, 60)), bound)
        if abs(mp.mpf(str(pari.numerator(q))) / mp.mpf(str(pari.denominator(q))) - x) > tol:
            return None
        fr.append(Fraction(int(pari.numerator(q)), int(pari.denominator(q))))
    from math import lcm, gcd
    den = 1
    for x in fr:
        den = lcm(den, x.denominator)
    v = [int(x * den) for x in fr]
    g = 0
    for x in v:
        g = gcd(g, x)
    v = [x // g for x in v]
    s = next(x for x in v if x)
    return tuple(x if s > 0 else -x for x in v)


def complete_spectrum(f, X):
    """Propose elliptic directions using a numerical form and the Hermite bound."""
    Pi, E, S = period_matrix(f)
    first, ev, det = riemann_check(Pi, E)
    if max(ev) < 0:
        E = [[-x for x in row] for row in E]
    Pr, Jl = real_frame(Pi)
    G = hermitian_gram(E, Jl)
    bound = 2 * X / mp.sqrt(3) * (1 + mp.mpf(10) ** -20)
    Gp = pari.matrix(4, 4, [pari(mp.nstr(G[i, j], 65)) for i in range(4) for j in range(4)])
    res = pari.qfminim(Gp, pari(mp.nstr(bound, 30)), 10 ** 8, 2)
    vecs = res[2]
    import numpy as np
    Jf = np.array([[float(Jl[i, j]) for j in range(4)] for i in range(4)])
    V = np.array([[int(vecs[i, j]) for j in range(len(vecs))] for i in range(4)], dtype=float)
    JV = Jf @ V
    found = {}
    for jcol in range(V.shape[1]):
        u, ju = V[:, jcol], JV[:, jcol]
        p = np.array([u[i] * ju[j] - u[j] * ju[i] for i, j in PAIRS])
        piv = np.argmax(np.abs(p))
        r = p / p[piv]
        ok = True
        for x in r:
            fx = Fraction(float(x)).limit_denominator(10 ** 5)
            if abs(float(fx) - x) > 1e-9:
                ok = False
                break
        if not ok:
            continue
        u = [int(vecs[i, jcol]) for i in range(4)]
        Ju = Jl * mp.matrix(u)
        pl = plucker_integral(u, [Ju[i] for i in range(4)])
        if pl is None or pl in found:
            continue
        # Plucker relation p01 p23 - p02 p13 + p03 p12 = 0 for a decomposable bivector
        assert pl[0] * pl[5] - pl[1] * pl[4] + pl[2] * pl[3] == 0
        d = abs(sum(E[i][j] * c for (i, j), c in zip(PAIRS, pl)))
        if d <= X:
            found[pl] = d
    return {"E": E, "Jl": Jl, "Pi": Pi, "S": S, "count_vectors": len(vecs), "subvarieties": found,
            "riemann_residual": first, "relation_residual": S["relation_residual"]}
