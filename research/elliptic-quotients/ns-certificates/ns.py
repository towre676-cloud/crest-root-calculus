"""Exact Neron-Severi certificate and complete elliptic spectrum of a genus-two Jacobian.

Classes.  For a certified map f_i : C -> E_i of degree d_i put N_i = f_i^* f_{i*} in End(J) and
alpha_i = d_i - N_i.  Under D -> phi_theta^{-1} phi_D, the symmetric endomorphism alpha_i is a
class D_i in NS(J_Qbar) with D_i^2 = 2 det(alpha_i) = 0 and D_i . theta = tr(alpha_i) = d_i
(analytic traces); theta corresponds to 1.  For i != j,
    D_i . D_j = tr(alpha_i) tr(alpha_j) - tr(alpha_i alpha_j) = d_i d_j - deg(f_{i*} f_j^*),
and deg mu for mu = f_{i*} f_j^* : E_j -> E_i is read off mu^ mu = [deg mu] on one point of
large order.  None of this assumes that f_i is optimal.

Rank.  If some D_i . D_j > 0, J ~ E_i x E_j and rho(J) = 2 + rank Hom(E_i, E_j).  deg mu > 0
certifies isogeny; unequal Frobenius fields of E at two ordinary primes certify non-CM, so the
rank of Hom is 1; equal j in the 13-element rational CM list certifies CM, rank 2.  A Frobenius
witness (unequal traces of pi^120) certifies Hom = 0.

Saturation.  D = sum c_k D_k lies in l NS iff alpha_D = sum c_k alpha_k kills J[l] (l != p):
NS = Hom^sym(J, J^), theta is principal, and [l] is separable.  Every prime with l^2 | disc
is tested on the reduction, where J[l] is rational.

Spectrum.  By Kani, elliptic curves in J correspond to primitive D in NS with D^2 = 0 and
D . theta > 0, and the optimal map attached to D has degree D . theta.  On NS the form
Q(D) = (D . theta)^2 - D^2 is positive definite (Hodge index), so all D with D . theta <= X are
found by qfminim with bound X^2.
"""
from __future__ import annotations

import itertools
import random
from fractions import Fraction
from math import gcd

import cypari2

from setup import Reduction, good_primes
from jac import jac_order

pari = cypari2.Pari()
pari.allocatemem(4 * 10 ** 9)

RATIONAL_CM_J = {0, 1728, -3375, 8000, -32768, 54000, 287496, -884736, -12288000, 16581375,
                 -884736000, -147197952000, -262537412640768000}


# ------------------------------------------------------------------ intersection numbers
def deg_mu(R, i, j, trials=6):
    """deg(f_{i*} f_j^*) via mu^ mu = [deg mu] on points of E_j."""
    Ej = R.maps[j].E
    bound = R.maps[i].deg * R.maps[j].deg
    mod, val = 1, 0
    for P in R.aux_points(j, 60):
        if P is None:
            continue
        o = point_order(Ej, P, 10 ** 6)
        if o <= 2:
            continue
        g = R.push(j, R.pull(i, R.push(i, R.pull(j, P))))
        t = next((t for t in range(o) if Ej.eq(Ej.mul(t, P), g)), None)
        assert t is not None, "mu^ mu is not a multiple of P"
        # combine t mod o with previous residues
        from sympy.ntheory.modular import crt
        r = crt([mod, o], [val, t])
        if r is None:
            raise RuntimeError("inconsistent residues for deg mu")
        val, mod = int(r[0]), int(r[1])
        if mod > bound:
            break
    cands = [t for t in range(bound + 1) if t % mod == val % mod]
    if len(cands) != 1:
        raise RuntimeError(f"deg mu not determined: modulus {mod}")
    return cands[0]


def point_order(E, P, cap):
    Q, n = P, 1
    while Q is not None:
        Q = E.add(Q, P)
        n += 1
        if n > cap:
            return cap
    return n


def gram(R, idx):
    d = [R.maps[i].deg for i in idx]
    n = len(idx)
    G = [[0] * (n + 1) for _ in range(n + 1)]
    G[0][0] = 2
    for a in range(n):
        G[0][a + 1] = G[a + 1][0] = d[a]
    for a in range(n):
        for b in range(a + 1, n):
            m = deg_mu(R, idx[a], idx[b])
            G[a + 1][b + 1] = G[b + 1][a + 1] = d[a] * d[b] - m
    return G


# ------------------------------------------------------------------ rank bound
def trace_Fp(E, p):
    return p + 1 - (sum(1 for _ in E.points()) + 1)


def frobenius_field(a, p):
    disc = a * a - 4 * p
    return int(pari.coredisc(disc)) if disc != 0 else 0


# ------------------------------------------------------------------ saturation on J[l]
def l_torsion(R, ell, need=4, tries=4000):
    J = R.J
    n = R.order
    v = 0
    while n % ell == 0:
        n //= ell
        v += 1
    if v < need:
        return None
    cof = R.order // ell ** v
    group = {J.key(J.zero): J.zero}
    gens = []
    rng = random.Random(11)
    for _ in range(tries):
        Q = J.mul(cof, J.random_divisor(rng))
        if J.is_zero(Q):
            continue
        while not J.is_zero(J.mul(ell, Q)):
            Q = J.mul(ell, Q)
        if J.key(Q) in group:
            continue
        gens.append(Q)
        new = [Q]
        # closure: add multiples of Q to every element
        elems = list(group.values())
        for e in elems:
            S = e
            for _k in range(ell - 1):
                S = J.add(S, Q)
                if J.key(S) not in group:
                    group[J.key(S)] = S
        if len(group) == ell ** 4:
            return gens, group
    return None


def coords(J, gens, ell):
    """Map key -> coordinate vector over F_l for the group generated by independent gens."""
    table = {J.key(J.zero): (0,) * len(gens)}
    frontier = [(J.zero, (0,) * len(gens))]
    for k, g in enumerate(gens):
        new = []
        for e, c in list(table.items()):
            pass
        items = list(table.items())
        table2 = dict(table)
        # rebuild by adding multiples of g
        for key, c in items:
            pass
    # direct enumeration
    table = {}
    for c in itertools.product(range(ell), repeat=len(gens)):
        S = J.zero
        for ci, g in zip(c, gens):
            if ci:
                S = J.add(S, J.mul(ci, g))
        table[J.key(S)] = c
    return table


def saturate(R, basis_vecs, idx, primes_to_test):
    """basis_vecs: integer combinations (over the generators theta, D_idx) spanning N0.
    Returns the list of saturation steps found at each prime."""
    J = R.J
    steps = []
    for ell in primes_to_test:
        T = l_torsion(R, ell)
        if T is None:
            raise RuntimeError(f"J[{ell}] not rational over this field")
        gens, group = T
        table = coords(J, gens, ell)
        # matrices of generators' endomorphisms on J[l]: column k = coords of A(gen_k)
        def endo_matrix(fn):
            cols = [table[J.key(fn(g))] for g in gens]
            return [[cols[c][r] for c in range(4)] for r in range(4)]
        mats = [endo_matrix(lambda g: g)]
        for i in idx:
            di = R.maps[i].deg
            mats.append(endo_matrix(lambda g, i=i, di=di: J.add(J.mul(di, g), J.neg(R.pull(i, R.push(i, g))))))
        # vectors of the basis elements of N0 as endomorphisms mod l
        found = []
        for c in itertools.product(range(ell), repeat=len(basis_vecs)):
            if not any(c):
                continue
            # combination of basis vectors -> combination of generators
            coef = [sum(c[t] * basis_vecs[t][k] for t in range(len(basis_vecs))) % ell for k in range(len(mats))]
            A = [[sum(coef[k] * mats[k][r][s] for k in range(len(mats))) % ell for s in range(4)] for r in range(4)]
            if all(A[r][s] == 0 for r in range(4) for s in range(4)):
                found.append(c)
        steps.append({"l": ell, "divisible_combinations": found})
    return steps


# ------------------------------------------------------------------ spectrum from the lattice
def spectrum_from_lattice(G, theta_vec, X):
    """G: Gram matrix of a Z-basis of NS; theta_vec: coordinates of theta."""
    n = len(G)
    tv = [sum(G[r][c] * theta_vec[c] for c in range(n)) for r in range(n)]   # D.theta = tv . D
    Q = [[tv[r] * tv[c] - G[r][c] for c in range(n)] for r in range(n)]
    M = pari.matrix(n, n, [x for row in Q for x in row])
    assert all(e > 0 for e in pari.qfjacobi(M)[0]), "Q is not positive definite"
    res = pari.qfminim(M, X * X, 10 ** 8)
    vecs = res[2]
    from collections import Counter
    spec = Counter()
    classes = []
    for j in range(len(vecs)):
        D = [int(vecs[r, j]) for r in range(n)]
        g = 0
        for x in D:
            g = gcd(g, x)
        if g != 1:
            continue
        sq = sum(D[r] * G[r][c] * D[c] for r in range(n) for c in range(n))
        if sq != 0:
            continue
        t = sum(tv[r] * D[r] for r in range(n))
        if t < 0:
            D, t = [-x for x in D], -t
        if 0 < t <= X:
            spec[t] += 1
            classes.append((t, D))
    return dict(sorted(spec.items())), classes


# ------------------------------------------------------------------ lattice bookkeeping
def lattice_basis(G):
    """Z-basis (lists in generator coordinates) of Z^n / ker(G), the image of the generators in NS."""
    n = len(G)
    M = pari.matrix(n, n, [x for r in G for x in r])
    K = pari.matkerint(M)
    k = len(K)
    if k == 0:
        return [[int(i == j) for i in range(n)] for j in range(n)]
    U, V, D = pari.matsnf(K, 1)[:3]
    Y = U * K
    support = {r for r in range(n) if any(Y[r, c] != 0 for c in range(k))}
    assert len(support) == k
    nz = [D[r, c] for r in range(n) for c in range(k) if D[r, c] != 0]
    assert len(nz) == k and all(abs(int(x)) == 1 for x in nz), "kernel not saturated"
    support = {r for r in range(n) if any(D[r, c] != 0 for c in range(k))}
    Ui = U ** -1
    return [[int(Ui[r, c]) for r in range(n)] for c in range(n) if c not in support]


def gram_of(G, B):
    n = len(G)
    return [[sum(B[a][r] * G[r][s] * B[b][s] for r in range(n) for s in range(n)) for b in range(len(B))]
            for a in range(len(B))]


def find_torsion_field(f, maps, ell, need_i, rmax=4, pstart=11, pstop=600):
    """A good (p, r) with p != l over which J[l] is rational; returns (Reduction, torsion data)."""
    for p in good_primes(f, maps, need_i, start=pstart, stop=pstop):
        if p == ell:
            continue
        for r in range(1, rmax + 1):
            try:
                R = Reduction(f, maps, p, r)
            except AssertionError:
                break
            n = R.order
            v = 0
            while n % ell == 0:
                n //= ell
                v += 1
            if v < 4 or ell ** v > 10 ** 7:
                continue
            T = l_torsion(R, ell)
            if T is not None:
                return R, T
    raise RuntimeError(f"no field with rational J[{ell}] found")


def saturate_lattice(f, maps, G, idx, need_i, log):
    """Saturate the image lattice of the generators (theta, D_idx) inside NS, exactly."""
    from fractions import Fraction
    B = [[Fraction(x) for x in v] for v in lattice_basis(G)]
    n = len(G)
    while True:
        GB = gram_of(G, B)
        det = int(pari.matdet(pari.matrix(len(B), len(B), [int(x) for r in GB for x in r])))
        primes = [int(q) for q in pari.factor(abs(det))[0]] if abs(det) > 1 else []
        cand = [q for q in primes if (det % (q * q) == 0)]
        changed = False
        for ell in cand:
            R, T = find_torsion_field(f, maps, ell, need_i)
            gens, group = T
            J = R.J
            table = coords(J, gens, ell)

            def endo_matrix(fn):
                cols = [table[J.key(fn(g))] for g in gens]
                return [[cols[c][r] for c in range(4)] for r in range(4)]
            mats = [endo_matrix(lambda g: g)]
            for i in idx:
                di = R.maps[i].deg
                mats.append(endo_matrix(lambda g, i=i, di=di: J.add(J.mul(di, g), J.neg(R.pull(i, R.push(i, g))))))
            # integer basis vectors (scale fractions: they are integral combos after saturation steps)
            den = 1
            for v in B:
                for x in v:
                    den = den * x.denominator // gcd(den, x.denominator)
            hit = None
            for c in itertools.product(range(ell), repeat=len(B)):
                if not any(c):
                    continue
                coef = [sum(c[t] * B[t][k] for t in range(len(B))) for k in range(n)]
                if any(Fraction(x).denominator % ell == 0 for x in coef):
                    continue
                coefm = [int(Fraction(x).numerator * pow(Fraction(x).denominator, -1, ell)) % ell for x in coef]
                A = [[sum(coefm[k] * mats[k][r][s] for k in range(n)) % ell for s in range(4)] for r in range(4)]
                if all(A[r][s] == 0 for r in range(4) for s in range(4)):
                    hit = c
                    break
            log.append({"l": ell, "field": [R.K.p, R.K.r], "divisible": hit})
            if hit is not None:
                new = [sum(Fraction(hit[t]) * B[t][k] for t in range(len(B))) / ell for k in range(n)]
                # replace one basis vector with nonzero coefficient mod l by the new vector
                pos = next(t for t in range(len(B)) if hit[t] % ell)
                B[pos] = new
                changed = True
                break
        if not changed:
            return B, gram_of(G, B)


def gram_of(G, B):
    from fractions import Fraction
    n = len(G)
    out = [[sum(Fraction(B[a][r]) * G[r][s] * Fraction(B[b][s]) for r in range(n) for s in range(n))
             for b in range(len(B))] for a in range(len(B))]
    for row in out:
        for x in row:
            assert Fraction(x).denominator == 1, "non-integral Gram after saturation step"
    return [[int(x) for x in row] for row in out]
