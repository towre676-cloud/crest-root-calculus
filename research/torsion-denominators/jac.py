"""Exact genus-two Jacobian arithmetic over F_{p^r} and integral torsion structure.

Curves are given as y^2 = f(x) with f squarefree of degree 5 or 6 over Z.  A degree-6
model is moved to an odd model by sending a root that is rational over F_p to infinity,
so the p-power Frobenius still acts coefficientwise on Mumford pairs.

Everything here is exact: field arithmetic is FLINT's F_q, group orders come from the
Weil polynomial, and the l-Sylow structure is computed from a complete Cayley graph
with a Smith normal form of its relation lattice.
"""
from __future__ import annotations

import random
from math import gcd

import flint


def _fp_roots(coeffs, p):
    """Roots in F_p of an integer polynomial (coefficients low to high)."""
    return [x for x in range(p) if sum(c * pow(x, i, p) for i, c in enumerate(coeffs)) % p == 0]


def odd_model(coeffs, p):
    """Return integer coefficients (mod p, low to high) of an odd-degree model over F_p.

    If deg f = 5 the model is f itself.  If deg f = 6 a root alpha in F_p is moved to
    infinity by x = alpha + 1/u, y = Y/u^3, giving Y^2 = u^6 f(alpha + 1/u).
    """
    c = [x % p for x in coeffs]
    while c and c[-1] == 0:
        c.pop()
    if len(c) == 6:
        return c, None
    if len(c) != 7:
        raise ValueError("degree must be 5 or 6 modulo p")
    roots = _fp_roots(c, p)
    if not roots:
        raise ValueError("no F_p-rational root to move to infinity")
    a = roots[0]
    # g(u) = sum c_i (a u + 1)^i u^(6-i)
    g = [0] * 7
    for i, ci in enumerate(c):
        # (a u + 1)^i = sum_k binom(i,k) a^k u^k
        from math import comb
        for k in range(i + 1):
            g[k + 6 - i] = (g[k + 6 - i] + ci * comb(i, k) * pow(a, k, p)) % p
    assert g[6] == 0, "root did not move to infinity"
    return g[:6], a


class Jacobian:
    """Jacobian of Y^2 = g(u), deg g = 5, over F_q with q = p^r."""

    def __init__(self, gcoeffs, p, r):
        self.p, self.r = p, r
        self.F = flint.fq_default_ctx(p, r)
        self.R = flint.fq_default_poly_ctx(self.F)
        self.x = self.R.gen()
        self.f = sum((self.F(int(c)) * self.x ** i for i, c in enumerate(gcoeffs)), self.R(0))
        assert self.f.degree() == 5 and self.f.is_squarefree()
        self.zero = (self.R(1), self.R(0))

    # --- Mumford arithmetic (Cantor) -------------------------------------------------
    def _reduce(self, u, v):
        f = self.f
        while u.degree() > 2:
            u = (f - v * v).exact_division(u)
            u = u.monic()
            v = (-v) % u
        u = u.monic()
        return (u, v % u)

    def add(self, D1, D2):
        u1, v1 = D1
        u2, v2 = D2
        d1, e1, e2 = u1.xgcd(u2)
        if d1.is_one():
            s1, s2, s3 = e1, e2, self.R(0)
            d = d1
        else:
            d, c1, c2 = d1.xgcd(v1 + v2)
            s1, s2, s3 = c1 * e1, c1 * e2, c2
        u = (u1 * u2).exact_division(d * d)
        v = (s1 * u1 * v2 + s2 * u2 * v1 + s3 * (v1 * v2 + self.f)).exact_division(d) % u
        return self._reduce(u, v)

    def neg(self, D):
        return (D[0], (-D[1]) % D[0])

    def mul(self, n, D):
        if n < 0:
            return self.mul(-n, self.neg(D))
        R = self.zero
        while n:
            if n & 1:
                R = self.add(R, D)
            D = self.add(D, D)
            n >>= 1
        return R

    def frob(self, D, k=1):
        """p-power Frobenius applied k times, coefficientwise."""
        e = self.p ** (k % self.r) if self.r else 1
        u, v = D
        return (self.R([c ** e for c in u.coeffs()]), self.R([c ** e for c in v.coeffs()]))

    def key(self, D):
        return (str(D[0]), str(D[1]))

    def is_zero(self, D):
        return D[0].is_one()

    # --- random points --------------------------------------------------------------
    def _rand_elem(self, rng):
        return self.F([rng.randrange(self.p) for _ in range(self.r)]) if self.r > 1 else self.F(rng.randrange(self.p))

    def _sqrt_mod(self, w, u, rng):
        """Square root of w in F_q[x]/(u), u monic irreducible of degree 2 (Tonelli-Shanks)."""
        Q = (self.p ** self.r) ** 2
        one = self.R(1)
        if w.pow_mod((Q - 1) // 2, u) != one:
            return None
        t, s = Q - 1, 0
        while t % 2 == 0:
            t //= 2
            s += 1
        while True:
            z = self.R([self._rand_elem(rng), self._rand_elem(rng)])
            if not z.is_zero() and z.pow_mod((Q - 1) // 2, u) != one:
                break
        M, c, T, Rr = s, z.pow_mod(t, u), w.pow_mod(t, u), w.pow_mod((t + 1) // 2, u)
        while T != one:
            i, TT = 0, T
            while TT != one:
                TT = TT.mul_mod(TT, u)
                i += 1
            bb = c
            for _ in range(M - i - 1):
                bb = bb.mul_mod(bb, u)
            M, c = i, bb.mul_mod(bb, u)
            T, Rr = T.mul_mod(c, u), Rr.mul_mod(bb, u)
        assert Rr.mul_mod(Rr, u) == w % u
        return Rr

    def random_divisor(self, rng):
        """Random reduced divisor: two rational points, or a conjugate pair over F_{q^2}."""
        if rng.random() < 0.5:
            return self.add(self.random_point_divisor(rng), self.random_point_divisor(rng))
        while True:
            u = self.R([self._rand_elem(rng), self._rand_elem(rng), 1])
            if not u.is_irreducible():
                continue
            v = self._sqrt_mod(self.f % u, u, rng)
            if v is not None:
                return (u, v)

    def random_point_divisor(self, rng):
        F = self.F
        while True:
            x0 = self._rand_elem(rng)
            y2 = self.f(x0)
            if y2.is_zero():
                continue
            if y2.is_square():
                y0 = y2.sqrt()
                return (self.x - x0, self.R(y0))


# --- Weil polynomials -------------------------------------------------------------------
def count_points(gcoeffs, p, r):
    """#C(F_{p^r}) for the odd model (one point at infinity)."""
    F = flint.fq_default_ctx(p, r)
    g = [F(int(c)) for c in gcoeffs]
    q = p ** r
    n = 1
    if r == 1:
        elems = [F(i) for i in range(p)]
    else:
        elems = _all_elements(F, p, r)
    for x in elems:
        val = F(0)
        for c in reversed(g):
            val = val * x + c
        if val.is_zero():
            n += 1
        elif val.is_square():
            n += 2
    return n


def _all_elements(F, p, r):
    import itertools
    gen = F.gen()
    basis = [gen ** i for i in range(r)]
    for t in itertools.product(range(p), repeat=r):
        s = F(0)
        for c, b in zip(t, basis):
            if c:
                s += c * b
        yield s


def weil_poly(gcoeffs, p):
    """Integer coefficients [1, -s1, s2, -p s1, p^2] of the Weil polynomial (high to low)."""
    N1 = count_points(gcoeffs, p, 1)
    N2 = count_points(gcoeffs, p, 2)
    t1 = p + 1 - N1            # sum of alpha
    t2 = p * p + 1 - N2        # sum of alpha^2
    s1 = t1
    s2 = (t1 * t1 - t2) // 2
    return [1, -s1, s2, -p * s1, p * p]


def power_sums(P, p, K):
    """Power sums p_0..p_K of the Frobenius roots, exact integers (Newton)."""
    s1, s2 = -P[1], P[2]
    e = [1, s1, s2, p * s1, p * p]
    pw = [4]
    for k in range(1, K + 1):
        val = sum((-1) ** (i - 1) * e[i] * pw[k - i] for i in range(1, min(k, 4) + 1) if k - i > 0)
        if k <= 4:
            val += (-1) ** (k - 1) * k * e[k]
        pw.append(val)
    return pw


def jac_order(P, p, r):
    """#J(F_{p^r}) = prod (1 - alpha_i^r), from power sums of beta = alpha^r."""
    pw = power_sums(P, p, 4 * r)
    b = [pw[r * k] for k in range(5)]
    E = [1]
    for k in range(1, 5):
        E.append(sum((-1) ** (i - 1) * E[k - i] * b[i] for i in range(1, k + 1)) // k)
    return 1 - E[1] + E[2] - E[3] + E[4]


def real_weil_roots(P, p):
    """Integer roots of h(s) = s^2 - s1 s + (s2 - 2p) if they exist, else None."""
    s1, s2 = -P[1], P[2]
    disc = s1 * s1 - 4 * (s2 - 2 * p)
    from math import isqrt
    if disc < 0 or isqrt(disc) ** 2 != disc:
        return None
    rt = isqrt(disc)
    return ((s1 + rt) // 2, (s1 - rt) // 2)


def frob_trace_power(a, p, r):
    """Trace of pi^r on an elliptic block with Frobenius trace a."""
    t0, t1 = 2, a
    if r == 0:
        return 2
    for _ in range(r - 1):
        t0, t1 = t1, a * t1 - p * t0
    return t1


def geometrically_distinct(a1, a2, p):
    """Blocks are non-isogenous over the algebraic closure iff traces of pi^120 differ.

    A ratio of Frobenius roots that is a root of unity lies in a field of degree at most 4,
    so its order n satisfies phi(n) <= 4 and divides 120.
    """
    return frob_trace_power(a1, p, 120) != frob_trace_power(a2, p, 120)


def vval(n, l):
    n = abs(n)
    if n == 0:
        return 10 ** 9
    k = 0
    while n % l == 0:
        n //= l
        k += 1
    return k


def prime_factors(n):
    n = abs(n)
    out, d = [], 2
    while d * d <= n:
        if n % d == 0:
            out.append(d)
            while n % d == 0:
                n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


# --- l-Sylow structure via Cayley graph and Smith normal form ---------------------------
def smith(A):
    """Smith normal form D = U A V for an integer matrix (lists). Returns (D, U, V)."""
    m, n = len(A), len(A[0])
    D = [row[:] for row in A]
    U = [[int(i == j) for j in range(m)] for i in range(m)]
    V = [[int(i == j) for j in range(n)] for i in range(n)]

    def swap_rows(M, i, j):
        M[i], M[j] = M[j], M[i]

    def swap_cols(M, i, j):
        for row in M:
            row[i], row[j] = row[j], row[i]

    t = 0
    while t < min(m, n):
        # find nonzero pivot of minimal absolute value in submatrix
        best = None
        for i in range(t, m):
            for j in range(t, n):
                if D[i][j] and (best is None or abs(D[i][j]) < abs(D[best[0]][best[1]])):
                    best = (i, j)
        if best is None:
            break
        i, j = best
        swap_rows(D, t, i); swap_rows(U, t, i)
        swap_cols(D, t, j); swap_cols(V, t, j)
        done = False
        while not done:
            done = True
            piv = D[t][t]
            for i in range(t + 1, m):
                q = D[i][t] // piv
                if q:
                    D[i] = [a - q * b for a, b in zip(D[i], D[t])]
                    U[i] = [a - q * b for a, b in zip(U[i], U[t])]
                if D[i][t]:
                    swap_rows(D, t, i); swap_rows(U, t, i)
                    done = False
                    break
            if not done:
                continue
            piv = D[t][t]
            for j in range(t + 1, n):
                q = D[t][j] // piv
                if q:
                    for row in D:
                        row[j] -= q * row[t]
                    for row in V:
                        row[j] -= q * row[t]
                if D[t][j]:
                    swap_cols(D, t, j); swap_cols(V, t, j)
                    done = False
                    break
            if not done:
                continue
            # divisibility condition
            piv = D[t][t]
            for i in range(t + 1, m):
                for j in range(t + 1, n):
                    if D[i][j] % piv:
                        D[t] = [a + b for a, b in zip(D[t], D[i])]
                        U[t] = [a + b for a, b in zip(U[t], U[i])]
                        done = False
                        break
                if not done:
                    break
        if D[t][t] < 0:
            D[t] = [-a for a in D[t]]
            U[t] = [-a for a in U[t]]
        t += 1
    return D, U, V


def _hnf_insert(basis, vec):
    """Keep a row-echelon integer basis of a lattice; insert vec."""
    vec = vec[:]
    for row in basis:
        piv = next(i for i, a in enumerate(row) if a)
        if vec[piv] == 0:
            continue
        # extended gcd combine
        a, b = row[piv], vec[piv]
        g, x, y = _egcd(a, b)
        new_row = [x * r + y * v for r, v in zip(row, vec)]
        vec = [(a // g) * v - (b // g) * r for r, v in zip(row, vec)]
        row[:] = new_row
    if any(vec):
        basis.append(vec)
        basis.sort(key=lambda r: next(i for i, a in enumerate(r) if a))


def _egcd(a, b):
    if b == 0:
        return (abs(a), 1 if a >= 0 else -1, 0)
    x0, y0, x1, y1 = 1, 0, 0, 1
    while b:
        q = a // b
        a, b = b, a - q * b
        x0, x1 = x1, x0 - q * x1
        y0, y1 = y1, y0 - q * y1
    if a < 0:
        a, x0, y0 = -a, -x0, -y0
    return a, x0, y0


class SylowStructure:
    """Complete l-Sylow subgroup S of J(F_q) with a Smith basis and coordinates."""

    def __init__(self, J, ell, order, seed=1, max_size=2 ** 20, ngens=7, tries=6):
        for t in range(tries):
            try:
                self._build(J, ell, order, seed + 1000 * t, max_size, ngens + t)
                return
            except _SpanError:
                continue
        raise RuntimeError("could not generate the Sylow subgroup")

    def _build(self, J, ell, order, seed, max_size, ngens):
        self.J, self.ell = J, ell
        v = vval(order, ell)
        self.size = ell ** v
        if self.size > max_size:
            raise RuntimeError(f"Sylow subgroup too large: {ell}^{v}")
        cof = order // self.size
        rng = random.Random(seed)
        gens = []
        while len(gens) < ngens:
            g = J.mul(cof, J.random_divisor(rng))
            gens.append(g)
        k = len(gens)
        coords = {J.key(J.zero): (J.zero, [0] * k)}
        frontier = [J.zero]
        lattice = []
        while frontier:
            nxt = []
            for D in frontier:
                c = coords[J.key(D)][1]
                for i, g in enumerate(gens):
                    E = J.add(D, g)
                    ce = c[:]
                    ce[i] += 1
                    key = J.key(E)
                    if key in coords:
                        rel = [a - b for a, b in zip(ce, coords[key][1])]
                        if any(rel):
                            _hnf_insert(lattice, rel)
                    else:
                        coords[key] = (E, ce)
                        nxt.append(E)
            frontier = nxt
        if len(coords) != self.size:
            raise _SpanError(f"generators span {len(coords)} of {self.size}")
        Dm, U, V = smith(lattice)
        self.invariants = [Dm[i][i] for i in range(k) if Dm[i][i] != 1]
        self.V = V
        self.k = k
        idx = [i for i in range(k) if Dm[i][i] != 1]
        self.idx = idx
        self.d = [Dm[i][i] for i in idx]
        self.coords = coords
        self.gens = gens
        # basis b_i = sum_j (V^{-1})_{ij} g_j for the nontrivial invariant positions
        Vinv = _inv_unimodular(V)
        self.basis = []
        for i in idx:
            B = J.zero
            for j in range(k):
                if Vinv[i][j]:
                    B = J.add(B, J.mul(Vinv[i][j], gens[j]))
            self.basis.append(B)

    def dlog(self, D):
        x = self.coords[self.J.key(D)][1]
        y = [sum(x[j] * self.V[j][i] for j in range(self.k)) for i in range(self.k)]
        return [y[i] % self.d[t] for t, i in enumerate(self.idx)]

    def torsion_rank(self, s):
        """Number of invariant factors divisible by ell^s."""
        return sum(1 for d in self.d if d % (self.ell ** s) == 0)

    def frobenius_matrix(self):
        """Columns are dlogs of F(b_i) in the Smith basis."""
        return [self.dlog(self.J.frob(b)) for b in self.basis]

    def torsion_generators(self, s):
        """Generators of S[ell^s] (requires S[ell^s] = J[ell^s], i.e. torsion_rank(s) == 4)."""
        return [self.J.mul(d // self.ell ** min(s, vval(d, self.ell)), b) for d, b in zip(self.d, self.basis)]


class _SpanError(Exception):
    pass


def _inv_unimodular(V):
    from fractions import Fraction
    n = len(V)
    M = [[Fraction(x) for x in row] + [Fraction(int(i == j)) for j in range(n)] for i, row in enumerate(V)]
    for c in range(n):
        piv = next(r for r in range(c, n) if M[r][c] != 0)
        M[c], M[piv] = M[piv], M[c]
        pv = M[c][c]
        M[c] = [a / pv for a in M[c]]
        for r in range(n):
            if r != c and M[r][c] != 0:
                fac = M[r][c]
                M[r] = [a - fac * b for a, b in zip(M[r], M[c])]
    out = [[int(M[i][n + j]) for j in range(n)] for i in range(n)]
    assert all(M[i][n + j].denominator == 1 for i in range(n) for j in range(n))
    return out
