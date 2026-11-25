"""Periods, homology intersection form, endomorphisms and elliptic subvarieties of a genus-two Jacobian.

Homology.  Fix a branch point e_1 and the star of straight segments from e_1 to the other
roots (a degree-5 model is first moved to degree 6).  The loop gamma_j lifts [e_1, e_j]
on the sheet y_+ and returns on y_-; its period is 2 * int_{e_1}^{e_j} omega_+.  In the
local coordinate z = sqrt(x - e_1) every gamma_j is a straight line through 0 oriented
along arg(e_j - e_1)/2, so with angles measured from a cut that misses every segment,
gamma_j . gamma_k = +1 whenever j precedes k.  The five loops span H_1 with the single
relation of alternating signs; any four of them form a basis.  Both the relation and
Riemann's bilinear relations are checked numerically as self-tests.

Branches.  On the segment x(t) = e_1 + (e_j - e_1) sin^2(t/2), t in [0, pi],
    dx / y = (e_j - e_1)^{1/2} / Q(t) dt,   Q(t)^2 = -(e_j - e_1) f(x) / ((x - e_1)(x - e_j)),
with Q smooth and nonvanishing.  Every factor sqrt(x - e_k) is taken relative to the nearest
point of the segment, where the principal branch is continuous; the global sign is fixed so
that y_+ ~ rho_0 sqrt(x - e_1) at e_1 for one rho_0 shared by all segments.
"""
from __future__ import annotations

import mpmath as mp

mp.mp.dps = 70


def poly_eval(c, x):
    s = mp.mpf(0)
    for a in reversed(c):
        s = s * x + a
    return s


def roots_of(c):
    return [mp.mpc(r) for r in mp.polyroots([mp.mpc(a) for a in reversed(c)], maxsteps=800, extraprec=800)]


Z0_CHOICES = [mp.mpc(a, 0) / 7 + mp.mpc(0, b) / 11 for a in (-5, -2, 1, 3, 6) for b in (2, 3, 5, 8, 13)]


def moved(f, z0):
    """Coefficients of x~^6 f(z0 + 1/x~): a sextic whose roots are 1/(r - z0)."""
    out = [mp.mpc(0)] * 7
    # f(z0 + u) expanded in u, then u = 1/x~ and multiply by x~^6
    from math import comb
    g = [mp.mpc(0)] * 7
    for i, c in enumerate(f):
        for k in range(i + 1):
            g[k] += c * comb(i, k) * z0 ** (i - k)
    for k in range(7):
        out[6 - k] = g[k]
    return out


def segment_clearance(e, i1):
    m = mp.inf
    for j in range(6):
        if j == i1:
            continue
        for k in range(6):
            if k in (i1, j):
                continue
            m = min(m, abs(nearest_point(e[i1], e[j], e[k]) - e[k]))
    return m


def nearest_point(a, b, p):
    d = b - a
    t = mp.re((p - a) * mp.conj(d)) / abs(d) ** 2
    t = min(max(t, 0), 1)
    return a + t * d


def star_periods(f):
    """Star periods of (dx/y, x dx/y) for a complex sextic f (coefficients low to high)."""
    e = roots_of(f)
    lc = mp.mpc(f[6])
    # choose e_1 and a cut direction in the largest angular gap
    best = None
    for i in range(6):
        angs = sorted(mp.arg(e[j] - e[i]) for j in range(6) if j != i)
        gaps = [(angs[(t + 1) % 5] - angs[t]) % (2 * mp.pi) for t in range(5)]
        g = max(gaps)
        cl = segment_clearance(e, i)
        if best is None or cl > best[0]:
            t = gaps.index(g)
            best = (cl, i, angs[t] + g / 2)
    _, i1, cut = best
    e1 = e[i1]
    clearance = segment_clearance(e, i1)
    others = [e[j] for j in range(6) if j != i1]

    def theta(z):
        th = mp.arg(z - e1)
        while th <= cut:
            th += 2 * mp.pi
        while th > cut + 2 * mp.pi:
            th -= 2 * mp.pi
        return th

    others.sort(key=theta)
    fprime_e1 = lc
    for z in others:
        fprime_e1 *= (e1 - z)
    rho0 = mp.sqrt(fprime_e1)
    periods = []
    for ej in others:
        th = theta(ej)
        sq = mp.sqrt(abs(ej - e1)) * mp.expj(th / 2)
        rest = [z for z in others if z != ej]
        mids = [nearest_point(e1, ej, z) for z in rest]
        pref = mp.sqrt(-(ej - e1) * lc)

        def Qraw(t):
            x = e1 + (ej - e1) * mp.sin(t / 2) ** 2
            q = pref
            for z, m in zip(rest, mids):
                q *= mp.sqrt(m - z) * mp.sqrt((x - z) / (m - z))
            return q

        # y_+ = (x-e1)^{1/2} R(x), R(e1) = rho0, and R = cos(t/2) Q
        sign = 1 if abs(Qraw(0) - rho0) < abs(Qraw(0) + rho0) else -1
        row = []
        for m in (0, 1):
            def integrand(t, m=m):
                x = e1 + (ej - e1) * mp.sin(t / 2) ** 2
                return x ** m * sq / (sign * Qraw(t))
            row.append(2 * mp.quad(integrand, mp.linspace(0, mp.pi, 9)))
        periods.append(row)
    # intersection matrix of the five loops in angular order
    I = [[0 if j == k else (1 if j < k else -1) for k in range(5)] for j in range(5)]
    relation = [periods[0][m] - periods[1][m] + periods[2][m] - periods[3][m] + periods[4][m] for m in range(2)]
    return {"roots": e, "base": e1, "periods": periods, "intersection": I, "clearance": clearance,
            "relation_residual": max(abs(r) for r in relation)}


def period_matrix(f):
    """Basis gamma_1..gamma_4 of H_1: 2x4 period matrix of (dx/y, x dx/y) for the ORIGINAL model."""
    best = None
    for z0 in Z0_CHOICES:
        e = roots_of(moved([mp.mpf(c) for c in f], z0))
        cl = max(segment_clearance(e, i) / max(abs(a - b) for a in e for b in e) for i in range(6))
        if best is None or cl > best[0]:
            best = (cl, z0)
    z0 = best[1]
    ft = moved([mp.mpf(c) for c in f], z0)
    S = star_periods(ft)
    # dx/y = -x~ dx~/y~,  x dx/y = -(z0 x~ + 1) dx~/y~
    P = []
    for w0, w1 in S["periods"]:
        P.append([-w1, -(z0 * w1 + w0)])
    S["z0"] = z0
    S["periods_original"] = P
    P = P[:4]
    Pi = mp.matrix(2, 4)
    for c in range(4):
        for r in range(2):
            Pi[r, c] = P[c][r]
    E = [row[:4] for row in S["intersection"][:4]]
    return Pi, E, S


def riemann_check(Pi, E):
    """First relation Pi E^{-1} Pi^T = 0 and positivity of i Pi E^{-1} conj(Pi)^T (up to a global sign)."""
    Em = mp.matrix(E)
    Ei = Em ** -1
    R1 = Pi * Ei * Pi.T
    H = mp.mpc(0, 1) * Pi * Ei * Pi.H
    first = max(abs(R1[i, j]) for i in range(2) for j in range(2))
    ev = mp.eighe(mp.matrix([[mp.re(H[i, j]) if i == j else H[i, j] for j in range(2)] for i in range(2)]))[0]
    return first, [mp.re(v) for v in ev], int(round(mp.det(Em)))
