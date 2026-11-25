"""Explicit optimal elliptic maps from a lattice direction, verified exactly.

Given a degree-d elliptic subvariety with saturated lattice L (Plucker vector pi), the
complementary subvariety has lattice L^perp = {v : E(v, L) = 0}.  The optimal quotient is
q : J -> J / B^c.  Analytically it is z -> l(z), with l = (1, lam) vanishing on C Pi(L^perp).
The target lattice is l(Pi Z^4) = Z l(Pi m1) + Z l(Pi m2) for a complement {m1, m2} of
L^perp.  With a Weierstrass base point r0 the map is P -> l(AJ_{r0}(P)), so its abscissa
X = wp(l(AJ(P))) is a rational function R(x) of degree d.  R is interpolated from sample
points, its coefficients and lam, g2, g3 are recognised as algebraic numbers, and the
result is verified exactly by substituting the horizontal and vertical coordinates.
"""
from __future__ import annotations

import json
import sys

from certificate import certify_record

import cypari2
import mpmath as mp

from endo import complete_spectrum, PAIRS
from periods import roots_of, nearest_point

pari = cypari2.Pari()
pari.allocatemem(10 ** 9)
# GP string conversion uses its display precision, independently of ellwp precision.
pari.set_real_precision_bits(320)
mp.mp.dps = 70


def lattice_from_plucker(pl):
    """Integer basis of {v in Z^4 : v wedge pi = 0}: solve the linear system v_i p_jk - ... = 0."""
    # v wedge pi = 0 gives 4 linear equations (one per triple i<j<k)
    p = dict(zip(PAIRS, pl))
    P = lambda i, j: p[(i, j)] if i < j else (-p[(j, i)] if i > j else 0)
    rows = []
    for (i, j, k) in [(0, 1, 2), (0, 1, 3), (0, 2, 3), (1, 2, 3)]:
        row = [0, 0, 0, 0]
        row[i] += P(j, k)
        row[j] -= P(i, k)
        row[k] += P(i, j)
        rows.append(row)
    M = pari.matrix(4, 4, [x for r in rows for x in r])
    K = pari.matkerint(M)
    return [[int(K[i, j]) for i in range(4)] for j in range(len(K))]


def orth_complement(E, L):
    rows = [[sum(E[i][j] * u[j] for j in range(4)) for i in range(4)] for u in L]
    M = pari.matrix(len(rows), 4, [x for r in rows for x in r])
    K = pari.matkerint(M)
    return [[int(K[i, j]) for i in range(4)] for j in range(len(K))]


def unimodular_complement(B):
    """Vectors m1, m2 with det[B0, B1, m1, m2] = +-1 (B saturated of rank 2)."""
    import itertools
    cand = []
    for c in itertools.product(range(-1, 2), repeat=4):
        if any(c):
            cand.append(list(c))
    for m1, m2 in itertools.combinations(cand, 2):
        M = pari.matrix(4, 4, [v[i] for i in range(4) for v in (B[0], B[1], m1, m2)])
        if abs(int(pari.matdet(M))) == 1:
            return m1, m2
    raise RuntimeError("no small unimodular complement")


def lattice_invariants(w1, w2):
    """g2, g3 of the lattice Z w1 + Z w2 via reduction of tau and Eisenstein q-series."""
    a, b, c, d = 1, 0, 0, 1   # (w1, w2) -> (a w1 + b w2, c w1 + d w2)
    W1, W2 = w1, w2
    for _ in range(200):
        tau = W1 / W2
        if mp.im(tau) < 0:
            W1 = -W1
            continue
        n = mp.nint(mp.re(tau))
        if n:
            W1 = W1 - n * W2
            continue
        if abs(W1) < abs(W2):
            W1, W2 = W2, -W1
            continue
        break
    tau = W1 / W2
    q = mp.exp(2j * mp.pi * tau)
    E4, E6 = mp.mpf(1), mp.mpf(1)
    qn = mp.mpf(1)
    for n in range(1, 400):
        qn *= q
        s3 = sum(dd ** 3 for dd in range(1, n + 1) if n % dd == 0)
        s5 = sum(dd ** 5 for dd in range(1, n + 1) if n % dd == 0)
        E4 += 240 * s3 * qn
        E6 -= 504 * s5 * qn
        if abs(qn) * n ** 6 < mp.mpf(10) ** -75:
            break
    k = 2 * mp.pi / W2
    return k ** 4 * E4 / 12, k ** 6 * E6 / 216


def aj(f, r0, xs, roots):
    """int_{r0}^{xs} (dx/y, x dx/y) along a straight path, continuous branch, x = r0 + (xs-r0) s^2."""
    lc = mp.mpf(f[6]) if len(f) == 7 and f[6] else mp.mpf(f[5])
    rest = [z for z in roots if abs(z - r0) > mp.mpf(10) ** -30]
    mids = [nearest_point(r0, xs, z) for z in rest]
    sq = mp.sqrt(xs - r0)

    def g_sqrt(x):
        q = mp.sqrt(lc)
        for z, m in zip(rest, mids):
            q *= mp.sqrt(m - z) * mp.sqrt((x - z) / (m - z))
        return q

    out = []
    for m in (0, 1):
        def integrand(s, m=m):
            x = r0 + (xs - r0) * s * s
            return x ** m * 2 * sq / g_sqrt(x)
        out.append(mp.quad(integrand, mp.linspace(0, 1, 13), method="gauss-legendre"))
    return out


def recognise(z, deg=2, bound=10 ** 30):
    """Minimal polynomial of a complex number (real part and imaginary part checked)."""
    zp = pari(f"{mp.nstr(mp.re(z), 65)} + I*({mp.nstr(mp.im(z), 65)})")
    for d in range(1, deg + 1):
        P = pari.algdep(zp, d)
        if abs(pari.subst(P, "x", zp)) < 1e-45 and max(abs(pari.polcoef(P, i)) for i in range(d + 1)) < bound:
            return P
    return None


def build_map(f, pl, E, Pi, X_samples=None, attempt=0):
    roots = roots_of([mp.mpf(c) for c in f])
    L = lattice_from_plucker(pl)
    Lp = orth_complement(E, L)
    m1, m2 = unimodular_complement(Lp)
    w = [sum(Pi[r, i] * Lp[0][i] for i in range(4)) for r in range(2)]
    # l(z) = a z1 + b z2 vanishing on C w, normalised so the larger coefficient is 1
    if abs(w[1]) >= abs(w[0]):
        ca, cb = mp.mpf(1), -w[0] / w[1]
    else:
        ca, cb = -w[1] / w[0], mp.mpf(1)
    lam = (ca, cb)
    ell = lambda v: sum((ca * Pi[0, i] + cb * Pi[1, i]) * v[i] for i in range(4))
    for v in Lp:
        assert abs(ell(v)) < mp.mpf(10) ** -40, mp.nstr(abs(ell(v)), 5)
    W1, W2 = ell(m1), ell(m2)
    if mp.im(W1 / W2) < 0:
        W1 = -W1
    lat = pari([pari(f"{mp.nstr(mp.re(W), 68)} + I*({mp.nstr(mp.im(W), 68)})") for W in (W1, W2)])
    g2, g3 = lattice_invariants(W1, W2)
    # base point: a root, preferably real
    r0 = min(roots, key=lambda z: abs(mp.im(z)))
    d = None
    samples = X_samples or [mp.mpc(mp.mpf(a) / 7, mp.mpf(b) / 5) for a, b in
                            [(1, 1), (2, 3), (-3, 2), (5, -1), (-4, -3), (6, 4), (-1, 5), (3, -4), (8, 2), (-6, 1),
                             (2, -6), (-2, -2), (9, -3), (4, 7), (-7, 4), (1, -8), (10, 1), (-9, -5)]]
    scale = max(abs(a - b) for a in roots for b in roots)

    def clearance(xs):
        rest = [z for z in roots if abs(z - r0) > mp.mpf(10) ** -30]
        return min(abs(nearest_point(r0, xs, z) - z) for z in rest) / scale

    if X_samples is None:
        pool = [mp.mpc(mp.mpf(a) / 7, mp.mpf(b) / 5) for a in range(-12, 13) for b in range(-9, 10) if b]
        pool.sort(key=lambda z: -clearance(z))
        samples = pool[18 * attempt:18 * attempt + 18]
    pts = []
    for xs in samples:
        z = aj(f, r0, xs, roots)
        t = ca * z[0] + cb * z[1]
        tp = pari(f"{mp.nstr(mp.re(t), 68)} + I*({mp.nstr(mp.im(t), 68)})")
        Xv = pari.ellwp(lat, tp, precision=640)
        pts.append((xs, mp.mpc(str(pari.real(Xv)).replace(" E", "e"), str(pari.imag(Xv)).replace(" E", "e"))))
    return {"lam": lam, "g2": g2, "g3": g3, "lattice": (W1, W2), "points": pts, "r0": r0, "L": L, "Lperp": Lp}


def interpolate(points, d):
    """Solve N(x) - X D(x) = 0 with deg N, deg D <= d (nullspace by SVD)."""
    rows = []
    for xs, Xv in points:
        rows.append([xs ** i for i in range(d + 1)] + [-Xv * xs ** i for i in range(d + 1)])
    A = mp.matrix(rows)
    U, S, V = mp.svd_c(A)
    v = V.H.column(V.cols - 1) if False else [mp.conj(V[V.rows - 1, j]) for j in range(V.cols)]
    sing = [S[i] for i in range(len(S))]
    N = v[:d + 1]
    D = v[d + 1:]
    k = max(i for i in range(d + 1) if abs(D[i]) > 1e-30)
    s = D[k]
    return [c / s for c in N], [c / s for c in D], sing


if __name__ == "__main__":
    f = json.loads(sys.argv[1])
    X = int(sys.argv[2]) if len(sys.argv) > 2 else 10
    r = complete_spectrum(f, X)
    dmin = min(r["subvarieties"].values())
    pl = next(p for p, dd in r["subvarieties"].items() if dd == dmin)
    m = build_map(f, pl, r["E"], r["Pi"])
    N, D, sing = interpolate(m["points"], dmin)
    print("degree", dmin, "lam", [mp.nstr(t, 30) for t in m["lam"]])
    print("g2", mp.nstr(m["g2"], 30), "g3", mp.nstr(m["g3"], 30))
    print("N", [mp.nstr(c, 25) for c in N])
    print("D", [mp.nstr(c, 25) for c in D])
    print("smallest singular values", [mp.nstr(x, 5) for x in sorted(sing)[:3]])


# ------------------------------------------------------------------ exact certificate
def gaussian_rational(z, bound=10 ** 12, tol=1e-25):
    out = []
    for part in (mp.re(z), mp.im(z)):
        q = pari.bestappr(pari(mp.nstr(part, 60)), bound)
        if abs(mp.mpf(str(pari.numerator(q))) / mp.mpf(str(pari.denominator(q))) - part) > tol:
            return None
        out.append(q)
    return out


def exact_certificate(f, N, D, g2, g3, lam):
    """Recognize coefficients in Q(i) and certify both map coordinates."""
    rec = [gaussian_rational(c) for c in list(N) + list(D) + [g2, g3, lam[0], lam[1]]]
    if any(r is None for r in rec):
        return {"verified": False, "reason": "coefficients not recognised in Q(i)"}

    def el(r):
        return pari(f"({r[0]}) + ({r[1]})*w")

    d = len(N) - 1
    Np = sum((el(rec[i]) * pari("x") ** i for i in range(d + 1)), pari(0))
    Dp = sum((el(rec[d + 1 + i]) * pari("x") ** i for i in range(d + 1)), pari(0))
    G2, G3, A, B = (el(value) for value in rec[-4:])
    record = {"N": str(Np), "D": str(Dp), "g2": str(G2), "g3": str(G3),
              "differential": f"({A}) dx/y + ({B}) x dx/y"}
    try:
        degree, field, _, vertical = certify_record(f, record)
    except ValueError as error:
        return {**record, "verified": False, "reason": str(error)}
    return {**record, "verified": True, "field": field, "Y_over_y": str(vertical),
            "degree_of_R": degree, "map_degree": degree}


def first_quotient(f, X=6):
    """Search numerically for minimal quotients and certify the resulting maps exactly."""
    r = complete_spectrum(f, X)
    if not r["subvarieties"]:
        return {"spectrum": {}, "note": f"no elliptic quotient of degree <= {X}"}
    from collections import Counter
    spec = dict(sorted(Counter(r["subvarieties"].values()).items()))
    dmin = min(spec)
    out = {"spectrum_through_X": {str(k): v for k, v in spec.items()}, "X": X, "minimal_degree": dmin,
           "relation_residual": mp.nstr(r["relation_residual"], 3), "riemann_residual": mp.nstr(r["riemann_residual"], 3),
           "maps": []}
    for pl, dd in r["subvarieties"].items():
        if dd != dmin:
            continue
        for attempt in range(4):
            m = build_map(f, pl, r["E"], r["Pi"], attempt=attempt)
            N, D, sing = interpolate(m["points"], dmin)
            cert = exact_certificate(f, N, D, m["g2"], m["g3"], m["lam"])
            cert["sample_set"] = attempt
            if cert["verified"] and cert["map_degree"] != dd:
                cert.update(verified=False, reason="recognized map has a different degree")
            if cert["verified"]:
                break
        cert["plucker"] = list(pl)
        cert["interpolation_singular_value"] = mp.nstr(sorted(sing)[0], 3)
        out["maps"].append(cert)
    return out
