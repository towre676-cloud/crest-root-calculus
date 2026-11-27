"""Discovery: the two degree-three elliptic covers of the determinant-1440 full-spectrum curves, over Q(j1).

For each Shimura conjugate (b Lambda, H / N b) of the target lattice (b over Pic(O_3), O_3 = Z[3 sqrt -10]):
  1. period matrix, theta constants and Rosenhain branch points (twin-curves/periods.py);
  2. the two degree-3 maps x -> X of the hyperelliptic quotient: f = P/Q with Q the cubic over the triple fibre
     and, over each of the other three Weierstrass points a_i, P - f(a_i) Q = c (x - a_i)(x - b_i)^2;
  3. the map phi1 whose target has j = j(b a1) is normalised canonically: its target is
     E1: Y^2 = X^3 + A X + B, A = 3 j (1728 - j), B = 2 j (1728 - j)^2, and the source coordinate x sends the fourth
     critical point to infinity, the other point of its fibre to 0, and makes Q monic with x^2-coefficient 1, so
     phi1 = v + kappa x / Q(x),  Q = x^3 + x^2 + q1 x + q0;
  4. the complementary map is written in the same coordinate as P2 / R with R the monic cubic over the other three
     Weierstrass points and P2 = x^2 + p1 x + p0.
The six numbers (v, kappa, q1, q0, p1, p0) are recognised in the power basis of theta, a root of the polredabs model of
H_{-360}, by solving the 8 x 8 Vandermonde system over the orbit.  This is discovery only: exact_models.py verifies.

Usage: TWIN_DPS=... python gluing33.py A|B
"""
import json, os, sys, itertools, random
sys.path.insert(0, "../twin-curves"); sys.path.insert(0, "../inverse-spectrum")
import orbit as O
import mpmath as mp
from fractions import Fraction as Fr
from gluing import K
from polarised import ideals

pari = O.pari
DPS = O.DPS
pari('default(format, "g.%d")' % (DPS - 20))                         # full-precision str() of PARI reals
D = -40


def cplx(Kf, e):
    w = (mp.mpf(D) + mp.sqrt(mp.mpf(D))) / 2
    return mp.mpf(e[0].numerator) / e[0].denominator + mp.mpf(e[1].numerator) / e[1].denominator * w


def jinv(Kf, basis):
    t = cplx(Kf, basis[1]) / cplx(Kf, basis[0])
    if t.imag < 0:
        t = -t
    return pari.ellj(pari("%s + I*(%s)" % (mp.nstr(t.real, DPS - 10), mp.nstr(t.imag, DPS - 10))))


def members(rec, f=3):
    Kf = K(D)
    f1, a1, f2, a2 = rec["product"]
    a1 = [tuple(Fr(x) for x in e) for e in a1]; a2 = [tuple(Fr(x) for x in e) for e in a2]
    H = tuple(Fr(x) for x in rec["H"])
    for b in ideals(D, f):
        b = [tuple(e) for e in b]
        Nb = O.covol(Kf, b) / O.covol(Kf, [(Fr(1), Fr(0)), (Fr(0), Fr(f))])
        A1, A2 = O.lattice_mult(Kf, b, a1), O.lattice_mult(Kf, b, a2)
        Bs = [(A1[0][0], A1[0][1], Fr(0), Fr(0)), (A1[1][0], A1[1][1], Fr(0), Fr(0)),
              (Fr(0), Fr(0), A2[0][0], A2[0][1]), (Fr(0), Fr(0), A2[1][0], A2[1][1])]
        tau = O.siegel_reduce(O.period_matrix(D, Bs, tuple(x / Nb for x in H)))
        lam, mu, nu = O.rosenhain(O.thetas(tau))
        yield {"branch": [mp.mpc(0), mp.mpc(1), lam, mu, nu], "j1": jinv(Kf, A1), "j2": jinv(Kf, A2),
               "shifted_j1": [jinv(Kf, O.lattice_mult(Kf, [tuple(e) for e in c], A1)) for c in ideals(D, f)]}


def tompc(g):
    return mp.mpc(mp.mpf(str(pari.real(g)).replace(" ", "")), mp.mpf(str(pari.imag(g)).replace(" ", "")))


def ev(p, x):
    s = mp.mpc(0)
    for c in reversed(p):
        s = s * x + c
    return s


def fromroots(rs):
    p = [mp.mpc(1)]
    for r in rs:
        q = [mp.mpc(0)] * (len(p) + 1)
        for i, c in enumerate(p):
            q[i] -= r * c; q[i + 1] += c
        p = q
    return p


def divlin(p, a):
    n = len(p) - 1; q = [mp.mpc(0)] * n; r = p[n]
    for k in range(n - 1, -1, -1):
        q[k] = r; r = p[k] + a * r
    return q


def mkP(a, b, mode):
    return [b, a, mp.mpc(1), mp.mpc(0)] if mode == 0 else [b, mp.mpc(1), a, mp.mpc(0)]


def residuals(Q, T, p1, p0, mode=0):
    P = mkP(p1, p0, mode)
    out = []
    for a in T:
        e = ev(P, a) / ev(Q, a)
        h = divlin([P[k] - e * Q[k] for k in range(4)], a)
        out.append(h[1] ** 2 - 4 * h[0] * h[2])
    return out


def find_maps(W):
    """All (triple, singles, P) at 50 digits, then refined at full precision; P normalised with x^3-coefficient 0
    and x^2-coefficient 1 (mode 0) or x-coefficient 1 (mode 1)."""
    found = {}
    with mp.workdps(50):
        Wl = [+w for w in W]
        for S in itertools.combinations(range(6), 3):
            T = [i for i in range(6) if i not in S]
            Q = fromroots([Wl[i] for i in S]); Ts = [Wl[i] for i in T]
            for mode in (0, 1):
                if S in found:
                    break
                for _ in range(60):
                    st = (mp.mpc(random.gauss(0, 3), random.gauss(0, 3)), mp.mpc(random.gauss(0, 3), random.gauss(0, 3)))
                    try:
                        r = mp.findroot(lambda a, b: residuals(Q, Ts, a, b, mode)[:2], st, tol=mp.mpf(10) ** -35, maxsteps=60)
                    except Exception:
                        continue
                    if abs(residuals(Q, Ts, r[0], r[1], mode)[2]) < mp.mpf(10) ** -25:
                        found[S] = (T, mode, r[0], r[1])
                        break
    out = []
    for S, (T, mode, a, b) in found.items():
        Q = fromroots([W[i] for i in S]); Ts = [W[i] for i in T]
        r = mp.findroot(lambda x, y: residuals(Q, Ts, x, y, mode)[:2], (a, b), tol=mp.mpf(10) ** -(DPS - 30), maxsteps=200)
        assert abs(residuals(Q, Ts, r[0], r[1], mode)[2]) < mp.mpf(10) ** -(DPS - 60)
        out.append((list(S), T, mkP(r[0], r[1], mode)))
    return out


def jcubic(es):
    s1 = sum(es); s2 = es[0] * es[1] + es[0] * es[2] + es[1] * es[2]; s3 = es[0] * es[1] * es[2]
    # depress X^3 - s1 X^2 + s2 X - s3
    a = s2 - s1 ** 2 / 3; b = -s3 + s1 * s2 / 3 - 2 * s1 ** 3 / 27
    return 1728 * 4 * a ** 3 / (4 * a ** 3 + 27 * b ** 2)


def polyroots(p):
    return mp.polyroots(list(reversed(p)), maxsteps=400, extraprec=4 * DPS)


def normalise(m):
    W = m["branch"]
    c0 = mp.mpc("0.31", "0.77")
    W = [1 / (w - c0) for w in W] + [mp.mpc(0)]                       # x = inf -> 0, all six finite
    maps = find_maps(W)
    tries = 0
    while len(maps) < 2 and tries < 5:
        tries += 1
        maps = find_maps(W)
    assert len(maps) == 2, f"expected exactly two degree-3 maps, found {len(maps)}"
    j1 = tompc(m["j1"])
    data = []
    for S, T, P in maps:
        Q = fromroots([W[i] for i in S])
        es = [ev(P, W[i]) / ev(Q, W[i]) for i in T]
        data.append((S, T, P, Q, es, jcubic(es)))
    data.sort(key=lambda d: abs(d[5] - j1))
    assert abs(data[0][5] - j1) < mp.mpf(10) ** -(DPS // 2) * abs(j1), "phi1 target j mismatch"
    S, T, P, Q, es, _ = data[0]
    S2, T2, P2, Q2, es2, j2c = data[1]
    assert sorted(S2) == sorted(T)
    # critical points of f = P/Q: roots of P'Q - P Q'
    dP = [P[1], 2 * P[2], 3 * P[3]]; dQ = [Q[1], 2 * Q[2], 3 * Q[3]]
    num = [mp.mpc(0)] * 6
    for i, x in enumerate(dP):
        for k, y in enumerate(Q):
            num[i + k] += x * y
    for i, x in enumerate(P):
        for k, y in enumerate(dQ):
            num[i + k] -= x * y
    while abs(num[-1]) < mp.mpf(10) ** -(DPS - 40):
        num.pop()
    crit = polyroots(num)
    f = lambda x: ev(P, x) / ev(Q, x)
    c = max(crit, key=lambda x: min(abs(f(x) - e) for e in es))
    vprime = f(c)
    # other point of the fibre over vprime: roots of P - v' Q (double root at c)
    cub = [P[k] - vprime * Q[k] for k in range(4)]
    cp = -cub[2] / cub[3] - 2 * c                                       # P - v'Q = lead (x - c)^2 (x - c')
    # affine identification of the target X-line with E1: X^3 + A X + B
    jj = j1
    A = 3 * jj * (1728 - jj); B = 2 * jj * (1728 - jj) ** 2
    ge = polyroots([B, A, mp.mpc(0), mp.mpc(1)])
    best = None
    for perm in itertools.permutations(range(3)):
        t = [ge[k] for k in perm]
        a = (t[1] - t[0]) / (es[1] - es[0]); d = t[0] - a * es[0]
        err = abs(a * es[2] + d - t[2])
        if best is None or err < best[0]:
            best = (err, a, d)
    assert best[0] < mp.mpf(10) ** -(DPS // 2) * (1 + abs(A))
    _, a, d = best
    v = a * vprime + d
    # source coordinate: x_new = lam (x - cp)/(x - c)
    mob = lambda x: (x - cp) / (x - c)
    qr = [mob(W[i]) for i in S]
    lam = -1 / sum(qr)
    qr = [lam * r for r in qr]
    Qn = fromroots(qr)                                                 # monic, x^2 coefficient 1
    pt = mp.mpc("0.123", "0.456")                                      # a test point in the new coordinate
    xold = (c * pt / lam - cp) / (pt / lam - 1)
    kappa = (a * f(xold) + d - v) * ev(Qn, pt) / pt
    # complementary map in the new coordinate
    rr = [lam * mob(W[i]) for i in T]
    Rn = fromroots(rr)
    f2 = lambda x: ev(P2, x) / ev(Q2, x)
    pts = [mp.mpc(k, 1) / 7 for k in range(1, 5)]
    vals = []
    for p in pts:
        xo = (c * p / lam - cp) / (p / lam - 1)
        vals.append(f2(xo) * ev(Rn, p))
    M = mp.matrix([[p ** k for k in range(4)] for p in pts])
    P2n = list(mp.lu_solve(M, mp.matrix(vals)))
    P2n = [P2n[k] - P2n[3] * Rn[k] for k in range(4)]
    P2n = [x / P2n[2] for x in P2n]
    return {"v": v, "kappa": kappa, "q1": Qn[1], "q0": Qn[0], "p1": P2n[1], "p0": P2n[0],
            "j1": m["j1"], "j2_from_map": j2c, "j2_lattice": m["j2"], "shifted_j1": m["shifted_j1"]}


if __name__ == "__main__":
    random.seed(7)
    name = sys.argv[1]
    rec = json.load(open("lattices_1440.json"))[name]
    import pickle
    cache = f"members_{name}_{DPS}.pkl"
    if os.path.exists(cache):
        ms = pickle.load(open(cache, "rb"))
    else:
        ms = list(members(rec))
        pickle.dump(ms, open(cache, "wb"))
    rows = []
    for m in ms:
        r = normalise(m)
        rows.append(r)
        print(name, "member", len(rows), "| j2 match %.1e" % float(abs(r["j2_from_map"] - tompc(r["j2_lattice"])) / abs(r["j2_from_map"])),
              "| |v| %.2e |kappa| %.2e" % (float(abs(r["v"])), float(abs(r["kappa"]))), flush=True)
    keys = ("v", "kappa", "q1", "q0", "p1", "p0")
    Hp = pari.polclass(-360)
    F, a2 = pari.polredabs(Hp, 1)
    phi = pari.lift(a2)                                                # j1 = phi(theta)
    rootsF = pari.polroots(F, precision=4 * DPS)
    P = lambda z: pari("%s + I*(%s)" % (mp.nstr(z.real, DPS - 30), mp.nstr(z.imag, DPS - 30)))
    sq = pari("sqrt(10)")
    def system(shift):
        th = []
        for r in rows:
            jt = r["shifted_j1"][shift]
            cand = sorted(rootsF, key=lambda t: abs(pari.subst(phi, "x", t) - jt) / abs(jt))
            th.append(cand[0])
        if len({str(pari.precision(t, 20)) for t in th}) != 8:
            return None
        ent = []
        for part in ("re", "im"):
            for i in range(8):
                row = []
                for k in range(8):
                    z = th[i] ** k
                    row.append(pari.real(z) if part == "re" else pari.imag(z))
                for k in range(8):
                    z = th[i] ** k
                    row.append(-sq * pari.imag(z) if part == "re" else sq * pari.real(z))
                ent += row
        return pari.matrix(16, 16, ent)
    def solve(Vm, key):
        rhs = pari.matrix(16, 1, [pari.real(P(r[key])) for r in rows] + [pari.imag(P(r[key])) for r in rows])
        c = pari.matsolve(Vm, rhs)
        coeffs, stable = [], 0
        for i in range(16):
            x = c[i, 0]
            q1 = pari.bestappr(x, pari("10^%d" % (DPS // 3)))
            q2 = pari.bestappr(x, pari("10^%d" % (DPS // 4)))
            stable += int(q1 == q2 and abs(x - q1) < pari("10^-%d" % int(0.55 * DPS)) * (1 + abs(x)))
            coeffs.append(q1)
        return coeffs, stable
    best = None
    for shift in range(8):
        Vm = system(shift)
        if Vm is None:
            continue
        _, st = solve(Vm, "p1")
        st2 = solve(Vm, "v")[1]
        print("shift", shift, "stable coefficients: p1", st, "v", st2, flush=True)
        if best is None or st + st2 > best[0]:
            best = (st + st2, shift, Vm)
    _, shift, Vm = best
    out = {"F": str(F), "j1_in_theta": str(phi), "shift_index": shift, "extension": "sqrt(-10)",
           "form": "value = alpha(theta) + beta(theta) * sqrt(-10)", "values": {}}
    for k in keys:
        coeffs, worst = solve(Vm, k)
        worst = f"{worst}/16 stable"
        X = pari("x")
        alpha = sum((coeffs[i] * X ** i for i in range(8)), pari(0))
        beta = sum((coeffs[8 + i] * X ** i for i in range(8)), pari(0))
        digits = max(len(str(q)) for q in coeffs)
        out["values"][k] = {"alpha": str(alpha), "beta": str(beta), "stability": worst, "max_coefficient_chars": digits}
        print(name, k, worst, "max coefficient chars", digits, flush=True)
    json.dump(out, open(f"recognised_{name}.json", "w"), indent=1)
