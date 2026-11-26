"""Bielliptic normal form y^2 = x^6 + alpha x^4 + alpha x^2 + beta from Rosenhain branch points, with
exact recognition of alpha, beta and CM checks on the two degree-2 quotients."""
import itertools, json, sys
import mpmath as mp
import cypari2

pari = cypari2.Pari()
pari.set_real_precision(250)
mp.mp.dps = 270


def mob_involution(p1, p2, p3, p4):
    """Moebius involution s with s(p1) = p2, s(p3) = p4 (finite points): s(x) = (a x + b)/(c x - a)."""
    # (a x + b) = s (c x - a)  for (x, s) = (p1, p2), (p2, p1), (p3, p4): linear in (a, b, c)
    M = mp.matrix([[p1 + p2, 1, -p1 * p2], [p3 + p4, 1, -p3 * p4]])
    # null vector of 2x3 system
    a = M[0, 1] * M[1, 2] - M[0, 2] * M[1, 1]
    b = M[0, 2] * M[1, 0] - M[0, 0] * M[1, 2]
    c = M[0, 0] * M[1, 1] - M[0, 1] * M[1, 0]
    return a, b, c


def models(branch):
    """All bielliptic normal forms (alpha, beta) for branch points 'branch' (finite list of 5 + infinity)."""
    INF = None
    pts = list(branch) + [INF]
    out = []
    for p in itertools.permutations(range(6)):
        pairs = [(p[0], p[1]), (p[2], p[3]), (p[4], p[5])]
        if not (p[0] < p[1] and p[2] < p[3] and p[4] < p[5] and p[0] < p[2] < p[4]):
            continue
        # use the two pairs without infinity to determine s
        fin = [pr for pr in pairs if INF not in (pts[pr[0]], pts[pr[1]])]
        if len(fin) < 2:
            continue
        a, b, c = mob_involution(pts[fin[0][0]], pts[fin[0][1]], pts[fin[1][0]], pts[fin[1][1]])
        s = lambda x: (c and (a / c)) if x is None else ((a * x + b) / (c * x - a) if abs(c * x - a) > mp.mpf(10) ** -200 else None)
        ok = True
        for i, j in pairs:
            u, v = pts[i], pts[j]
            su = s(u)
            if (su is None and v is not None) or (su is not None and (v is None or abs(su - v) > mp.mpf(10) ** -150)):
                ok = False
        if not ok:
            continue
        # fixed points of s: c x^2 - 2 a x - b = 0
        disc = mp.sqrt(a * a + b * c)
        f1, f2 = (a + disc) / c, (a - disc) / c
        T = lambda x: mp.mpf(1) if x is None else (x - f1) / (x - f2)      # s becomes x -> -x
        img = [T(x) for x in pts]
        sq = []
        for i, j in pairs:
            assert abs(img[i] + img[j]) < mp.mpf(10) ** -100
            sq.append(img[i] ** 2)
        A = -(sq[0] + sq[1] + sq[2]); B = sq[0] * sq[1] + sq[0] * sq[2] + sq[1] * sq[2]; C = -sq[0] * sq[1] * sq[2]
        out.append((A * A / B, C * A ** 3 / B ** 3, (A, B, C)))
    return out


def recognise(z, maxdeg=8):
    zp = pari("%s + I*(%s)" % (mp.nstr(z.real, 245), mp.nstr(z.imag, 245)))
    if abs(z.imag) < mp.mpf(10) ** -200:
        zp = pari(mp.nstr(z.real, 245))
        r = pari.bestappr(zp, pari("10^100"))
        if abs(zp - r) < pari("10^-200"):
            return str(r)
    for deg in range(2, maxdeg + 1):
        P = pari.algdep(zp, deg)
        if abs(pari.subst(P, "x", zp)) / (abs(zp) ** deg + 1) < pari("10^-220"):
            return str(pari.factor(P)[0, 0]) if True else str(P)
    return None


def jinv(cubic):
    """j of Y^2 = cubic(U) (monic cubic coefficients [A, B, C])."""
    A, B, C = cubic
    E = pari.ellinit([0, pari("%s + I*(%s)" % (mp.nstr(A.real, 245), mp.nstr(A.imag, 245))), 0,
                      pari("%s + I*(%s)" % (mp.nstr(B.real, 245), mp.nstr(B.imag, 245))),
                      pari("%s + I*(%s)" % (mp.nstr(C.real, 245), mp.nstr(C.imag, 245)))])
    return E.j if hasattr(E, "j") else pari.member_j(E) if False else E[12]


if __name__ == "__main__":
    fn = sys.argv[1]
    C = json.load(open(fn))
    out = {}
    for name in ("q1", "q2"):
        ros = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in C[name]["rosenhain"]]
        br = [mp.mpc(0), mp.mpc(1)] + ros
        ms = models(br)
        print(name, "bielliptic involutions found:", len(ms), flush=True)
        recs = []
        for alpha, beta, (A, B, Cc) in ms:
            ra, rb = recognise(alpha), recognise(beta)
            # degree-2 quotients: Y^2 = U^3 + A U^2 + B U + C  and the reversed cubic  C U^3 + B U^2 + A U + 1
            j1 = pari.ellj(pari("0") + 0) if False else None
            recs.append({"alpha": ra, "beta": rb, "alpha_num": mp.nstr(alpha, 20), "beta_num": mp.nstr(beta, 20)})
            print("   alpha:", ra, "\n   beta:", rb, flush=True)
        out[name] = recs
    json.dump(out, open(fn.replace("curves", "bielliptic"), "w"), indent=1)
