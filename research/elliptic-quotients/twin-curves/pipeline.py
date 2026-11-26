"""Full lattice -> curve pipeline for every located twin form of one discriminant.
Steps (each exact or with explicit acceptance thresholds):
  1. tau (Siegel-reduced) from (Lambda, H); theta constants; Rosenhain model; absolute Igusa invariants (260 digits).
  2. field of moduli F: minimal polynomial of z = J2^5/J10 by algdep (degree <= 16), accepted only if the
     residual < 10^-200 AND the coefficient height is far below the random-relation height; F = polredbest.
  3. the three absolute invariants as exact elements of F (lindep over a small generator, height check).
  4. if q represents 4 (bielliptic): dihedral invariants u, v in F, model y^2 = x^6 - s1 x^4 + s2 x^2 - 1 over
     F(s1), s1^6 - u s1^3 + v^3 = 0, and an EXACT check that its Igusa invariants equal step 3."""
import json, sys
sys.path.insert(0, "."); sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
import mpmath as mp
import cypari2
from periods import period_matrix, thetas, rosenhain, siegel_reduce
from bielliptic_model import models

import os
DPS = int(os.environ.get('TWIN_DPS', '280'))
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9); pari.set_real_precision(DPS - 20)
mp.mp.dps = DPS
RECOGNITION_TOLERANCE = pari(f"10^-{DPS - 80}")
P = lambda z: pari("%s + I*(%s)" % (mp.nstr(z.real, DPS - 30), mp.nstr(z.imag, DPS - 30)))
GEN = pari("genus2igusa(x^6 - a*x^4 + b*x^2 - 1)")


def lindep_F(t, theta, n, maxh=None):
    basis = [t] + [theta ** k for k in range(n)]
    v = pari.lindep(basis)
    res = abs(sum((v[i] * basis[i] for i in range(n + 1)), pari(0)))
    h = max(abs(int(x)) for x in v)
    rand = 2 * (DPS - 30) / (n + 1)                  # random-relation height (decimal digits)
    if v[0] != 0 and res < RECOGNITION_TOLERANCE and len(str(h)) < 0.6 * rand:
        return str(pari("-(%s)/(%s)" % (" + ".join(f"({int(v[i])})*x^{i-1}" for i in range(1, n + 1)), int(v[0]))))
    return None


def run(d, name, rec):
    f1, a1, f2, a2 = rec["product"]
    a1 = [[Fr(x) for x in e] for e in a1]; a2 = [[Fr(x) for x in e] for e in a2]
    B = [(a1[0][0], a1[0][1], Fr(0), Fr(0)), (a1[1][0], a1[1][1], Fr(0), Fr(0)),
         (Fr(0), Fr(0), a2[0][0], a2[0][1]), (Fr(0), Fr(0), a2[1][0], a2[1][1])]
    H = tuple(Fr(x) for x in rec["H"])
    tau = siegel_reduce(period_matrix(d, B, H))
    lam, mu, nu = rosenhain(thetas(tau))
    s = lambda z: "(%s + I*(%s))" % (mp.nstr(z.real, DPS - 15), mp.nstr(z.imag, DPS - 15))
    J2, J4, J6, J8, J10 = pari.genus2igusa(pari("x*(x-1)*(x-%s)*(x-%s)*(x-%s)" % (s(lam), s(mu), s(nu))))
    absinv = [J2 ** 5 / J10, J2 ** 3 * J4 / J10, J2 ** 2 * J6 / J10]
    z = absinv[0]
    out = {"product": rec["product"], "H": rec["H"], "q": rec["q"]}
    F = None
    for deg in range(2, 17):
        Q = pari.algdep(z, deg)
        res = abs(pari.subst(Q, "x", z)) / (abs(z) ** deg + 1)
        if res < RECOGNITION_TOLERANCE:
            Fq = pari.factor(Q)
            for k in range(int(pari.matsize(Fq)[0])):
                f = Fq[k, 0]
                if abs(pari.subst(f, "x", z)) / (abs(z) ** int(pari.poldegree(f)) + 1) < RECOGNITION_TOLERANCE:
                    # genuine only if the field is small: reduced polynomial with small coefficients
                    if int(pari.poldegree(f)) <= 16:
                        red = pari.polredabs(f) if int(pari.poldegree(f)) <= 8 else pari.polredbest(f)
                        if max(abs(int(pari.polcoef(red, j))) for j in range(int(pari.poldegree(red)) + 1)) < 10 ** 8:
                            F = f
                    break
            if F is not None:
                break
    if F is None:
        out["status"] = "no accepted field candidate through degree 16"; return out
    Qr, amod = pari.polredbest(F, 1)
    a = pari.lift(amod)
    theta = min(pari.polroots(Qr, precision=4 * DPS), key=lambda r: abs(pari.subst(a, "x", r) - z))
    out["F"] = str(Qr); out["F_disc"] = str(pari.factor(pari.nfdisc(Qr)))
    names = ["J2^5/J10", "J2^3J4/J10", "J2^2J6/J10"]
    nF = int(pari.poldegree(Qr))
    exact = [lindep_F(t, theta, nF) for t in absinv]
    out["abs_igusa_in_F"] = dict(zip(names, exact))
    biell = int(pari.qfminim(pari.matrix(3, 3, [x for r in rec["q"] for x in r]), 4, 0, 2)[0]) > 0
    out["bielliptic"] = biell
    if not (biell and all(exact)):
        out["status"] = "exact invariants only" if all(exact) else "recognition incomplete"; return out
    ms = models([mp.mpc(0), mp.mpc(1), lam, mu, nu])
    (al, be, (A, Bc, C)) = ms[0]
    u, v = lindep_F(P(A ** 3 / C + Bc ** 3 / C ** 2), theta, nF), lindep_F(P(A * Bc / C), theta, nF)
    out["u"], out["v"] = u, v
    if not (u and v):
        out["status"] = "dihedral invariants not recognised"; return out
    Ft = Qr.subst("x", "t"); nf = pari.nfinit(Ft)
    el = lambda e: pari.Mod(pari(e).subst("x", "t"), Ft)
    uu, vv = el(u), el(v)
    h = pari.nffactor(nf, pari.lift(pari("x^6") - uu * pari("x^3") + vv ** 3))[0, 0]
    zpol, aa, k = pari.rnfequation(nf, h, 1)
    aa = pari.Mod(pari.lift(aa), zpol)
    s1 = pari.Mod(pari("x"), zpol) - int(k) * aa
    toz = lambda X: pari.Mod(pari.subst(pari.lift(X), "t", pari.lift(aa)), zpol)
    s2 = toz(vv) / s1
    Jm = [pari.substvec(g, ["a", "b"], [s1, s2]) for g in GEN]
    model_abs = [Jm[0] ** 5 / Jm[4], Jm[0] ** 3 * Jm[1] / Jm[4], Jm[0] ** 2 * Jm[2] / Jm[4]]
    target = [pari.Mod(pari.subst(pari(e), "x", pari.lift(aa)), zpol) for e in exact]
    out["model_field"] = str(pari.polredbest(zpol)); out["model_field_degree"] = int(pari.poldegree(zpol))
    out["model"] = "y^2 = x^6 - s1 x^4 + s2 x^2 - 1"
    out["s1"] = str(pari.lift(s1)); out["s2"] = str(pari.lift(s2)); out["model_field_raw"] = str(zpol)
    out["igusa_exact_match"] = bool(all(pari.lift(m - t) == 0 for m, t in zip(model_abs, target)))
    out["status"] = "explicit model verified exactly" if out["igusa_exact_match"] else "model mismatch"
    return out


if __name__ == "__main__":
    d = int(sys.argv[1])
    L = json.load(open(f"twin_lattices_{d}.json"))
    res = {}
    for name in sorted(L):
        r = run(d, name, L[name])
        res[name] = r
        print(name, r.get("status"), "| F:", r.get("F"), "| disc", r.get("F_disc"), "| bielliptic:", r.get("bielliptic"),
              "| model field:", r.get("model_field"), flush=True)
        json.dump(res, open(f"twin_curves_{d}.json", "w"), indent=1)
