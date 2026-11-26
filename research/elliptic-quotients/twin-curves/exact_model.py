"""Canonical bielliptic model y^2 = x^6 + alpha x^4 + alpha x^2 + beta, oriented so that Y^2 = U^3 + A U^2 + B U + C
is the quotient with real j; alpha = A^2/B, beta = C A^3/B^3; exact recognition."""
import json, sys
import mpmath as mp
import cypari2
from bielliptic_model import models

pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9); pari.set_real_precision(250)
mp.mp.dps = 270
P = lambda z: pari("%s + I*(%s)" % (mp.nstr(z.real, 245), mp.nstr(z.imag, 245)))


def recog(zp, maxdeg=16, tol="10^-215"):
    if abs(pari.imag(zp)) < pari("10^-200"):
        r = pari.bestappr(pari.real(zp), pari("10^110"))
        if abs(pari.real(zp) - r) < pari("10^-200"):
            return {"degree": 1, "value": str(r)}
    for deg in range(2, maxdeg + 1):
        Q = pari.algdep(zp, deg)
        if abs(pari.subst(Q, "x", zp)) / (abs(zp) ** deg + 1) < pari(tol):
            F = pari.factor(Q)
            for k in range(int(pari.matsize(F)[0])):
                f = F[k, 0]
                if abs(pari.subst(f, "x", zp)) / (abs(zp) ** int(pari.poldegree(f)) + 1) < pari(tol):
                    return {"degree": int(pari.poldegree(f)), "minpoly": str(f)}
    return None


d = int(sys.argv[1]); fn = sys.argv[2]
C = json.load(open(fn))
out = {}
for name in ("q1", "q2"):
    ros = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in C[name]["rosenhain"]]
    (alpha, beta, (A, B, Cc)), = models([mp.mpc(0), mp.mpc(1)] + ros)
    E1 = pari.ellinit([0, P(A), 0, P(B), P(Cc)])
    if abs(pari.imag(E1[12])) > pari("10^-100"):            # reorient: reverse the cubic
        A, B, Cc = B / Cc, A / Cc, 1 / Cc
    alpha, beta = A * A / B, Cc * A ** 3 / B ** 3
    ra, rb = recog(P(alpha)), recog(P(beta))
    out[name] = {"alpha": ra, "beta": rb, "alpha_num": [mp.nstr(alpha.real, 60), mp.nstr(alpha.imag, 60)],
                 "beta_num": [mp.nstr(beta.real, 60), mp.nstr(beta.imag, 60)]}
    print(name, "alpha:", ra, "\n    beta:", rb, flush=True)
json.dump(out, open(fn.replace("curves", "exact_model"), "w"), indent=1)
