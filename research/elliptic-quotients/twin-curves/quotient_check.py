"""j-invariants of the two degree-2 elliptic quotients of each twin, tested against Hilbert class polynomials
of orders of Q(sqrt d) (polclass), and the exact j-pair data."""
import json, sys
import mpmath as mp
import cypari2
from bielliptic_model import models

pari = cypari2.Pari()
pari.allocatemem(4 * 10 ** 9)
pari.set_real_precision(250)
mp.mp.dps = 270
d = int(sys.argv[1]); fn = sys.argv[2]
C = json.load(open(fn))
P = lambda z: pari("%s + I*(%s)" % (mp.nstr(z.real, 245), mp.nstr(z.imag, 245)))
discs = [d * f * f for f in (1, 2, 3, 4, 6, 9)]
H = {D: pari.polclass(D) for D in discs}
out = {}
for name in ("q1", "q2"):
    ros = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in C[name]["rosenhain"]]
    for alpha, beta, (A, B, Cc) in models([mp.mpc(0), mp.mpc(1)] + ros):
        E1 = pari.ellinit([0, P(A), 0, P(B), P(Cc)])
        E2 = pari.ellinit([0, P(B / Cc), 0, P(A / Cc), P(1 / Cc)])   # reversed cubic, made monic
        js = [pari.member_j(E1) if hasattr(pari, "member_j") else E1[12], E2[12]]
        hits = []
        for j in (E1[12], E2[12]):
            hit = None
            for D, h in H.items():
                if abs(pari.subst(h, "x", j)) < pari("10^-100") * (1 + abs(j)) ** int(pari.poldegree(h)):
                    hit = D; break
            hits.append(hit)
        out[name] = {"quotient_cm_discriminants": hits, "j": [str(pari.precision(E1[12], 20)), str(pari.precision(E2[12], 20))]}
        print(name, out[name], flush=True)
json.dump(out, open(fn.replace("curves", "quotients"), "w"), indent=1)
