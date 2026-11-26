"""Lattice + polarisation -> reduced tau -> theta constants -> Rosenhain model -> Igusa invariants."""
import json, sys
sys.path.insert(0, "."); sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
import mpmath as mp
import cypari2
from periods import period_matrix, thetas, rosenhain, siegel_reduce

pari = cypari2.Pari()
pari.set_real_precision(260)
mp.mp.dps = 280


def load(name, fn):
    t = json.load(open(fn))[name]
    f1, a1, f2, a2 = t["product"]
    a1 = [[Fr(x) for x in e] for e in a1]; a2 = [[Fr(x) for x in e] for e in a2]
    B = [(a1[0][0], a1[0][1], Fr(0), Fr(0)), (a1[1][0], a1[1][1], Fr(0), Fr(0)),
         (Fr(0), Fr(0), a2[0][0], a2[0][1]), (Fr(0), Fr(0), a2[1][0], a2[1][1])]
    return B, tuple(Fr(x) for x in t["H"])


def curve(d, B, H):
    tau = siegel_reduce(period_matrix(d, B, H))
    th = thetas(tau)
    lam, mu, nu = rosenhain(th)
    s = lambda z: "(%s + I*(%s))" % (mp.nstr(z.real, 265), mp.nstr(z.imag, 265))
    P = pari("x*(x-1)*(x-%s)*(x-%s)*(x-%s)" % (s(lam), s(mu), s(nu)))
    J = pari.genus2igusa(P)
    return tau, th, (lam, mu, nu), P, J


if __name__ == "__main__":
    d, fn = int(sys.argv[1]), sys.argv[2]
    out = {}
    for name in ("q1", "q2"):
        B, H = load(name, fn)
        tau, th, ros, P, J = curve(d, B, H)
        # absolute invariants (weights: J2:2, J4:4, J6:6, J8:8, J10:10)
        J2, J4, J6, J8, J10 = J
        abs_inv = [J2 ** 5 / J10, J2 ** 3 * J4 / J10, J2 ** 2 * J6 / J10]
        print(name, "reduced tau:", mp.nstr(tau, 8))
        print("   rosenhain:", [mp.nstr(z, 12) for z in ros])
        print("   absolute Igusa (J2^5/J10, J2^3 J4/J10, J2^2 J6/J10):", [str(pari.precision(x, 20)) for x in abs_inv])
        out[name] = {"tau": [[[str(tau[i, j].real), str(tau[i, j].imag)] for j in range(2)] for i in range(2)],
                     "rosenhain": [[str(z.real), str(z.imag)] for z in ros],
                     "abs_igusa": [str(x) for x in abs_inv]}
    json.dump(out, open(fn.replace("twin_lattices", "curves"), "w"), indent=1)
