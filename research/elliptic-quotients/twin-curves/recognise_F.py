"""Numerical discovery of invariants in a quartic power basis.

For each form, the orbit gives (z_i, w_i), i = 1..4.  theta_i is the conjugate of theta with phi(theta_i) = z_i, where
phi expresses a root of the minimal polynomial of z in theta.  Then w = sum c_k theta^k with c solving the 4 x 4
Vandermonde system in the theta_i (well conditioned: theta has small height).  Each c_k is rounded by bestappr and
accepted as a candidate only if it meets the working error tolerance and
its denominator is smooth (all prime factors below 10^4). The candidate
filter does not prove equality with the analytic orbit values."""
import json, os, sys
sys.path.insert(0, ".")
import mpmath as mp
import cypari2
import orbit as O

pari = O.pari
DPS = O.DPS
mp.mp.dps = DPS
ACC = int(0.6 * DPS)                # accepted accuracy; measured delivery is ~0.9 DPS once precisions are set
from bielliptic_model import models
from fractions import Fraction as Fr
from gluing import K
from polarised import ideals
from periods import period_matrix, thetas, rosenhain, siegel_reduce
# imported modules set their own precisions (bielliptic_model: 270 / 250); restore the requested one
mp.mp.dps = DPS
pari.set_real_precision(DPS - 20)


def member_values(d, rec, f=3):
    """[z, J2^3J4/J10, J2^2J6/J10, u, v] for each orbit member (u, v only when bielliptic)."""
    Kf = K(d)
    f1, a1, f2, a2 = rec["product"]
    a1 = [tuple(Fr(x) for x in e) for e in a1]; a2 = [tuple(Fr(x) for x in e) for e in a2]
    H = tuple(Fr(x) for x in rec["H"])
    P = lambda z: pari("%s + I*(%s)" % (mp.nstr(z.real, DPS - 30), mp.nstr(z.imag, DPS - 30)))
    s = lambda z: "(%s + I*(%s))" % (mp.nstr(z.real, DPS - 15), mp.nstr(z.imag, DPS - 15))
    out = []
    for b in ideals(d, f):
        b = [tuple(e) for e in b]
        Nb = O.covol(Kf, b) / O.covol(Kf, [(Fr(1), Fr(0)), (Fr(0), Fr(f))])
        A1, A2 = O.lattice_mult(Kf, b, a1), O.lattice_mult(Kf, b, a2)
        B = [(A1[0][0], A1[0][1], Fr(0), Fr(0)), (A1[1][0], A1[1][1], Fr(0), Fr(0)),
             (Fr(0), Fr(0), A2[0][0], A2[0][1]), (Fr(0), Fr(0), A2[1][0], A2[1][1])]
        tau = siegel_reduce(period_matrix(d, B, tuple(x / Nb for x in H)))
        lam, mu, nu = rosenhain(thetas(tau))
        J2, J4, J6, J8, J10 = pari.genus2igusa(pari("x*(x-1)*(x-%s)*(x-%s)*(x-%s)" % (s(lam), s(mu), s(nu))))
        row = {"z": J2 ** 5 / J10, "J2^3J4/J10": J2 ** 3 * J4 / J10, "J2^2J6/J10": J2 ** 2 * J6 / J10}
        ms = models([mp.mpc(0), mp.mpc(1), lam, mu, nu])
        row["n_involutions"] = len(ms)
        if len(ms) == 1:
            (_, _, (A, Bc, C)) = ms[0]
            row["u"] = P(A ** 3 / C + Bc ** 3 / C ** 2); row["v"] = P(A * Bc / C)
        out.append(row)
    return out


def smooth(n, bound=10 ** 4):
    n = abs(int(n))
    for p in [int(p) for p in pari.primes(1229)]:
        while n % p == 0:
            n //= p
    return n == 1


def recognise(F, thetas_, ws):
    V = pari.matrix(4, 4, [thetas_[i] ** k for i in range(4) for k in range(4)])
    c = pari.matsolve(V, pari.matrix(4, 1, ws))
    coeffs, ok = [], True
    for k in range(4):
        x = c[k, 0]
        if abs(pari.imag(x)) > pari("10^-%d" % ACC) * (1 + abs(x)):
            ok = False
        q = pari.bestappr(pari.real(x), pari("10^%d" % (DPS // 4)))
        if abs(x - q) > pari("10^-%d" % ACC) * (1 + abs(x)) or not smooth(pari.denominator(q)):
            ok = False
        coeffs.append(q)
    expr = sum((coeffs[k] * pari("x") ** k for k in range(4)), pari(0))
    return (str(expr) if ok else None), expr


if __name__ == "__main__":
    d = int(sys.argv[1])
    L = json.load(open(f"twin_lattices_{d}.json"))
    Oz = json.load(open(f"orbit_{d}.json"))
    names = sys.argv[2].split(",") if len(sys.argv) > 2 else sorted(L)
    fn = f"invariants_F_{d}.json"
    out = json.load(open(fn)) if os.path.exists(fn) else {}
    for name in names:
        rows = member_values(d, L[name])
        Qz = pari(Oz[name]["minpoly_z"])
        P1, a1 = pari.polredbest(Qz, 1); F, a2 = pari.polredabs(P1, 1)
        phi = pari.lift(pari.Mod(pari.subst(pari.lift(a1), "x", pari.lift(a2)), F))
        rootsF = pari.polroots(F, precision=4 * DPS)
        th = []
        for r in rows:
            cand = sorted(rootsF, key=lambda t: abs(pari.subst(phi, "x", t) - r["z"]) / (1 + abs(r["z"])))
            assert abs(pari.subst(phi, "x", cand[0]) - r["z"]) / (1 + abs(r["z"])) < pari("10^-%d" % ACC)
            th.append(cand[0])
        assert len({str(pari.precision(t, 20)) for t in th}) == 4, "orbit members must hit four distinct conjugates"
        rec = {"F": str(F), "z_in_theta": str(phi), "n_involutions": [r["n_involutions"] for r in rows]}
        keys = ["J2^3J4/J10", "J2^2J6/J10"] + (["u", "v"] if all("u" in r for r in rows) else [])
        for key in keys:
            e, _ = recognise(F, th, [r[key] for r in rows])
            rec[key] = e
        rec["J2^5/J10"] = str(phi)
        out[name] = rec
        print(name, "F", F, "| recognised:", {k: (rec[k] is not None) for k in keys}, flush=True)
        json.dump(out, open(fn, "w"), indent=1)
