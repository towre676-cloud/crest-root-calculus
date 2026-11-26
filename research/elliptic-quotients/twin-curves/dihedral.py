"""Dihedral invariants (u, v) of the bielliptic structure (Shaska-Voelklein): for y^2 = g(x^2), g = U^3 + A U^2 + B U + C,
u = A^3/C + B^3/C^2, v = A B / C (invariant under x -> t x and under reversal of g).  Test u, v in F = Q(theta)."""
import json, sys
import mpmath as mp
import cypari2
from bielliptic_model import models

pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9); pari.set_real_precision(260)
mp.mp.dps = 270
fn = sys.argv[1]
Cv = json.load(open(fn)); Mf = json.load(open(fn.replace("curves", "moduli")))
P = lambda z: pari("%s + I*(%s)" % (mp.nstr(z.real, 250), mp.nstr(z.imag, 250)))
out = {}
for name in ("q1", "q2"):
    ros = [mp.mpc(mp.mpf(r), mp.mpf(i)) for r, i in Cv[name]["rosenhain"]]
    (alpha, beta, (A, B, C)), = models([mp.mpc(0), mp.mpc(1)] + ros)
    u = A ** 3 / C + B ** 3 / C ** 2
    v = A * B / C
    Q = pari(Mf[name]["F"])
    a = pari(Mf[name]["z_in_theta"])
    z = pari(Cv[name]["abs_igusa"][0])
    theta = next(r for r in pari.polroots(Q, precision=900) if abs(pari.subst(a, "x", r) - z) < pari("10^-150") * (1 + abs(z)))
    out[name] = {}
    for key, t in (("u", P(u)), ("v", P(v))):
        basis = [t, pari(1), theta, theta ** 2, theta ** 3]
        vec = pari.lindep(basis)
        res = abs(sum((vec[i] * basis[i] for i in range(5)), pari(0)))
        h = max(abs(int(x)) for x in vec)
        ok = vec[0] != 0 and float(res) < 1e-200 and len(str(h)) < 60
        expr = str(pari("-(%s)/(%s)" % (" + ".join(f"({int(vec[i])})*x^{i-1}" for i in range(1, 5)), int(vec[0])))) if ok else None
        out[name][key] = {"in_F": bool(ok), "expr_in_theta": expr, "height_digits": len(str(h)), "residual": float(res)}
        print(name, key, "in F:", ok, "height 10^%d" % len(str(h)), "residual %.1e" % float(res), "\n   ", expr, flush=True)
json.dump(out, open(fn.replace("curves", "dihedral"), "w"), indent=1)
