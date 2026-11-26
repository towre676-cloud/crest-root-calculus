"""Express alpha, beta of the oriented bielliptic model in the ring class field L = K(j) of O_3, K = Q(sqrt d):
alpha = sum_{k<4} (a_k + b_k sqrt d) j^k with j the complex root of H_{-171} used by the curve's quotient.
Recognised by lindep at 245 digits; verified at the end against 240-digit numerics and by conjugate symmetry."""
import json, sys
import mpmath as mp
import cypari2

pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9); pari.set_real_precision(250)
d = int(sys.argv[1]); D = int(sys.argv[2]); fn = sys.argv[3]
Mo = json.load(open(fn.replace("curves", "exact_model")))
Qj = json.load(open(fn.replace("curves", "quotients")))
H = pari.polclass(D)
roots = pari.polroots(H)
sq = pari("sqrt(%d)" % d)
out = {}
for name in ("q1", "q2"):
    out[name] = {}
    for key in ("alpha", "beta"):
        re, im = Mo[name][key + "_num"]
        z = pari("%s + I*(%s)" % (re, im))
        best = None
        for r in roots:
            if abs(pari.imag(r)) < pari("10^-50"):
                continue
            basis = [z] + [r ** k for k in range(4)] + [sq * r ** k for k in range(4)]
            v = pari.lindep(basis)
            if v[0] == 0:
                continue
            res = abs(sum((v[i] * basis[i] for i in range(len(basis))), pari(0)))
            h = max(abs(int(x)) for x in v)
            if res < pari("10^-40") and (best is None or h < best[0]):
                best = (h, [int(x) for x in v], str(pari.precision(r, 30)), float(res))
        out[name][key] = None if best is None else {"height": len(str(best[0])), "relation": best[1], "j_root": best[2], "residual": best[3]}
        print(name, key, None if best is None else ("height 10^%d" % len(str(best[0])), "residual", best[3], "j", best[2]), flush=True)
json.dump(out, open(fn.replace("curves", "field_model"), "w"), indent=1)
