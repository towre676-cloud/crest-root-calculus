"""Exact data over the field of moduli F: theta = root of polredbest(minpoly of z = J2^5/J10), embedded so that
z = a(theta).  Each quantity t is tested for t = sum_{k<4} c_k theta^k (lindep at 250 digits).  Acceptance:
residual < 10^-200 with coefficient height far below the random-relation height (~10^50 for 5 complex terms)."""
import json, sys
import cypari2

pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9); pari.set_real_precision(260)
fn = sys.argv[1]
Cv = json.load(open(fn)); Mo = json.load(open(fn.replace("curves", "exact_model")))
A = json.load(open(fn.replace("curves", "algdep")))
out = {}
for name in ("q1", "q2"):
    z = pari(Cv[name]["abs_igusa"][0])
    P = pari(A[name])
    Q, amod = pari.polredbest(P, 1)
    a = pari.lift(amod)
    theta = next(r for r in pari.polroots(Q, precision=900) if abs(pari.subst(a, "x", r) - z) < pari("10^-200") * (1 + abs(z)))
    out[name] = {"F": str(Q), "theta_numeric": str(pari.precision(theta, 30)), "z_in_theta": str(a)}
    targets = {"J2^5/J10": z, "J2^3J4/J10": pari(Cv[name]["abs_igusa"][1]), "J2^2J6/J10": pari(Cv[name]["abs_igusa"][2]),
               "alpha": pari("%s + I*(%s)" % tuple(Mo[name]["alpha_num"])),
               "beta": pari("%s + I*(%s)" % tuple(Mo[name]["beta_num"]))}
    for key, t in targets.items():
        basis = [t, pari(1), theta, theta ** 2, theta ** 3]
        v = pari.lindep(basis)
        res = abs(sum((v[i] * basis[i] for i in range(5)), pari(0)))
        h = max(abs(int(x)) for x in v)
        ok = v[0] != 0 and float(res) < 1e-200 and len(str(h)) < 45
        expr = str(pari("-(%s)/(%s)" % (" + ".join(f"({int(v[i])})*x^{i-1}" for i in range(1, 5)), int(v[0])))) if ok else None
        out[name][key] = {"in_F": bool(ok), "expr_in_theta": expr, "height_digits": len(str(h)), "residual": float(res)}
        print(name, key, "in F:", ok, "height 10^%d" % len(str(h)), "residual %.1e" % float(res), expr if ok else "", flush=True)
json.dump(out, open(fn.replace("curves", "moduli"), "w"), indent=1)
