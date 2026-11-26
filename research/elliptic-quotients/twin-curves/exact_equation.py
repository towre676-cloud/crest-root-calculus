"""Explicit equations for the twins with exact verification.

F = Q(t), t^4 - t^3 - 3t^2 + 2t + 4 = 0 (field of moduli); exact dihedral invariants u, v in F.
Model: y^2 = x^6 - s1 x^4 + s2 x^2 - 1 with s1 a root of S^6 - u S^3 + v^3 over F and s2 = v/s1.
The sextic is factored over F, a root adjoined (rnfequation), and the absolute Igusa invariants of the model are
computed EXACTLY in that number field and compared with the exact F-invariants recognised from the period matrix,
under every embedding of F into the model field."""
import json, sys
import cypari2

pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9)
fn = sys.argv[1]
Mf = json.load(open(fn.replace("curves", "moduli"))); Dh = json.load(open(fn.replace("curves", "dihedral")))
out = {}
for name in ("q1", "q2"):
    Ft = pari(Mf[name]["F"]).subst("x", "t")
    nf = pari.nfinit(Ft)
    el = lambda s: pari.Mod(pari(s).subst("x", "t"), Ft)
    u, v = el(Dh[name]["u"]["expr_in_theta"]), el(Dh[name]["v"]["expr_in_theta"])
    sext = pari("x^6") - u * pari("x^3") + v ** 3
    fac = pari.nffactor(nf, pari.lift(sext))
    degs = [int(pari.poldegree(fac[k, 0], "x")) for k in range(int(pari.matsize(fac)[0]))]
    h = fac[0, 0]
    z, a, k = pari.rnfequation(nf, h, 1)
    a = pari.Mod(pari.lift(a), z)
    s1 = pari.Mod(pari("x"), z) - int(k) * a
    tz = a                                                   # image of t in Q[x]/z
    uz = pari.Mod(pari.subst(pari.lift(u), "t", pari.lift(tz)), z)
    vz = pari.Mod(pari.subst(pari.lift(v), "t", pari.lift(tz)), z)
    assert pari.lift(s1 ** 6 - uz * s1 ** 3 + vz ** 3) == 0
    s2 = vz / s1
    GEN = pari("genus2igusa(x^6 - a*x^4 + b*x^2 - 1)")          # generic invariants, exact polynomials in a, b
    J2, J4, J6, J8, J10 = [pari.substvec(g, ["a", "b"], [s1, s2]) for g in GEN]
    model_abs = [J2 ** 5 / J10, J2 ** 3 * J4 / J10, J2 ** 2 * J6 / J10]
    target = [pari.Mod(pari.subst(pari(Mf[name][kk]["expr_in_theta"]), "x", pari.lift(tz)), z)
              for kk in ("J2^5/J10", "J2^3J4/J10", "J2^2J6/J10")]
    match = all(pari.lift(m - t) == 0 for m, t in zip(model_abs, target))
    red = pari.polredbest(z)
    out[name] = {"F": str(Ft), "u": str(pari.lift(u)), "v": str(pari.lift(v)),
                 "sextic_factor_degrees_over_F": degs, "model_field": str(z), "model_field_degree": int(pari.poldegree(z)),
                 "model_field_polredbest": str(red), "s1": str(pari.lift(s1)), "s2": str(pari.lift(s2)),
                 "igusa_exact_match": bool(match)}
    print(name, "S^6 - uS^3 + v^3 factors over F with degrees", degs, "| model field degree", int(pari.poldegree(z)),
          "| polredbest", red, "| exact Igusa match:", match, flush=True)
json.dump(out, open(fn.replace("curves", "equations"), "w"), indent=1)
