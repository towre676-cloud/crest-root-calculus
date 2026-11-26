"""Archived construction of quartic models for bielliptic twins.

The ball computation below uses a supplied list of CM ideal classes.
Use verify_d43.py to prove the list is exhaustive for p3 before matching
the lattice. The D=-19 equations are independently checked by certify_link.py.

1. Exact: from (u, v) in F, the absolute Igusa invariants (dihedral_igusa.json polynomials) equal the recognised
   invariants in invariants_F_<d>.json (equality in F).
2. Exact model over F (descend.py's formula):  Y^2 = a X^6 + b X^5 + c X^4 + d X^3 + Delta c X^2 + Delta^2 b X + Delta^3 a,
   Delta = u^2 - 4 v^3, a = u - v^2, b = 3 Delta + 3u^2 - 2u v^2, c = Delta (15u + v^2), d = 2 Delta (5 Delta + 5u^2 + 2u v^2).
   It is geometrically isomorphic to y^2 = x^6 - s1 x^4 + s2 x^2 - 1 whenever s1^3 + s2^3 = u, s1 s2 = v (identity of descend.py).
3. Certificate (Arb, at every embedding of F): u, v are enclosed, s1 = cube root of (u + sqrt(u^2 - 4 v^3))/2, s2 = v/s1,
   and certify_link_any's certified gluing computation returns the refined Humbert form, compared by qfisom with the twins.
"""
import json, sys
sys.path.insert(0, ".")
import cypari2
from flint import acb, arb, acb_poly, ctx

pari = cypari2.Pari(); pari.allocatemem(4 * 10 ** 9)
d = int(sys.argv[1])
N = json.load(open(f"invariants_F_{d}.json"))
Jp = [pari(s) for s in json.load(open("dihedral_igusa.json"))["J2,J4,J6,J8,J10 in (u,v)"]]
sys.argv = [sys.argv[0], str(d)]
import certify_link_any as C
from polarised import ideals
from gluing import isometric
T = json.load(open(f"twin_lattices_{d}.json"))
out = {}
for name, r in sorted(N.items()):
    if not r.get("u"):
        continue
    F = pari(r["F"])
    md = lambda e: pari.Mod(pari(e), F)
    u, v = md(r["u"]), md(r["v"])
    J2, J4, J6, J8, J10 = [pari.substvec(p, ["u", "v"], [u, v]) for p in Jp]
    mine = [J2 ** 5 / J10, J2 ** 3 * J4 / J10, J2 ** 2 * J6 / J10]
    rec = [md(r["J2^5/J10"]), md(r["J2^3J4/J10"]), md(r["J2^2J6/J10"])]
    exact_ok = all(pari.lift(a - b) == 0 for a, b in zip(mine, rec))
    De = u ** 2 - 4 * v ** 3
    A_ = u - v ** 2; B_ = 3 * De + 3 * u ** 2 - 2 * u * v ** 2; C_ = De * (15 * u + v ** 2)
    D_ = 2 * De * (5 * De + 5 * u ** 2 + 2 * u * v ** 2)
    coeffs = [De ** 3 * A_, De ** 2 * B_, De * C_, D_, C_, B_, A_]          # constant first
    model = [str(pari.lift(c)) for c in coeffs]
    # Arb certificate at every embedding of F
    ctx.prec = 3000
    Fc = [int(pari.polcoef(F, k)) for k in range(5)]
    roots = acb_poly([acb(c) for c in Fc]).roots(tol=arb(2) ** -2900, maxprec=4 * ctx.prec)
    cands = [(f, a) for f in (1, 3) for a in ideals(d, f)]
    ev = lambda e, th: C.eval_L(e, th)
    certs = []
    for k, th in enumerate(roots):
        uu, vv = ev(r["u"], th), ev(r["v"], th)
        t = (uu + (uu * uu - 4 * vv ** 3).sqrt()) / 2
        s1 = t ** (acb(1) / 3)
        s2 = vv / s1
        assert (s1 ** 3 + s2 ** 3).overlaps(uu) and (s1 * s2).overlaps(vv)
        orig = C.eval_L
        C.eval_L = lambda expr, eta: expr          # pass the enclosures straight through
        try:
            res = C.certify(name, s1, s2, None, cands)
        finally:
            C.eval_L = orig
        iso = [tn for tn in sorted(T) if isometric(res["humbert_form"], T[tn]["q"])]
        certs.append({"embedding": k, "humbert_form": res["humbert_form"], "isometric_to": iso,
                      "E1_conductor": res["E1_class"]["conductor"], "E2_conductor": res["E2_class"]["conductor"]})
        print(name, "embedding", k, "Humbert form", res["humbert_form"], "isometric to", iso, flush=True)
    out[name] = {"F": r["F"], "u": r["u"], "v": r["v"], "invariants_from_uv_equal_recognised": bool(exact_ok),
                 "model_constant_first": model, "certificates": certs,
                 "certified": all(c["isometric_to"] == [name] for c in certs)}
    print(name, "| exact (u,v) -> invariants check:", exact_ok, "| certified at all", len(certs), "embeddings:", out[name]["certified"], flush=True)
json.dump(out, open(f"models_certified_{d}.json", "w"), indent=1)
