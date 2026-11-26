"""Exact bielliptic models over the quartic field of moduli and their lattice identification, for any level-three D.

Generalises verify_d43.py.  For every bielliptic record in invariants_F_<D>.json the quartic-field sextic is built by
the exact descent formula descend.from_dihedral(u, v), its three absolute Igusa invariants are checked exactly against
the record, the two quotient j-invariants are checked exactly to be roots of the ring class polynomial H_{9D}, and at
every embedding of F the Arb gluing certificate of certify_link.certify_gluing gives the refined Humbert form, which
must be isometric to exactly the intended archive form.

Usage: python verify_bielliptic.py D [--write]
"""
import json
import sys

from flint import acb, ctx

from certify_link import HERE, PARI, certify_gluing, checked_roots, eval_L, gram, ideal_candidates, read, require
from descend import elliptic_j_pair, from_dihedral, multiplication_polynomial, poly
from igusa_exact import absolute_invariants

PRECISION = 3000
PARI.allocatemem(256_000_000)
KEYS = ("J2^5/J10", "J2^3J4/J10", "J2^2J6/J10")


def build(D):
    moduli = read(HERE / f"invariants_F_{D}.json")
    names = sorted(n for n, r in moduli.items() if r.get("u") and r.get("v"))
    fields = {moduli[n]["F"] for n in moduli}
    require(len(fields) == 1, "All records must share one quartic field.")
    Fs = fields.pop()
    F = poly(Fs)
    require(F.is_irreducible and F.degree() == 4, "The field polynomial must be an irreducible quartic.")
    H = list(map(int, PARI.Vec(PARI.polclass(9 * D))))
    mod = PARI(Fs).subst("x", "w")
    models = {}
    for name in names:
        r = moduli[name]
        u, v = (poly(r[k]).rem(F) for k in ("u", "v"))
        coefficients = from_dihedral(u, v, F)
        disc = (v * v - 4 * u + 18 * v - 27).rem(F)
        require(bool(disc), "The even sextic must be smooth.")
        js = elliptic_j_pair(u, v, disc, F, H)
        values = [PARI.Mod(PARI(str(c.as_expr()).replace("**", "^")).subst("x", "w"), mod) for c in coefficients]
        for key, value in zip(KEYS, absolute_invariants(PARI, values), strict=True):
            require(value == PARI.Mod(PARI(r[key]).subst("x", "w"), mod), f"{name}: invariant {key} differs.")
        z = poly(r["J2^5/J10"]).rem(F)
        require(multiplication_polynomial(z, F).degree() == 4, "The first invariant must generate F.")
        models[name] = {"F": Fs, "u": r["u"], "v": r["v"],
                        "model_constant_first": [str(c.as_expr()).replace("**", "^") for c in coefficients],
                        "quotient_j_sum_and_product": js}
    return Fs, H, models


def certify(D, Fs, H, models):
    archive = read(HERE / f"twin_lattices_{D}.json")
    F = poly(Fs)
    out = {n: dict(m, embeddings=[]) for n, m in models.items()}
    with ctx.workprec(PRECISION):
        candidates = ideal_candidates(H, D)
        for th in checked_roots([int(F.nth(i)) for i in range(5)]):
            for name, m in models.items():
                u, v = (eval_L(m[k], th) for k in ("u", "v"))
                s1 = ((u + (u * u - 4 * v ** 3).sqrt()) / 2) ** (acb(1) / 3)
                require(s1.is_finite() and not s1.contains(0), "Could not enclose a nonzero dihedral parameter.")
                glued = certify_gluing(s1, v / s1, candidates, gram(archive[name]["q"]), D)
                form = gram(glued["humbert_form"])
                matches = sorted(o for o in archive if PARI.qfisom(form, gram(archive[o]["q"])))
                require(matches == [name], f"{name} matches {matches}.")
                out[name]["embeddings"].append(glued)
    for name in out:
        out[name]["embeddings"].sort(key=lambda row: json.dumps(row, sort_keys=True))
    pairs = {}
    for name in out:
        p, side = name.split("_")
        other = f"{p}_{'q2' if side == 'q1' else 'q1'}"
        if side == "q1" and other in out:
            counts = {}
            for n in (name, other):
                q = gram(archive[n]["q"]); reps = PARI.qfrep(q, 400); prim = []
                for m in range(1, 21):
                    prim.append(2 * int(reps[m * m - 1]) - sum(prim[j - 1] for j in range(1, m) if m % j == 0))
                counts[n] = prim
            a, b = counts[name], counts[other]
            first = next(m + 1 for m in range(20) if a[m] != b[m])
            pairs[p] = {"first_degree": first, "primitive_counts_to_first": {n: c[:first] for n, c in counts.items()}}
    return {"D": D, "field": Fs, "ring_class_polynomial": H, "models": out, "pairs": pairs}


if __name__ == "__main__":
    D = int(sys.argv[1])
    Fs, H, models = build(D)
    print("exact quartic-field bielliptic models:", sorted(models), flush=True)
    result = certify(D, Fs, H, models)
    path = HERE / f"certified_bielliptic_{D}.json"
    if "--write" in sys.argv:
        path.write_text(json.dumps(result, indent=1) + "\n", encoding="utf-8")
    else:
        require(result == read(path), "Replay differs from the committed certificate.")
    for p, r in result["pairs"].items():
        print(f"pair {p}: first distinguishing degree {r['first_degree']}, counts {r['primitive_counts_to_first']}")
    for name, r in result["models"].items():
        print(f"{name}: certified at {len(r['embeddings'])} embeddings, isometric only to its archive form.")
