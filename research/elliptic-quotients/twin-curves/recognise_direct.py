"""One-pass variant of orbit.py + recognise_F.py: orbit values, minimal polynomial of z, and recognition in F.

Writes invariants_F_<d>_<tag>.json incrementally.  Same acceptance filter as recognise_F.py."""
import json, os, sys
sys.path.insert(0, ".")
import recognise_F as R
pari, DPS = R.pari, R.DPS

d = int(sys.argv[1]); names = sys.argv[2].split(","); tag = sys.argv[3]
L = json.load(open(f"twin_lattices_{d}.json"))
fn = f"invariants_F_{d}_{tag}.json"
out = json.load(open(fn)) if os.path.exists(fn) else {}
for name in names:
    if name in out:
        continue
    rows = R.member_values(d, L[name])
    zs = R.O.distinct([r["z"] for r in rows] + [pari.conj(r["z"]) for r in rows], pari("10^-%d" % (DPS // 2)))
    Qz, worst = R.O.rational_poly(zs)
    P1, a1 = pari.polredbest(Qz, 1); F, a2 = pari.polredabs(P1, 1)
    phi = pari.lift(pari.Mod(pari.subst(pari.lift(a1), "x", pari.lift(a2)), F))
    rootsF = pari.polroots(F, precision=4 * DPS)
    th = []
    for r in rows:
        cand = sorted(rootsF, key=lambda t: abs(pari.subst(phi, "x", t) - r["z"]) / (1 + abs(r["z"])))
        assert abs(pari.subst(phi, "x", cand[0]) - r["z"]) / (1 + abs(r["z"])) < pari("10^-%d" % R.ACC)
        th.append(cand[0])
    assert len({str(pari.precision(t, 20)) for t in th}) == 4
    rec = {"F": str(F), "minpoly_z": str(Qz), "z_in_theta": str(phi), "n_involutions": [r["n_involutions"] for r in rows]}
    keys = ["J2^3J4/J10", "J2^2J6/J10"] + (["u", "v"] if all("u" in r for r in rows) else [])
    for key in keys:
        rec[key] = R.recognise(F, th, [r[key] for r in rows])[0]
    rec["J2^5/J10"] = str(phi)
    out[name] = rec
    print(name, "F", F, "| involutions", rec["n_involutions"], "| recognised:", {k: (rec[k] is not None) for k in keys}, flush=True)
    json.dump(out, open(fn, "w"), indent=1)
