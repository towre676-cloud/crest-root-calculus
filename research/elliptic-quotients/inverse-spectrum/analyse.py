"""Genus and square-equivalence analysis of the isometry classes produced by inverse.py."""
import json
import sys

import cypari2

from discform import isomorphic as disc_iso
from inverse import level, square_counts, sturm, mat, run
from gluing import isometric

pari = cypari2.Pari()


def analyse(classes):
    # genera (exact, Nikulin)
    genera = []
    for c in classes:
        for g in genera:
            if g[0]["det"] == c["det"] and disc_iso(g[0]["q"], c["q"]):
                g.append(c)
                break
        else:
            genera.append([c])
    report = []
    for g in genera:
        N = level(g[0]["q"])
        B = sturm(N)
        sq = [square_counts(c["q"], B) for c in g]
        # square-equivalence classes within the genus (proved by Sturm)
        sqcl = []
        for c, s in zip(g, sq):
            for k in sqcl:
                if k["counts"] == s:
                    k["members"].append(c["q"])
                    break
            else:
                sqcl.append({"counts": s, "members": [c["q"]]})
        full = []
        # full theta-series equality would force isometry (Schiemann); record it as a check
        report.append({"det": g[0]["det"], "level_2q": N, "sturm_bound": B, "levels_l": sorted({c["level"] for c in g}),
                       "isometry_classes": len(g), "forms": [c["q"] for c in g],
                       "gluings": [c["gluings"] for c in g],
                       "square_classes": [{"members": k["members"], "first_counts": k["counts"][:12]} for k in sqcl]})
    return report


if __name__ == "__main__":
    d = int(sys.argv[1])
    ls = [int(x) for x in sys.argv[2].split(",")]
    cl = run([d], ls)[d]
    rep = analyse(cl)
    for r in rep:
        sizes = [len(k["members"]) for k in r["square_classes"]]
        print(f"d={d} l={r['levels_l']} det={r['det']} N={r['level_2q']} B={r['sturm_bound']} "
              f"classes={r['isometry_classes']} square-classes={len(sizes)} sizes={sizes}", flush=True)
    json.dump(rep, open(f"analysis_{d}.json", "w"), indent=1)


def cross_witnesses(report, Bmax=400):
    """For forms in different genera: the first m with r(m^2) different (an exact separating witness)."""
    forms = [(gi, q) for gi, r in enumerate(report) for q in r["forms"]]
    B = Bmax
    counts = {}
    unresolved = []
    for gi, q in forms:
        counts[str(q)] = square_counts(q, B)
    seps = 0
    for (g1, q1), (g2, q2) in __import__("itertools").combinations(forms, 2):
        if g1 == g2:
            continue
        a, b = counts[str(q1)], counts[str(q2)]
        m = next((k + 1 for k in range(B) if a[k] != b[k]), None)
        if m is None:
            unresolved.append([q1, q2])
        else:
            seps += 1
    return {"pairs_separated": seps, "unresolved_to_B": unresolved, "B": B}
