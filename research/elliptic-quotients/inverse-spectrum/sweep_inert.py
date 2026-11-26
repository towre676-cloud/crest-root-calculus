"""Inert-level inverse test: within each genus of nonscalar CM self-gluing Humbert forms, are inequivalent
forms separated by the quotient spectrum?  Separation is proved by an exact witness degree m (first m with
different numbers of primitive representations of m^2, i.e. of degree-m quotient classes); a pair with no
witness up to the Sturm bound B would be square-equivalent (proved by the Shimura lift + Sturm).
"""
import itertools
import json
import math
import sys
import time

import cypari2

from discform import isomorphic as disc_iso
from gluing import ns_and_q, reduce_form, isometric
from inverse import level, sturm, mat, is_scalar

pari = cypari2.Pari()
pari.allocatemem(8 * 10 ** 9)


def prim_counts(q, B):
    """primitive representation numbers of m^2 (mod +-), m = 1..B, by Moebius from qfrep."""
    r = pari.qfrep(mat(q), B * B)
    tot = [0] + [int(r[m * m - 1]) for m in range(1, B + 1)]
    prim = [0] * (B + 1)
    for m in range(1, B + 1):
        prim[m] = tot[m] - sum(prim[m // k] for k in range(2, m + 1) if m % k == 0)
    return prim[1:]


def classes_at(d, l):
    buckets = {}
    n_psi = n_prod = 0
    for a, b, c, e in itertools.product(range(l), repeat=4):
        if (a * e - b * c + 1) % l:
            continue
        psi = [[a, b], [c, e]]
        if is_scalar(d, l, psi):
            continue
        r = ns_and_q(d, l, psi)
        if r is None:
            continue
        n_psi += 1
        q = reduce_form(r["q"])
        if int(pari.qfminim(mat(q), 1, 1)[0]) > 0:
            n_prod += 1
            continue
        key = tuple(int(x) for x in pari.qfrep(mat(q), 60))
        lst = buckets.setdefault(key, [])
        for cl in lst:
            if isometric(cl["q"], q):
                cl["n"] += 1
                break
        else:
            lst.append({"q": q, "n": 1, "psi": psi})
    cls = [c for v in buckets.values() for c in v]
    return cls, n_psi, n_prod


def study(d, l):
    t = time.time()
    cls, n_psi, n_prod = classes_at(d, l)
    genera = []
    for c in cls:
        for g in genera:
            if disc_iso(g[0]["q"], c["q"]):
                g.append(c)
                break
        else:
            genera.append([c])
    res = {"d": d, "l": l, "nonscalar_glued": n_psi, "products": n_prod, "isometry_classes": len(cls), "genera": []}
    for g in genera:
        N = level(g[0]["q"])
        B = sturm(N)
        pc = [prim_counts(c["q"], B) for c in g]
        wit, unresolved = [], []
        for i, j in itertools.combinations(range(len(g)), 2):
            m = next((k + 1 for k in range(B) if pc[i][k] != pc[j][k]), None)
            if m is None:
                unresolved.append([g[i]["q"], g[j]["q"]])
            else:
                wit.append({"pair": [g[i]["q"], g[j]["q"]], "degree": m, "counts": [pc[i][m - 1], pc[j][m - 1]]})
        res["genera"].append({"size": len(g), "level_2q": N, "sturm_bound": B,
                              "forms": [{"q": c["q"], "gluings": c["n"], "psi": c["psi"]} for c in g],
                              "max_witness_degree": max((w["degree"] for w in wit), default=None),
                              "witnesses": wit, "square_equivalent_pairs": unresolved})
    res["seconds"] = round(time.time() - t, 1)
    return res


if __name__ == "__main__":
    out = []
    for spec in sys.argv[1:]:
        d, l = map(int, spec.split(":"))
        r = study(d, l)
        out.append(r)
        gs = [(g["size"], g["max_witness_degree"], len(g["square_equivalent_pairs"])) for g in r["genera"]]
        print(f"RES d={d} l={l} glued={r['nonscalar_glued']} products={r['products']} classes={r['isometry_classes']} "
              f"genera(size,maxwit,sq-eq)={gs} {r['seconds']}s", flush=True)
        json.dump(out, open(f"inert_{sys.argv[1].replace(':','_')}.json", "w"), indent=1)
