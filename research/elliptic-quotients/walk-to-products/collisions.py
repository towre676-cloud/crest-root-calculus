"""Locate the coprime-spectrum collisions found by the sweep and decide whether gluing classes are involved.
All gluing classes are enumerated (every admissible nonscalar psi); collisions = pairs with observed fraction 0
(equivalently, identical vectors in every Hecke-cyclic coordinate of 1_R)."""
import itertools, json, sys
sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
import cypari2
from gluing import ns_and_q, reduce_form
from inverse import mat, is_scalar
from hecke import neighbours, classify
from discform import isomorphic as disc_iso

pari = cypari2.Pari()
pari.allocatemem(8 * 10 ** 9)


def run(d, l, primes):
    classes, buckets = [], {}
    glue = set()
    for a, b, c, e in itertools.product(range(l), repeat=4):
        if (a * e - b * c + 1) % l or is_scalar(d, l, [[a, b], [c, e]]):
            continue
        r = ns_and_q(d, l, [[a, b], [c, e]])
        if r is None:
            continue
        q = reduce_form(r["q"])
        glue.add(classify(q, classes, buckets))
    rows = {p: {} for p in primes}
    changed = True
    while changed:
        changed = False
        for p in primes:
            for i in range(len(classes)):
                if i not in rows[p]:
                    rows[p][i] = [classify(q, classes, buckets) for q in neighbours(classes[i], p)]
                    changed = True
    n = len(classes)
    A = {p: pari.matrix(n, n, [sum(1 for j in rows[p][i] if j == k) for i in range(n) for k in range(n)]) for p in primes}
    r1 = [int(pari.qfrep(mat(q), 1)[0]) for q in classes]
    v = pari.matrix(n, 1, [1 if r1[x] else 0 for x in range(n)])
    K = [v]
    M = v
    rk = 1
    while True:
        B = pari.matimage(pari.matconcat([M] + [A[p] * M for p in primes]))
        r2 = int(pari.matsize(B)[1])
        if r2 == rk:
            break
        M, rk = B, r2
    # Hecke-word values: rows of M (basis of the cyclic module) are the observation vectors of each class
    obs = [tuple(str(M[x, k]) for k in range(rk)) for x in range(n)]
    gen = [0 if disc_iso(classes[0], q) else 1 for q in classes]
    groups = {}
    for x in range(n):
        if not r1[x]:
            groups.setdefault((gen[x], obs[x]), []).append(x)
    coll = [g for g in groups.values() if len(g) > 1]
    return {"d": d, "l": l, "closure": n, "gluing_classes": len(glue), "product_classes": sum(1 for x in r1 if x),
            "collision_groups": [[{"q": classes[x], "gluing": x in glue,
                                   "aut": int(pari.qfauto(mat(classes[x]))[0])} for x in g] for g in coll],
            "collisions_involving_gluing": sum(1 for g in coll if any(x in glue for x in g))}


if __name__ == "__main__":
    out = []
    for spec in sys.argv[1:]:
        d, l, ps = spec.split(":")
        r = run(int(d), int(l), [int(x) for x in ps.split(",")])
        out.append(r)
        print(f"d={d} l={l} closure={r['closure']} gluing={r['gluing_classes']} products={r['product_classes']} "
              f"collision groups={len(r['collision_groups'])} sizes={[len(g) for g in r['collision_groups']]} "
              f"involving gluing={r['collisions_involving_gluing']}", flush=True)
    json.dump(out, open("collisions.json", "w"), indent=1)
