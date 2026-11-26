"""Genus-resolved Hecke decomposition of the square kernel.

For each level: closure of the gluing classes under p-neighbours, split into genera (Nikulin test),
characteristic polynomials of T(p^2) on each genus block, and the square kernel decomposed as
  K = K_within (combinations inside one genus)  +  K_cross (the rest).
Hypothesis tested: every cuspidal Hecke eigensystem occurs in both genera of the level, the square
coefficients see only the eigensystem (Shimura lift), so each shared eigensystem contributes one
cross-genus relation; the within-genus kernel is the genuine Waldspurger vanishing.
"""
import json, sys
import cypari2
from hecke import neighbours, classify, mat
from discform import isomorphic as disc_iso
import glob

pari = cypari2.Pari()


def run(d, l, p, Bsq=300):
    forms = []
    for f in glob.glob("inert_*.json"):
        for R in json.load(open(f)):
            if (R["d"], R["l"]) == (d, l):
                for g in R["genera"]:
                    forms += [x["q"] for x in g["forms"]]
    classes, buckets = [], {}
    for q in forms:
        classify(q, classes, buckets)
    ng = len(classes)
    rows, i = {}, 0
    while i < len(classes):
        rows[i] = [classify(q, classes, buckets) for q in neighbours(classes[i], p)]
        i += 1
    n = len(classes)
    A = [[sum(1 for j in rows[i] if j == k) for k in range(n)] for i in range(n)]
    gen, reps = [], []
    for q in classes:
        for gi, r in enumerate(reps):
            if disc_iso(r, q):
                gen.append(gi); break
        else:
            reps.append(q); gen.append(len(reps) - 1)
    S = []
    for q in classes:
        r = pari.qfrep(mat(q), Bsq * Bsq)
        S.append([int(r[m * m - 1]) for m in range(1, Bsq + 1)])
    def kerdim(idx):
        if not idx:
            return 0
        M = pari.matrix(Bsq, len(idx), [S[j][i] for i in range(Bsq) for j in idx])
        return int(pari.matsize(pari.matker(M))[1])
    allidx = list(range(n))
    kall = kerdim(allidx)
    kwithin = sum(kerdim([i for i in allidx if gen[i] == g]) for g in range(len(reps)))
    glue = list(range(ng))
    kglue = kerdim(glue)
    kglue_within = sum(kerdim([i for i in glue if gen[i] == g]) for g in range(len(reps)))
    cps = []
    for g in range(len(reps)):
        idx = [i for i in allidx if gen[i] == g]
        M = pari.matrix(len(idx), len(idx), [A[i][j] for i in idx for j in idx])
        cps.append(pari.factor(pari.charpoly(M)))
    def fdict(F):
        return {str(F[k, 0]): int(F[k, 1]) for k in range(int(pari.matsize(F)[0]))}
    f0, f1 = (fdict(cps[0]), fdict(cps[1])) if len(cps) == 2 else (fdict(cps[0]), {})
    shared = sum(min(f0[k], f1.get(k, 0)) * int(pari.poldegree(pari(k))) for k in f0)
    res = {"d": d, "l": l, "p": p, "gluing_classes": ng, "closure_classes": n,
           "genus_sizes_closure": [gen.count(g) for g in range(len(reps))],
           "genus_sizes_gluing": [sum(1 for i in glue if gen[i] == g) for g in range(len(reps))],
           "charpoly_factors_by_genus": [f0, f1],
           "shared_eigenvalue_dim": shared,
           "K_closure": kall, "K_closure_within_genus": kwithin, "K_closure_cross": kall - kwithin,
           "K_gluing": kglue, "K_gluing_within_genus": kglue_within, "K_gluing_cross": kglue - kglue_within}
    return res


if __name__ == "__main__":
    out = []
    for spec in sys.argv[1:]:
        d, l, p = map(int, spec.split(":"))
        r = run(d, l, p)
        out.append(r)
        print({k: v for k, v in r.items() if k != "charpoly_factors_by_genus"}, flush=True)
    json.dump(out, open("hecke2_results.json", "w"), indent=1)
