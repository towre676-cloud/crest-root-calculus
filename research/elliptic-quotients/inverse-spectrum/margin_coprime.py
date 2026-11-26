"""Hecke margin for coprime-degree injectivity.  Per genus of the p-neighbour closure:
  k   = dim of the kernel of the coprime square-count matrix (m <= B, gcd(m, 2 l |D|) = 1), an upper bound,
  h   = min over gluing pairs (i, j) of dim Z[T(p^2), T(p'^2)] (e_i - e_j)  (exact),
and the margin h - k.  h > k certifies coprime-degree separation of every pair with no witness degree,
because the coprime square kernel is stable under T(p^2) for good p."""
import itertools, json, math, sys
import cypari2
from twins import closure_and_matrices, krylov_dim
from discform import isomorphic as disc_iso
from inverse import mat

pari = cypari2.Pari()


def run(d, l, primes, B=300):
    classes, ng, A = closure_and_matrices(d, l, primes)
    n = len(classes)
    gen, reps = [], []
    for q in classes:
        for gi, r in enumerate(reps):
            if disc_iso(r, q):
                gen.append(gi); break
        else:
            reps.append(q); gen.append(len(reps) - 1)
    bad = 2 * l * abs(d)
    ms = [m for m in range(1, B + 1) if math.gcd(m, bad) == 1]
    S = []
    for q in classes:
        r = pari.qfrep(mat(q), B * B)
        S.append([int(r[m * m - 1]) for m in ms])
    ops = [pari.mattranspose(A[p]) for p in primes]
    out = {"d": d, "l": l, "primes": primes, "closure": n, "gluing": ng, "coprime_degrees_used": len(ms), "genera": []}
    for g in range(len(reps)):
        idx = [i for i in range(n) if gen[i] == g]
        M = pari.matrix(len(ms), len(idx), [S[j][i] for i in range(len(ms)) for j in idx])
        k = int(pari.matsize(pari.matker(M))[1])
        dims = []
        for i, j in itertools.combinations([x for x in idx if x < ng], 2):
            v = pari.matrix(n, 1, [1 if t == i else (-1 if t == j else 0) for t in range(n)])
            dims.append(krylov_dim([v], ops))
        h = min(dims) if dims else None
        out["genera"].append({"closure_size": len(idx), "gluing_in_genus": sum(1 for x in idx if x < ng),
                              "coprime_kernel": k, "min_pair_module": h,
                              "margin": (h - k) if h is not None else None})
    return out


if __name__ == "__main__":
    res = []
    for spec in sys.argv[1:]:
        d, l, ps = spec.split(":")
        r = run(int(d), int(l), [int(x) for x in ps.split(",")])
        res.append(r)
        print(r["d"], r["l"], "closure", r["closure"], [(g["closure_size"], g["coprime_kernel"], g["min_pair_module"], g["margin"]) for g in r["genera"]], flush=True)
    json.dump(res, open("margin_coprime.json", "w"), indent=1)
