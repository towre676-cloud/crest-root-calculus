"""Exact p-power observability of the gluing Humbert forms at one level.

Verified relation on the p-neighbour closure (eichler_check2.py; correction from the parallel branch):
    r(p^2 n) + (D n / p) r(n) + p r(n/p^2) = (A_p r)(n),
(the earlier middle term -(-n/p) agrees with it only when (-D/p) = -1; the Krylov space is the same either way)
so v_k = (r_j(p^{2k}))_j = P_k(A_p) v_0 with v_0 = (r_j(1))_j and monic P_k of degree k.  Hence the degrees
1, p, p^2, ... separate gluing classes i, j iff (e_i - e_j)^T A_p^k v_0 != 0 for some k < n.  Exact over Z.
Reports, per prime, the number of unseparated gluing pairs and whether p splits or is inert in K."""
import itertools, json, sys
import cypari2
from twins import closure_and_matrices
from inverse import mat

pari = cypari2.Pari()


def observe(d, l, p):
    classes, ng, A = closure_and_matrices(d, l, [p])
    n = len(classes)
    Ap = A[p]
    v0 = pari.matrix(n, 1, [2 * int(pari.qfrep(mat(q), 1)[0]) for q in classes])
    K = [v0]
    for _ in range(n):
        K.append(Ap * K[-1])
    Kr = pari.matconcat(K)
    rows = [[int(Kr[i, k]) for k in range(n + 1)] for i in range(ng)]
    bad = [(i, j) for i, j in itertools.combinations(range(ng), 2) if rows[i] == rows[j]]
    kind = {1: "split", -1: "inert", 0: "ramified"}[int(pari.kronecker(d, p))]
    return {"d": d, "l": l, "p": p, "kind_in_K": kind, "closure": n, "gluing": ng,
            "krylov_rank": int(pari.matrank(Kr)), "unseparated_pairs": len(bad),
            "pairs": [[classes[i], classes[j]] for i, j in bad]}


if __name__ == "__main__":
    out = []
    d, l = int(sys.argv[1]), int(sys.argv[2])
    for p in map(int, sys.argv[3].split(",")):
        r = observe(d, l, p)
        out.append(r)
        print({k: v for k, v in r.items() if k != "pairs"}, flush=True)
    json.dump(out, open(f"observe_{d}_{l}.json", "w"), indent=1)
