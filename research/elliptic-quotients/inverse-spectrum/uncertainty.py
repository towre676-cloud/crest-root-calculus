"""Uncertainty certificate for coprime-degree injectivity.

Coefficient vectors c (combinations sum c_x theta_x) carry the inner product <c, c'> = sum c_x c'_x / m_x,
m_x = 1 / |O(L_x)|; the transposed neighbour matrices are self-adjoint for it (m_i A_ij = m_j A_ji).
Let V be the coprime square kernel of a genus (computed with m <= B: an upper bound, and the diagonal of
the projection only decreases on a subspace) and rho_V(x) = m_x <P_V e_x, e_x> in [0, 1], sum_x rho_V(x) = dim V.
If e_i - e_j lies in V then sqrt(rho_i/m_i) + sqrt(rho_j/m_j) >= sqrt(1/m_i + 1/m_j), forcing
max(rho_i, rho_j) >= 1/2.  Hence max_x rho_V(x) < 1/2 certifies that no two classes of the genus share
their coprime spectrum.  The symmetry check m_i A_ij = m_j A_ji is asserted."""
import json, math, sys
import cypari2
from twins import closure_and_matrices
from discform import isomorphic as disc_iso
from inverse import mat
from fractions import Fraction as Fr

pari = cypari2.Pari()


def run(d, l, primes, B=300):
    classes, ng, A = closure_and_matrices(d, l, primes)
    n = len(classes)
    aut = [int(pari.qfauto(mat(q))[0]) for q in classes]
    for p in primes:
        for i in range(n):
            for j in range(n):
                assert Fr(int(A[p][i, j]), aut[i]) == Fr(int(A[p][j, i]), aut[j]), "mass symmetry fails"
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
    out = {"d": d, "l": l, "genera": []}
    for g in range(len(reps)):
        idx = [i for i in range(n) if gen[i] == g]
        M = pari.matrix(len(ms), len(idx), [S[j][i] for i in range(len(ms)) for j in idx])
        K = pari.matker(M)
        k = int(pari.matsize(K)[1])
        if k == 0:
            out["genera"].append({"size": len(idx), "dim_V": 0, "max_rho": 0.0, "certified": True}); continue
        # Gram of V in the inner product sum c c' / m_x  (= sum c c' * aut_x)
        W = pari.matrix(len(idx), len(idx), [aut[idx[i]] if i == j else 0 for i in range(len(idx)) for j in range(len(idx))])
        G = pari.mattranspose(K) * W * K
        Gi = pari.matsolve(G, pari.matid(k))
        # P_V e_x = K Gi K^T W e_x ;  rho_x = m_x <P e_x, e_x> = (1/aut_x) * (e_x^T W K Gi K^T W e_x) = aut_x (K Gi K^T)_{xx}
        P = K * Gi * pari.mattranspose(K)
        rho = [Fr(str(P[t, t])) * aut[idx[t]] for t in range(len(idx))]
        assert sum(rho) == k
        mr = max(rho)
        out["genera"].append({"size": len(idx), "dim_V": k, "max_rho": float(mr), "max_rho_exact": str(mr),
                              "certified": mr < Fr(1, 2)})
    return out


if __name__ == "__main__":
    res = []
    for spec in sys.argv[1:]:
        d, l, ps = spec.split(":")
        r = run(int(d), int(l), [int(x) for x in ps.split(",")])
        res.append(r)
        print(r["d"], r["l"], [(g["size"], g["dim_V"], round(g["max_rho"], 4), g["certified"]) for g in r["genera"]], flush=True)
    json.dump(res, open("uncertainty.json", "w"), indent=1)
