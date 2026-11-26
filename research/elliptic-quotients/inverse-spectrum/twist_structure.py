"""Joint Hecke eigensystems on the p-neighbour closure, their observability, and their twist structure.

The A_p (p good) commute; a generic combination is diagonalised (numerically, then each eigenvalue is
confirmed as an exact integer root of the exact characteristic polynomial of A_p, or the system is
reported as algebraic).  A system chi is observable when the projection of v_0 = (r_j(1)) on its
eigenspace is nonzero.  Two observable systems are coincidence-dangerous at p when chi(p) = chi'(p).
Rajan (IMRN 1998): coincidence on a set of positive upper density forces twist-equivalence, which for
these weight-2 eigensystems we test by |chi(p)| = |chi'(p)| at all sampled p (necessary condition).
"""
import itertools, json, sys
import numpy as np
import cypari2
from twins import closure_and_matrices
from inverse import mat

pari = cypari2.Pari()


def run(d, l, primes):
    classes, ng, A = closure_and_matrices(d, l, primes)
    n = len(classes)
    Ms = {p: np.array([[int(A[p][i, j]) for j in range(n)] for i in range(n)], dtype=float) for p in primes}
    rng = np.random.default_rng(1)
    C = sum(rng.normal() * Ms[p] for p in primes)
    w, V = np.linalg.eig(C)
    V = np.real_if_close(V)
    systems = []
    for k in range(n):
        v = V[:, k]
        lam = {p: (Ms[p] @ v)[np.argmax(np.abs(v))] / v[np.argmax(np.abs(v))] for p in primes}
        lam = {p: complex(x).real for p, x in lam.items()}
        systems.append((tuple(round(lam[p], 6) for p in primes), k))
    v0 = np.array([2 * int(pari.qfrep(mat(q), 1)[0]) for q in classes], dtype=float)
    # observability: coordinates of v0 in eigenbasis
    coeff = np.linalg.solve(V, v0)
    groups = {}
    for (key, k) in systems:
        groups.setdefault(key, []).append(k)
    obs = []
    for key, ks in groups.items():
        # projection of v0 onto the eigenspace
        proj = sum(coeff[k] * V[:, k] for k in ks)
        if np.linalg.norm(proj[:ng]) > 1e-8:
            obs.append({"eigs": dict(zip(primes, key)), "mult": len(ks)})
    eis = [s for s in obs if all(abs(s["eigs"][p] - (p + 1)) < 1e-6 for p in primes)]
    cusp = [s for s in obs if s not in eis]
    twist_pairs, sporadic = [], []
    for a, b in itertools.combinations(cusp, 2):
        same_abs = all(abs(abs(a["eigs"][p]) - abs(b["eigs"][p])) < 1e-6 for p in primes)
        coinc = [p for p in primes if abs(a["eigs"][p] - b["eigs"][p]) < 1e-6]
        if same_abs:
            twist_pairs.append({"a": a["eigs"], "b": b["eigs"], "coincide_at": coinc})
        elif coinc:
            sporadic.append({"a": a["eigs"], "b": b["eigs"], "coincide_at": coinc})
    bad_primes = sorted({p for t in twist_pairs + sporadic for p in t["coincide_at"]})
    return {"d": d, "l": l, "primes": primes, "closure": n, "observable_systems": len(obs),
            "observable_cusp_systems": len(cusp), "cusp_eigenvalues": [s["eigs"] for s in cusp],
            "twist_equivalent_pairs": twist_pairs, "sporadic_coincidences": sporadic,
            "dangerous_primes": bad_primes, "safe_primes": [p for p in primes if p not in bad_primes]}


LEVELS = {(-3, 5): [7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
          (-4, 7): [3, 5, 11, 13, 17, 19, 23, 29, 31, 37, 41],
          (-7, 5): [3, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43],
          (-8, 5): [3, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41],
          (-8, 7): [3, 5, 11, 13, 17, 19, 23, 29, 31, 37, 41],
          (-11, 7): [3, 5, 13, 17, 19, 23, 29, 31, 37, 41, 43],
          (-3, 11): [7, 13, 17, 19, 23, 29, 31, 37, 41, 43, 47],
          (-4, 11): [3, 5, 7, 13, 17, 19, 23, 29, 31, 37, 41]}

if __name__ == "__main__":
    out = []
    for (d, l), ps in LEVELS.items():
        r = run(d, l, ps)
        out.append(r)
        print(f"d={d} l={l} closure={r['closure']} observable cusp systems={r['observable_cusp_systems']} "
              f"twist pairs={len(r['twist_equivalent_pairs'])} sporadic={len(r['sporadic_coincidences'])} "
              f"dangerous={r['dangerous_primes']} safe={r['safe_primes']}", flush=True)
        for t in r["twist_equivalent_pairs"]:
            print("   twist pair:", t, flush=True)
    json.dump(out, open("twist_structure.json", "w"), indent=1)
