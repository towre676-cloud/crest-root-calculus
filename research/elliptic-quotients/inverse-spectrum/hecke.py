"""Hecke structure of the all-degree (square) relations among gluing Humbert forms at one level.

T(p^2) acts on classes of ternary lattices through Kneser p-neighbours (p odd, p not dividing the
discriminant): for each isotropic line v mod p, lifted to q(v) = 0 mod p^2,
    L' = { x in L : B(x, v) = 0 mod p } + Z v/p.
Eichler's commutation relation gives theta_L | T(p^2) = sum_{L'} theta_{L'} + c_p theta_L, so the space of
coefficient vectors c with sum c_i r_{q_i}(m^2) = 0 for all m (both genera mixed) must be stable under
c -> A_p^T c, A_p the neighbour matrix.  We (1) compute the closure of the gluing classes under
p-neighbours, (2) the square kernel K on that closure, (3) test A_p^T K subset K, (4) decompose into
simultaneous eigenspaces of the A_p and record which are square-null.
"""
import itertools, json, sys, glob, math
import cypari2
from gluing import reduce_form, isometric
from inverse import mat

pari = cypari2.Pari()
pari.allocatemem(8 * 10 ** 9)


def qv(B, x):
    return sum(x[i] * B[i][j] * x[j] for i in range(3) for j in range(3))


def bv(B, x, y):
    return sum(x[i] * B[i][j] * y[j] for i in range(3) for j in range(3))


def neighbours(B, p):
    out = []
    pts = [(1, y, z) for y in range(p) for z in range(p)] + [(0, 1, z) for z in range(p)] + [(0, 0, 1)]
    for v in pts:
        if qv(B, v) % p:
            continue
        v = list(v)
        Bv = [sum(B[i][j] * v[j] for j in range(3)) for i in range(3)]
        # lift: v + p w with q(v) + 2 p B(v, w) = 0 mod p^2
        k = (qv(B, v) // p) % p
        i0 = next(i for i in range(3) if Bv[i] % p)
        t = (-k * pow(2 * Bv[i0], -1, p)) % p
        w = [0, 0, 0]; w[i0] = t
        v = [v[i] + p * w[i] for i in range(3)]
        assert qv(B, v) % (p * p) == 0
        Bv = [sum(B[i][j] * v[j] for j in range(3)) for i in range(3)]
        # L_v = {x : Bv . x = 0 mod p}: generators p e_i and kernel vectors
        gens = [[p if i == j else 0 for i in range(3)] for j in range(3)]
        for j in range(3):
            e = [1 if i == j else 0 for i in range(3)]
            if Bv[j] % p == 0:
                gens.append(e)
            else:
                for k2 in range(3):
                    if k2 != j:
                        c = (-Bv[k2] * pow(Bv[j], -1, p)) % p
                        g = [0, 0, 0]; g[j] = c; g[k2] = 1
                        gens.append(g)
        gens = [[p * x for x in g] for g in gens] + [v]          # scale by p, add v (i.e. v/p after /p)
        M = pari.matrix(3, len(gens), [gens[j][i] for i in range(3) for j in range(len(gens))])
        H = pari.mathnf(M)
        basis = [[pari.frac(0) + int(H[i, j]) for i in range(3)] for j in range(3)]
        Bq = [[0] * 3 for _ in range(3)]
        for a in range(3):
            for b in range(3):
                val = bv(B, [int(H[i, a]) for i in range(3)], [int(H[i, b]) for i in range(3)])
                assert val % (p * p) == 0
                Bq[a][b] = val // (p * p)
        out.append(reduce_form(Bq))
    return out


def classify(q, classes, buckets):
    key = tuple(int(x) for x in pari.qfrep(mat(q), 40))
    for idx in buckets.get(key, []):
        if isometric(classes[idx], q):
            return idx
    classes.append(q)
    buckets.setdefault(key, []).append(len(classes) - 1)
    return len(classes) - 1


def study(d, l, primes, Bsq=300):
    forms = []
    for f in glob.glob("inert_*.json"):
        for R in json.load(open(f)):
            if (R["d"], R["l"]) == (d, l):
                for g in R["genera"]:
                    forms += [x["q"] for x in g["forms"]]
    classes, buckets = [], {}
    for q in forms:
        classify(q, classes, buckets)
    n_glue = len(classes)
    A = {}
    for p in primes:
        rows = {}
        i = 0
        while i < len(classes):
            rows[i] = [classify(q, classes, buckets) for q in neighbours(classes[i], p)]
            i += 1
        n = len(classes)
        A[p] = [[sum(1 for j in rows[i] if j == k) for k in range(n)] for i in range(n)]
    n = len(classes)
    for p in primes:                                   # pad rows computed before the closure grew
        A[p] = [r + [0] * (n - len(r)) for r in A[p]]
    # square kernel on the closure
    S = []
    for q in classes:
        r = pari.qfrep(mat(q), Bsq * Bsq)
        S.append([int(r[m * m - 1]) for m in range(1, Bsq + 1)])
    Sm = pari.matrix(Bsq, n, [S[j][i] for i in range(Bsq) for j in range(n)])
    K = pari.matker(Sm)
    kdim = int(pari.matsize(K)[1])
    stable = True
    for p in primes:
        At = pari.mattranspose(pari.matrix(n, n, [x for r in A[p] for x in r]))
        img = At * K
        # img columns must lie in span(K): rank([K | img]) == rank(K)
        if kdim and int(pari.matrank(pari.matconcat([K, img]))) != kdim:
            stable = False
    # kernel restricted to gluing classes only (the user's relation space)
    Sg = pari.matrix(Bsq, n_glue, [S[j][i] for i in range(Bsq) for j in range(n_glue)])
    kg = int(pari.matsize(pari.matker(Sg))[1])
    # eigen-structure: characteristic polynomial factorisation of A_p0 and kernel dimension per factor
    p0 = primes[0]
    Ap = pari.matrix(n, n, [x for r in A[p0] for x in r])
    cp = pari.charpoly(pari.mattranspose(Ap))
    fac = pari.factor(cp)
    per_factor = []
    for k in range(int(pari.matsize(fac)[0])):
        P, e = fac[k, 0], int(fac[k, 1])
        Esp = pari.matker(pari.subst(P, "x", pari.mattranspose(Ap)) ** e)
        dimE = int(pari.matsize(Esp)[1])
        inter = int(pari.matsize(pari.matintersect(Esp, K))[1]) if kdim else 0
        per_factor.append({"poly": str(P), "mult": e, "dim": dimE, "square_null_dim": inter})
    return {"d": d, "l": l, "gluing_classes": n_glue, "closure_classes": n, "closed": n == n_glue,
            "primes": primes, "square_kernel_dim_closure": kdim, "square_kernel_dim_gluing": kg,
            "kernel_Hecke_stable": stable, "eigen_factors_T_p0": per_factor,
            "kernel_is_sum_of_eigenspaces": sum(f["square_null_dim"] for f in per_factor) == kdim}


if __name__ == "__main__":
    out = []
    for spec in sys.argv[1:]:
        d, l, ps = spec.split(":")
        r = study(int(d), int(l), [int(x) for x in ps.split(",")])
        out.append(r)
        print({k: v for k, v in r.items() if k != "eigen_factors_T_p0"}, flush=True)
        for f in r["eigen_factors_T_p0"]:
            if f["square_null_dim"]:
                print("   square-null eigen-factor:", f, flush=True)
    json.dump(out, open("hecke_" + "_".join(s.replace(":", "_").replace(",", "-") for s in sys.argv[1:])[:60] + ".json", "w"), indent=1)
