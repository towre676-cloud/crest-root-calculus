"""Exhaustive observed-fraction sweep from the product formula (no representation counts beyond norm 1).

Per level: one admissible nonscalar gluing per genus as seeds, p-neighbour closure under the given good
primes (each genus closure is a union of p-neighbour components, which for the theorem only needs to be
Hecke-stable), product classes R (q represents 1), O_g = H . 1_{R_g}.  For EVERY pair of non-product classes
in one genus, the observed fraction ||P_{O_g} g||^2 / ||g||^2 is computed exactly; the minimum is a rigorous
lower bound for the true coprime observed fraction (finitely many Hecke operators span a submodule)."""
import itertools, json, random, sys, time
sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
import cypari2
from gluing import ns_and_q, reduce_form
from inverse import mat, is_scalar
from hecke import neighbours, classify
from discform import isomorphic as disc_iso

pari = cypari2.Pari()
pari.allocatemem(12 * 10 ** 9)


def seeds(d, l, rng):
    reps = []
    tried = 0
    while len(reps) < 2 and tried < 4000:
        tried += 1
        a, b, c, e = (rng.randrange(l) for _ in range(4))
        if (a * e - b * c + 1) % l or is_scalar(d, l, [[a, b], [c, e]]):
            continue
        r = ns_and_q(d, l, [[a, b], [c, e]])
        if r is None:
            continue
        q = reduce_form(r["q"])
        if int(pari.qfminim(mat(q), 1, 1)[0]) > 0:
            continue
        if not any(disc_iso(x, q) for x in reps):
            reps.append(q)
    return reps


def run(d, l, primes, maxh=600):
    t = time.time()
    rng = random.Random(l * 1000 + abs(d))
    S = seeds(d, l, rng)
    classes, buckets = [], {}
    for q in S:
        classify(q, classes, buckets)
    rows = {p: {} for p in primes}
    changed = True
    while changed:
        changed = False
        for p in primes:
            for i in range(len(classes)):
                if i not in rows[p]:
                    rows[p][i] = [classify(q, classes, buckets) for q in neighbours(classes[i], p)]
                    changed = True
            if len(classes) > maxh:
                return {"d": d, "l": l, "skipped": f"closure > {maxh}"}
    n = len(classes)
    A = {p: pari.matrix(n, n, [sum(1 for j in rows[p][i] if j == k) for i in range(n) for k in range(n)]) for p in primes}
    aut = [int(pari.qfauto(mat(q))[0]) for q in classes]
    r1 = [int(pari.qfrep(mat(q), 1)[0]) for q in classes]
    gen = []
    for q in classes:
        gen.append(0 if disc_iso(S[0], q) else 1)
    res = {"d": d, "l": l, "primes": primes, "closure": n, "genera": []}
    worst_all = Fr(1)
    for g in (0, 1):
        idx = [x for x in range(n) if gen[x] == g]
        pos = {x: a for a, x in enumerate(idx)}
        Ag = [pari.matrix(len(idx), len(idx), [int(A[p][x, y]) for x in idx for y in idx]) for p in primes]
        Rg = [x for x in idx if r1[x]]
        v = pari.matrix(len(idx), 1, [1 if x in Rg else 0 for x in idx])
        M, rk = v, 1
        while True:
            B = pari.matimage(pari.matconcat([M] + [op * M for op in Ag]))
            r2 = int(pari.matsize(B)[1])
            if r2 == rk:
                break
            M, rk = B, r2
        O = M
        W = pari.matrix(len(idx), len(idx), [0] * (len(idx) ** 2))
        for a, x in enumerate(idx):
            W[a, a] = pari("1/%d" % aut[x])
        Gi = pari.matsolve(pari.mattranspose(O) * W * O, pari.matid(rk))
        OW = pari.mattranspose(O) * W
        nonprod = [x for x in idx if not r1[x]]
        worst, cnt = None, 0
        for i, j in itertools.combinations(nonprod, 2):
            gvec = pari.matrix(len(idx), 1, [aut[x] if x == i else (-aut[x] if x == j else 0) for x in idx])
            b = OW * gvec
            fr = Fr(str((pari.mattranspose(b) * Gi * b)[0, 0])) / (aut[i] + aut[j])
            cnt += 1
            if worst is None or fr < worst:
                worst = fr
        worst_all = min(worst_all, worst) if worst is not None else worst_all
        res["genera"].append({"size": len(idx), "products": len(Rg), "dim_O": rk, "nonproduct_pairs": cnt,
                              "min_observed": str(worst), "min_float": float(worst) if worst is not None else None})
    res["min_observed_level"] = str(worst_all)
    res["all_separated"] = worst_all > 0
    res["seconds"] = round(time.time() - t, 1)
    return res


if __name__ == "__main__":
    out = []
    fn = sys.argv[1]
    for spec in sys.argv[2:]:
        d, l, ps = spec.split(":")
        r = run(int(d), int(l), [int(x) for x in ps.split(",")])
        out.append(r)
        if "skipped" in r:
            print(r, flush=True)
        else:
            print(f"d={d} l={l} closure={r['closure']} " + " | ".join(
                f"h={g['size']} R={g['products']} dimO={g['dim_O']} pairs={g['nonproduct_pairs']} min={g['min_observed']}"
                for g in r["genera"]) + f"  {r['seconds']}s", flush=True)
        json.dump(out, open(fn, "w"), indent=1)
