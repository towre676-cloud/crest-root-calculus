"""Family sweep: Richelot neighbours of y^2 = (x^2-a^2)(x^2-b^2)(x^2-c^2).

Each neighbour has a degree-four map from its Richelot correspondence.
If a good odd prime has geometrically distinct elliptic blocks, the common
optimal degree divides four and the torsion certificate determines it there.  An independent oracle
decides whether a non-hyperelliptic involution exists (d = 2 iff it does, for genus two).
"""
import json
import sys
import time
from multiprocessing import Pool

import mpmath as mp
import sympy as sp

from certify import certify_curve
from richelot import pairings, richelot, x

mp.mp.dps = 60


def involutions(coeffs):
    """Count Mobius involutions permuting the roots of f (high-precision numerical oracle)."""
    c = [mp.mpf(int(t)) for t in coeffs]
    roots = list(mp.polyroots(c[::-1], maxsteps=400, extraprec=400))
    if len(roots) == 5:
        roots.append(mp.inf)
    shift = mp.mpf("0.3141592653589793") + mp.mpf("0.2718281828") * 1j
    R = [mp.mpf(0) if r == mp.inf else 1 / (r - shift) for r in roots]
    count = 0
    for P in pairings(list(range(6))):
        (a, b), (cc, d), (e, f) = [(R[i], R[j]) for i, j in P]
        r1 = [a + b, 1, -a * b]
        r2 = [cc + d, 1, -cc * d]
        A = r1[1] * r2[2] - r1[2] * r2[1]
        B = r1[2] * r2[0] - r1[0] * r2[2]
        Cc = r1[0] * r2[1] - r1[1] * r2[0]
        scale = max(abs(A), abs(B), abs(Cc))
        A, B, Cc = A / scale, B / scale, Cc / scale
        if abs(A * A + B * Cc) < mp.mpf(10) ** -30:
            continue
        if abs(A * (e + f) + B - Cc * e * f) < mp.mpf(10) ** -30:
            count += 1
    return count


def good_primes(start=23, stop=400):
    p = start
    while p < stop:
        if sp.isprime(p):
            yield p
        p += 1


def certify_degree(coeffs, need=1, budget=60):
    """Find good multiplicity-one primes until `need` coprime determinations agree or conflict."""
    t0 = time.time()
    rows, degrees = [], []
    for p in good_primes():
        if time.time() - t0 > budget or len(degrees) >= need:
            break
        try:
            res, _ = certify_curve(coeffs, p, rmax=6, sylow_limit=2 ** 16)
        except Exception as e:  # bad reduction, no rational root, etc.
            continue
        if res["status"] != "certified":
            continue
        ds = [B["prime_to_p_degree"] for B in res["blocks"]]
        rows.append({"p": p, "traces": [B["trace"] for B in res["blocks"]], "N": res["blocks"][0]["N"], "degrees": ds})
        if None not in ds:
            degrees.append((p, ds))
    return rows, degrees


def job(args):
    triple, P = args
    G = [sp.expand((x - a) * (x - b)) for a, b in P]
    D, H = richelot(G)
    if D == 0:
        return None
    poly = sp.Poly(sp.expand(H[0] * H[1] * H[2]) * sp.sign(D) * abs(D), x)  # y^2 = D*H1H2H3, same curve over Q
    coeffs = [int(c) for c in poly.all_coeffs()[::-1]]
    g = int(sp.gcd_list(coeffs))
    coeffs = [c // g for c in coeffs]  # scaling by a constant only twists; blocks' degrees are geometric
    t = time.time()
    rows, degrees = certify_degree(coeffs)
    inv = involutions(coeffs)
    # Corollary 4: the optimal degree divides 4 and is a power of two, so one odd
    # multiplicity-one prime determines it.
    final = None
    if degrees:
        vals = {tuple(d) for _, d in degrees}
        final = list(list(vals)[0]) if len(vals) == 1 else "conflict"
    return {"source": list(triple), "pairing": [[int(a) for a in t] for t in P], "Delta": int(D), "coeffs": coeffs,
            "local": rows, "certified_degrees": final, "involutions": inv,
            "seconds": round(time.time() - t, 1)}


if __name__ == "__main__":
    triples = [(1, 2, 3), (1, 2, 4), (1, 3, 4), (1, 2, 5), (2, 3, 5)]
    tasks = [(T, P) for T in triples for P in pairings([T[0], -T[0], T[1], -T[1], T[2], -T[2]])]
    out = []
    with Pool(2) as pool, open("sweep.jsonl", "w") as fh:
        for r in pool.imap_unordered(job, tasks, chunksize=1):
            if r:
                out.append(r)
                fh.write(json.dumps(r) + "\n")
                fh.flush()
    out.sort(key=lambda r: (r["source"], r["pairing"]))
    json.dump(out, open("sweep.json", "w"), indent=1)
    agree = sum(1 for r in out if r["certified_degrees"] not in (None, "conflict")
                and ((r["certified_degrees"][0] == 2) == (r["involutions"] > 0))
                and r["certified_degrees"][0] == r["certified_degrees"][1])
    print(len(out), "curves;", agree, "agree with the involution oracle")
    for r in out:
        print(r["source"], r["pairing"], r["certified_degrees"], "inv", r["involutions"], [(l["p"], l["degrees"]) for l in r["local"]], r["seconds"])
