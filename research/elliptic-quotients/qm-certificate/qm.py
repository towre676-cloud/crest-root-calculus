"""Sample, interpolate and recognise psi_alpha(P) = alpha(P - W0) as an explicit morphism C -> Sym^2 C.

Mumford form of psi(P): u = T^2 - s1 T + s0, v = y (g1 T + g0) with s1, s0, g1, g0 in L(x) (psi o iota = -psi).
"""
from __future__ import annotations

import json
import random
import sys

import cypari2
import mpmath as mp
import numpy as np

from invert import make, mp_aj_routed, solve_divisor

pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9)


def tangent(C, R):
    Pi = C.Pi
    PR = Pi * mp.matrix(R)
    A = mp.matrix([[PR[i, j] for j in range(2)] for i in range(2)])
    B = mp.matrix([[Pi[i, j] for j in range(2)] for i in range(2)])
    M = A * B ** -1
    # check M Pi = Pi R on all four columns
    err = max(abs((M * Pi - PR)[i, j]) for i in range(2) for j in range(4))
    return M, err


def sample(C, F, M, n, rng, box=1.6):
    pts = []
    tries = 0
    while len(pts) < n and tries < 20 * n:
        tries += 1
        x = mp.mpc(rng.uniform(-box, box), rng.uniform(-box, box))
        try:
            A, y = mp_aj_routed(C, F.router, x)
        except (ValueError, ZeroDivisionError):
            continue
        z = [M[0, 0] * A[0] + M[0, 1] * A[1], M[1, 0] * A[0] + M[1, 1] * A[1]]
        r = solve_divisor(C, F, z, rng)
        if r is None:
            continue
        (x1, x2), (y1, y2) = r
        if abs(x1 - x2) < mp.mpf(10) ** -5:
            continue
        s1, s0 = x1 + x2, x1 * x2
        r1 = (y1 - y2) / (x1 - x2)
        r0 = y1 - r1 * x1
        pts.append({"x": x, "y": y, "s1": s1, "s0": s0, "g1": r1 / y, "g0": r0 / y})
        print("sample", len(pts), flush=True)
    return pts


def fit(pts, key, deg):
    """N(x) - s D(x) = 0, deg N, deg D <= deg, D monic of degree <= deg: nullspace by SVD."""
    rows = []
    for p in pts:
        x, s = p["x"], p[key]
        rows.append([x ** i for i in range(deg + 1)] + [-s * x ** i for i in range(deg + 1)])
    A = mp.matrix(rows)
    U, S, V = mp.svd_c(A)
    sv = sorted([S[i] for i in range(len(S))])
    v = [mp.conj(V[V.rows - 1, j]) for j in range(V.cols)]
    N, D = v[:deg + 1], v[deg + 1:]
    k = max(i for i in range(deg + 1) if abs(D[i]) > mp.mpf(10) ** -20)
    s = D[k]
    return [c / s for c in N], [c / s for c in D], sv


if __name__ == "__main__":
    f = [-1, 5, -8, 4, -1, 1, 0]
    which = sys.argv[1]
    R = {"alpha": [[1, 0, 1, 0], [0, 0, 0, 2], [1, 0, -1, 0], [0, 1, 0, 0]],
         "beta": [[1, 2, -2, -2], [-2, -1, 2, 0], [0, 0, -1, -2], [0, 0, 2, 1]],
         "gamma": [[0, 0, 1, 1], [-1, -1, -1, 0], [0, -1, -1, -1], [1, 1, 0, 0]]}[which]
    n = int(sys.argv[2]) if len(sys.argv) > 2 else 40
    C, F = make(f, dps=int(sys.argv[3]) if len(sys.argv) > 3 else 40)
    M, err = tangent(C, R)
    print("tangent", [[mp.nstr(M[i, j], 25) for j in range(2)] for i in range(2)], "check", mp.nstr(err, 3), flush=True)
    rng = random.Random(7)
    pts = sample(C, F, M, n, rng)
    print("samples", len(pts), flush=True)
    json.dump([{k: [str(mp.re(v)), str(mp.im(v))] for k, v in p.items()} for p in pts], open(f"samples_{which}.json", "w"))
    for key in ("s1", "s0", "g1", "g0"):
        for deg in range(1, 16):
            if 2 * deg + 3 > len(pts):
                break
            N, D, sv = fit(pts, key, deg)
            if sv[0] < mp.mpf(10) ** (-(mp.mp.dps - 12)) and sv[1] > mp.mpf(10) ** -8:
                print(key, "degree", deg, "sv", mp.nstr(sv[0], 3), mp.nstr(sv[1], 3))
                print("  N", [mp.nstr(c, 20) for c in N])
                print("  D", [mp.nstr(c, 20) for c in D])
                break
