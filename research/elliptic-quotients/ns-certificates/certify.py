"""Orchestrator: exact NS lattice, rank certificate, saturation and complete spectrum from providers."""
from __future__ import annotations

import itertools
import json
import random
from fractions import Fraction
from math import gcd, isqrt

import cypari2

from setup import Reduction, good_primes
from providers import MapProvider, SumProvider, IsogenyPair
from ns import lattice_basis, gram_of, spectrum_from_lattice, l_torsion, coords, point_order, RATIONAL_CM_J
from sympy.ntheory.modular import crt

pari = cypari2.Pari()


def scalar_on(E, fn, points, bound):
    """Integer t in [0, bound] with fn(P) = tP for all P (exact, by CRT on point orders)."""
    mod, val = 1, 0
    for P in points:
        if P is None:
            continue
        o = point_order(E, P, 10 ** 6)
        if o <= 2:
            continue
        g = fn(P)
        t = next((t for t in range(o) if E.eq(E.mul(t, P), g)), None)
        if t is None:
            raise AssertionError("not a scalar")
        r = crt([mod, o], [val, t])
        if r is None:
            raise AssertionError("inconsistent scalar")
        val, mod = int(r[0]), int(r[1])
        if mod > bound:
            break
    cands = [t for t in range(bound + 1) if t % mod == val % mod]
    if len(cands) != 1:
        raise RuntimeError(f"scalar not determined (modulus {mod}, bound {bound})")
    return cands[0]


def provider_degree(pv, bound):
    return scalar_on(pv.E, lambda P: pv.push(pv.pull(P)), pv.aux(80), bound)


def gram_from_providers(pvs, degs):
    n = len(pvs)
    G = [[0] * (n + 1) for _ in range(n + 1)]
    G[0][0] = 2
    for a in range(n):
        G[0][a + 1] = G[a + 1][0] = degs[a]
    for a in range(n):
        for b in range(a + 1, n):
            pa, pb = pvs[a], pvs[b]
            bound = degs[a] * degs[b]
            m = scalar_on(pb.E, lambda P: pb.push(pa.pull(pa.push(pb.pull(P)))), pb.aux(80), bound)
            G[a + 1][b + 1] = G[b + 1][a + 1] = degs[a] * degs[b] - m
    return G


def lpart(n, ell):
    m = 1
    while n % ell == 0:
        n //= ell
        m *= ell
    return m


def torsion_module(R, ell, s, limit=2 ** 18):
    """Smith-basis generators of J[l^s] (requires J[l^s] rational) and a coordinate function mod l^s."""
    from jac import SylowStructure
    if ell ** (4 * s) > limit:
        return None
    v = 0
    n = R.order
    while n % ell == 0:
        n //= ell
        v += 1
    if v < 4 * s or ell ** v > limit:
        return None
    S = SylowStructure(R.J, ell, R.order, max_size=limit)
    if S.torsion_rank(s) != 4:
        return None
    m = ell ** s
    gens = S.torsion_generators(s)
    scale = [d // m for d in S.d]

    def coord(D):
        y = S.dlog(D)
        out = []
        for yi, sc in zip(y, scale):
            assert yi % sc == 0, "point not in J[l^s]"
            out.append((yi // sc) % m)
        return out
    return gens, coord


def find_torsion_module(build, ell, s, need_i, f, recs, rmax=6):
    for p in good_primes(f, recs, need_i, start=11, stop=1500):
        if p == ell:
            continue
        for r in range(1, rmax + 1):
            try:
                R, pvs = build(p, r)
            except (AssertionError, ValueError, NotImplementedError):
                break
            T = torsion_module(R, ell, s)
            if T is not None:
                return R, pvs, T
    raise RuntimeError(f"no field with rational J[{ell}^{s}] within limits")


def saturate(build, n_gen, G, log, need_i, f, recs, rmax=6, replay=None):
    """Saturate the image lattice of (theta, D_1, ..., D_n) in NS, exactly.

    D = sum c_t B_t lies in l NS iff M alpha_D kills J[l m], where M is a common denominator of the
    coefficients of D in the generators and m its l-part (the prime-to-l part acts invertibly).
    """
    B = [[Fraction(x) for x in v] for v in lattice_basis(G)]
    n = len(G)
    degs = G[0][1:]
    while True:
        GB = gram_of(G, B)
        det = int(pari.matdet(pari.matrix(len(B), len(B), [x for r in GB for x in r])))
        cand = [int(q) for q in pari.factor(abs(det))[0] if det % (int(q) ** 2) == 0] if abs(det) > 1 else []
        changed = False
        for ell in cand:
            M = 1
            for v in B:
                for x in v:
                    M = M * x.denominator // gcd(M, x.denominator)
            m = lpart(M, ell)
            s = 1
            while ell ** (s - 1) < m:
                s += 1
            if replay is None:
                R, pvs, T = find_torsion_module(build, ell, s, need_i, f, recs, rmax)
            else:
                recorded = replay[len(log)]
                p, r = recorded["field"]
                if p == ell or (recorded["l"], recorded["level"]) != (ell, s):
                    raise ArithmeticError("Invalid recorded torsion field.")
                R, pvs = build(p, r)
                T = torsion_module(R, ell, s)
                if T is None:
                    raise ArithmeticError("The recorded field does not supply full torsion.")
            gens, coord = T
            J = R.J
            mod = ell ** s

            def mat(fn):
                cols = [coord(fn(g)) for g in gens]
                return [[cols[c][r] for c in range(4)] for r in range(4)]
            mats = [mat(lambda g: g)]
            for a, pv in enumerate(pvs):
                d = degs[a]
                mats.append(mat(lambda g, pv=pv, d=d: J.add(J.mul(d, g), J.neg(pv.pull(pv.push(g))))))
            hit = None
            for c in itertools.product(range(ell), repeat=len(B)):
                if not any(c):
                    continue
                coef = [sum(c[t] * B[t][k] for t in range(len(B))) * M for k in range(n)]
                assert all(x.denominator == 1 for x in coef)
                w = [int(x) % mod for x in coef]
                A = [[sum(w[k] * mats[k][r][q] for k in range(n)) % mod for q in range(4)] for r in range(4)]
                if all(A[r][q] == 0 for r in range(4) for q in range(4)):
                    hit = c
                    break
            log.append({"l": ell, "level": s, "field": [R.K.p, R.K.r], "divisible_combination": hit})
            if hit is not None:
                pos = next(t for t in range(len(B)) if hit[t] % ell)
                unit = pow(hit[pos], -1, ell)
                normalized = [(unit*c) % ell for c in hit]
                new = [sum(Fraction(normalized[t]) * B[t][k] for t in range(len(B))) / ell for k in range(n)]
                B[pos] = new
                changed = True
                break
        if not changed:
            return B, gram_of(G, B), det


def theta_coords(B):
    n = len(B[0])
    Bm = pari.matrix(n, len(B), [pari(str(B[c][r])) for r in range(n) for c in range(len(B))])
    e0 = pari.mattranspose(pari.matrix(1, n, [1] + [0] * (n - 1)))
    sol = pari.matsolve(pari.mattranspose(Bm) * Bm, pari.mattranspose(Bm) * e0)
    tv = [Fraction(str(sol[r, 0])) for r in range(len(B))]
    assert all(t.denominator == 1 for t in tv)
    return [int(t) for t in tv]
