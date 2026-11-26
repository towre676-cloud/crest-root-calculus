"""Torsion-denominator certificates for elliptic blocks of genus-two Jacobians.

For a good prime p at which the real Weil polynomial splits as (s-a)(s-b) with the two
elliptic blocks geometrically non-isogenous, put phi = F + V - b and N = |a - b|.  The
symmetric projector onto the a-block is e = phi/(a-b).  The optimal degree d of the
special elliptic quotient in the a-block satisfies, for every prime l != p dividing N,

    v_l(d) = v_l(N) - s_l,   s_l = max{ s <= v_l(N) : phi kills J[l^s] }.

In residue characteristic the same test uses geometric p-power torsion when
p divides N: both blocks are ordinary and phi is Rosati self-adjoint. For a
genus-two Jacobian over F_p in this case N=p, so d=p already follows from
d | N and d >= 2. An optional witness checks the torsion explanation.

A certificate for s_l consists of
  * a level-s_l kill record: a field F_{p^r} over which J[l^s] is rational, its Smith
    generators, and phi(generator) = 0 for each; and
  * unless s_l = v_l(N), a witness: a point W over some F_{p^r'} with l^(s+1) W = 0 and
    phi(W) != 0.
"""
from __future__ import annotations

import argparse
import json
import random
import time

from sympy import isprime

from jac import (Jacobian, SylowStructure, odd_model, weil_poly, jac_order, real_weil_roots,
                 geometrically_distinct, vval, prime_factors)


def phi_apply(J, D, b):
    """(F + V - b)(D) on J(F_q), using V = [p] F^{r-1} on F_q-points."""
    FD = J.frob(D, 1)
    VD = J.mul(J.p, J.frob(D, J.r - 1)) if J.r > 1 else J.mul(J.p, D)
    return J.add(J.add(FD, VD), J.mul(-b, D))


def mumford_record(D):
    """Coefficients low to high; each F_q coefficient is its list over the FLINT power basis."""
    enc = lambda poly: [[int(t) for t in c.to_list()] for c in poly.coeffs()]
    return {"u": enc(D[0]), "v": enc(D[1])}


def primary_order(J, prime, point):
    """Logarithm of the order of a point known to have prime-power order."""
    level = 0
    while not J.is_zero(point):
        point = J.mul(prime, point)
        level += 1
    return level


def cyclic_log(J, prime, generator, level, point):
    """Extract base-prime digits, or return None outside the cyclic subgroup."""
    first_digit = J.mul(prime**(level-1), generator)
    digits = {J.key(J.mul(digit, first_digit)): digit for digit in range(prime)}
    result = 0
    for i in range(level):
        digit = digits.get(J.key(J.mul(prime**(level-1-i), point)))
        if digit is None:
            return None
        contribution = digit*prime**i
        point = J.add(point, J.mul(-contribution, generator))
        result += contribution
    return result if J.is_zero(point) else None


class Certifier:
    def __init__(self, fcoeffs, p, rmax=6, sylow_limit=2 ** 17, seed=1):
        if p <= 2 or not isprime(p):
            raise ValueError("Use an odd prime for the hyperelliptic model.")
        self.p = p
        self.g, self.root = odd_model(fcoeffs, p)
        Jacobian(self.g, p, 1)  # Validate smoothness before using point counts.
        self.P = weil_poly(self.g, p)
        self.rmax, self.limit, self.seed = rmax, sylow_limit, seed
        self._cache = {}

    def sylow(self, ell, r):
        key = (ell, r)
        if key not in self._cache:
            n = jac_order(self.P, self.p, r)
            if ell ** vval(n, ell) > self.limit:
                self._cache[key] = None
            else:
                J = Jacobian(self.g, self.p, r)
                self._cache[key] = (J, SylowStructure(J, ell, n, seed=self.seed, max_size=self.limit))
        return self._cache[key]

    def quick_witness(self, ell, k, b, trials=24, rmax=8):
        """Randomized witness search without building the Sylow subgroup.

        For random R in J(F_q), P = (n / l^v) R has l-power order l^j; if j >= k then
        Q = l^(j-k) P lies in J[l^k] and is a witness when phi(Q) != 0.  Only positive
        outcomes are used; failure to find a witness proves nothing.
        """
        rng = random.Random(7 + k)
        for r in range(1, rmax + 1):
            n = jac_order(self.P, self.p, r)
            vv = vval(n, ell)
            if vv < k:
                continue
            J = Jacobian(self.g, self.p, r)
            cof = n // ell ** vv
            for _ in range(trials):
                P = J.mul(cof, J.random_divisor(rng))
                j = primary_order(J, ell, P)
                if j < k:
                    continue
                Q = J.mul(ell ** (j - k), P)
                I = phi_apply(J, Q, b)
                if not J.is_zero(I):
                    return {"level": k, "r": r, "point": mumford_record(Q), "phi_point": mumford_record(I),
                            "method": "random projection"}
        return None

    def residue_witness(self, a, b):
        """Find a p-torsion witness in the ordinary N=p case."""
        p = self.p
        # The two unit roots have the same residue a, so the Frobenius
        # matrix on the rank-two geometric p-torsion has this eigenvalue twice.
        r = p*next(i for i in range(1, p) if pow(a, i, p) == 1)
        n = jac_order(self.P, p, r)
        cofactor = n//p**vval(n, p)
        J = Jacobian(self.g, p, r)
        rng = random.Random(self.seed)
        generator, level = J.zero, 0

        def witness(point):
            image = phi_apply(J, point, b)
            if not J.is_zero(image):
                modulus = [int(c) for c in J.F.modulus().coeffs()]
                return {"level": 1, "r": r, "field_modulus": modulus,
                        "point": mumford_record(point), "phi_point": mumford_record(image),
                        "method": "ordinary torsion projection"}
            return None

        for _ in range(24):
            point = J.mul(cofactor, J.random_divisor(rng))
            height = primary_order(J, p, point)
            if not height:
                continue
            found = witness(J.mul(p**(height-1), point))
            if found:
                return found
            if height > level:
                generator, level = point, height
                continue

            # Take the last nonzero step modulo the pivot subgroup. Removing
            # a multiple of the pivot then gives a second p-torsion direction.
            while cyclic_log(J, p, generator, level, point) is None:
                next_point = J.mul(p, point)
                log = cyclic_log(J, p, generator, level, next_point)
                if log is not None:
                    if log % p:
                        raise ArithmeticError("The pivot should have maximal sampled order.")
                    return witness(J.add(point, J.mul(-(log//p), generator)))
                point = next_point
        return None

    def blocks(self):
        h = real_weil_roots(self.P, self.p)
        if h is None or h[0] == h[1]:
            return None
        if not geometrically_distinct(h[0], h[1], self.p):
            return None
        return h

    def certify_block(self, a, b, *, residue_witness=False):
        N = abs(a - b)
        out = {"trace": a, "complement_trace": b, "N": N, "primes": {}}
        d = 1
        for ell in prime_factors(N):
            if ell == self.p:
                if N != ell or a % ell == 0 or b % ell == 0:
                    raise ArithmeticError("Unexpected traces for a distinct pair with p dividing N.")
                witness = self.residue_witness(a, b) if residue_witness else None
                out["primes"][str(ell)] = {
                    "v_N": 1, "s": 0, "v_d": 1, "kills": [], "witness": witness,
                    "closure": "N=p and degree is at least 2",
                }
                continue
            v = vval(N, ell)
            rec = {"v_N": v, "kills": [], "witness": None}
            s = 0
            for k in range(1, v + 1):
                verdict = None
                qw = self.quick_witness(ell, k, b)
                if qw is not None:
                    verdict = ("witness", qw)
                for r in ([] if verdict else range(1, self.rmax + 1)):
                    n = jac_order(self.P, self.p, r)
                    if n % ell:
                        continue
                    syl = self.sylow(ell, r)
                    if syl is None:
                        continue
                    J, S = syl
                    gens = S.torsion_generators(k)
                    for G in gens:
                        assert J.is_zero(J.mul(ell ** k, G))
                    imgs = [phi_apply(J, G, b) for G in gens]
                    bad = [G for G, I in zip(gens, imgs) if not J.is_zero(I)]
                    if bad:
                        W = bad[0]
                        verdict = ("witness", {"level": k, "r": r, "point": mumford_record(W),
                                               "phi_point": mumford_record(phi_apply(J, W, b))})
                        break
                    if S.torsion_rank(k) == 4:
                        verdict = ("kill", {"level": k, "r": r, "sylow_invariants": S.d,
                                            "generators": [mumford_record(G) for G in gens]})
                        break
                if verdict is None:
                    rec["undetermined_level"] = k
                    break
                if verdict[0] == "witness":
                    rec["witness"] = verdict[1]
                    break
                rec["kills"].append(verdict[1])
                s = k
            rec["s"] = s
            if "undetermined_level" not in rec:
                rec["v_d"] = v - s
                d *= ell ** (v - s)
            out["primes"][str(ell)] = rec
        # Closure: a smooth curve of genus >= 2 has no degree-one map, so d >= 2.  If the only
        # undetermined level is the top level v of a single prime and everything else gives 1,
        # that level cannot be a kill.
        und = [(l, r) for l, r in out["primes"].items() if "undetermined_level" in r]
        if und:
            if len(und) == 1 and d == 1:
                l, r = und[0]
                if r["undetermined_level"] == r["v_N"] and r["s"] == r["v_N"] - 1:
                    r["v_d"] = 1
                    del r["undetermined_level"]
                    r["closure"] = "level v_N excluded because the optimal degree is at least 2"
                    d *= int(l)
                    und = []
            if und:
                out["prime_to_p_degree"] = None
                out["degree"] = None
                out["prime_to_p_degree_bound"] = {"multiple_of": d}
                return out
        out["prime_to_p_degree"] = d
        out["degree"] = d * (self.p if N % self.p == 0 else 1)
        return out

    def frobenius_mod(self, ell, k):
        """Matrix of F on J[ell^k] in Smith-derived generators (requires rationality)."""
        for r in range(1, self.rmax + 1):
            n = jac_order(self.P, self.p, r)
            if vval(n, ell) < 4 * k:
                continue
            syl = self.sylow(ell, r)
            if syl is None:
                continue
            J, S = syl
            if S.torsion_rank(k) != 4:
                continue
            m = ell ** k
            scale = [dd // m for dd in S.d]
            gens = S.torsion_generators(k)
            cols = []
            for G in gens:
                y = S.dlog(J.frob(G))
                col = []
                for yi, sc in zip(y, scale):
                    assert yi % sc == 0
                    col.append((yi // sc) % m)
                cols.append(col)
            M = [[cols[j][i] for j in range(4)] for i in range(4)]
            return r, M
        return None


def certify_curve(fcoeffs, p, *, residue_witness=False, **kw):
    c = Certifier(fcoeffs, p, **kw)
    h = c.blocks()
    res = {"p": p, "weil": c.P, "odd_model": c.g, "moved_root": c.root}
    if h is None:
        res["status"] = "not multiplicity one"
        return res, c
    a1, a2 = h
    res["status"] = "certified"
    res["blocks"] = [c.certify_block(a, b, residue_witness=residue_witness)
                     for a, b in ((a1, a2), (a2, a1))]
    return res, c


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coefficients", type=json.loads, help="JSON coefficient list, ascending")
    parser.add_argument("prime", type=int)
    parser.add_argument("--residue-witness", action="store_true", help="Try to find an explicit p-torsion witness")
    args = parser.parse_args()
    t = time.monotonic()
    res, _ = certify_curve(args.coefficients, args.prime, residue_witness=args.residue_witness)
    res["seconds"] = round(time.monotonic() - t, 2)
    print(json.dumps(res, indent=2))
