"""Independent replay of stored torsion-denominator certificates.

Checks, for each stored certificate: the Weil polynomial by direct point counts, geometric block distinctness,
Mumford validity, annihilation on a spanning set at each kill level, and the
nonzero images at witness levels. The geometric torsion has rank four away
from p and rank two in the ordinary p-primary case. Each full degree is
recomputed from the local levels and the recorded degree closure.  It also checks the explicit degree-4 map identically.
"""
import json
from pathlib import Path

import sympy as sp

from jac import Jacobian, odd_model, weil_poly, real_weil_roots, geometrically_distinct, vval, prime_factors
from certify import phi_apply

HERE = Path(__file__).resolve().parent


def decode(J, rec):
    F, R = J.F, J.R
    g = F.gen()
    elem = lambda cs: sum((int(c) * g ** i for i, c in enumerate(cs)), F(0))
    poly = lambda L: R([elem(c) for c in L]) if L else R(0)
    u, v = poly(rec["u"]), poly(rec["v"])
    assert 0 <= u.degree() <= 2 and u.leading_coefficient().is_one() and v.degree() < u.degree()
    assert ((J.f - v*v) % u).is_zero()
    return u, v


def span_size(J, gens, bound):
    seen = {J.key(J.zero)}
    frontier = [J.zero]
    while frontier:
        nxt = []
        for D in frontier:
            for G in gens:
                E = J.add(D, G)
                k = J.key(E)
                if k not in seen:
                    seen.add(k)
                    nxt.append(E)
                    if len(seen) > bound:
                        return len(seen)
        frontier = nxt
    return len(seen)


def check_certificate(cert):
    p = cert["p"]
    assert p > 2 and sp.isprime(p)
    g, _ = odd_model(cert["curve"], p)
    Jacobian(g, p, 1)
    P = weil_poly(g, p)
    assert P == cert["weil"], "Weil polynomial mismatch"
    a1, a2 = real_weil_roots(P, p)
    assert a1 != a2 and geometrically_distinct(a1, a2, p)
    assert sorted(B["trace"] for B in cert["blocks"]) == sorted((a1, a2))
    degrees = []
    for B in cert["blocks"]:
        a, b, N = B["trace"], B["complement_trace"], B["N"]
        assert {a, b} == {a1, a2} and N == abs(a - b)
        assert set(map(int, B["primes"])) == set(prime_factors(N))
        d = 1
        for ls, rec in B["primes"].items():
            ell = int(ls)
            assert rec["v_N"] == vval(N, ell)
            assert "undetermined_level" not in rec or rec.get("closure")
            if ell == p:
                assert a % p and b % p
            rank = 2 if ell == p else 4
            for K in rec["kills"]:
                J = Jacobian(g, p, K["r"])
                gens = [decode(J, G) for G in K["generators"]]
                k = K["level"]
                for G in gens:
                    assert J.is_zero(J.mul(ell ** k, G)) and J.is_zero(phi_apply(J, G, b))
                order = ell ** (rank * k)
                assert span_size(J, gens, order) == order, "kill generators do not span the geometric torsion"
            levels = [K["level"] for K in rec["kills"]]
            assert levels == list(range(1, len(levels) + 1))
            s = len(levels)
            assert rec["s"] == s and 0 <= s <= rec["v_N"]
            witness = rec.get("witness")
            if witness:
                assert witness["level"] == s + 1 <= rec["v_N"]
                J = Jacobian(g, p, witness["r"])
                W = decode(J, witness["point"])
                image = phi_apply(J, W, b)
                assert J.is_zero(J.mul(ell ** witness["level"], W)) and not J.is_zero(image)
                assert image == decode(J, witness["phi_point"])
            if rec.get("closure"):
                assert N == ell ** rec["v_N"] and s == rec["v_N"] - 1
            elif not witness:
                assert s == rec["v_N"]
            v_d = rec["v_N"] - s
            assert rec["v_d"] == v_d
            d *= ell ** v_d
        assert d // p ** vval(d, p) == B["prime_to_p_degree"]
        assert d == B.get("degree", B["prime_to_p_degree"])
        degrees.append(d)
    assert len(set(degrees)) == 1
    return degrees


def check_map():
    z = sp.symbols("z")
    m = json.loads((HERE / "map4.json").read_text())
    u = sp.sympify(m["u"], locals={"z": z})
    Yw = sp.sympify(m["Y_times_w"], locals={"z": z})
    w2 = 11 * (z**2 - 10*z + 13) * (2*z**2 + 2*z - 7) * (7*z**2 + 18*z + 3)
    assert sp.cancel(Yw**2 / w2 - (u - 1) * (u - 4) * (u - 9)) == 0
    num, den = sp.fraction(sp.cancel(u))
    return max(sp.degree(num, z), sp.degree(den, z))


if __name__ == "__main__":
    for filename in ("certificates.json", "residue_examples.json"):
        certs = json.loads((HERE / filename).read_text())
        for c in certs:
            if filename == "residue_examples.json":
                assert all(block["primes"][str(c["p"])]["witness"] for block in c["blocks"])
            print(c["name"], "p =", c["p"], "degrees", check_certificate(c))
    print("explicit map C' -> E1: u has degree", check_map(), "in z, so the map has degree 4")
    print("all checks passed")
