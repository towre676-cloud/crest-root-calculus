"""Independent replay of stored torsion-denominator certificates.

Checks, for each stored certificate: the Weil polynomial by direct point counts, block
traces and geometric distinctness, that every kill generator is l^k-torsion with
phi(G) = 0, that the kill generators span a group of order l^(4k) (so they generate all
of J[l^k]), and that each witness W satisfies l^k W = 0 and phi(W) != 0.  The degree is
then recomputed from the levels.  It also checks the explicit degree-4 map identically.
"""
import json
import sys

import sympy as sp

from jac import Jacobian, odd_model, weil_poly, real_weil_roots, geometrically_distinct, vval, prime_factors
from certify import phi_apply


def decode(J, rec):
    F, R = J.F, J.R
    g = F.gen()
    elem = lambda cs: sum((int(c) * g ** i for i, c in enumerate(cs)), F(0))
    poly = lambda L: R([elem(c) for c in L]) if L else R(0)
    return (poly(rec["u"]), poly(rec["v"]))


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
    g, _ = odd_model(cert["curve"], p)
    P = weil_poly(g, p)
    assert P == cert["weil"], "Weil polynomial mismatch"
    a1, a2 = real_weil_roots(P, p)
    assert a1 != a2 and geometrically_distinct(a1, a2, p)
    degrees = []
    for B in cert["blocks"]:
        a, b, N = B["trace"], B["complement_trace"], B["N"]
        assert {a, b} == {a1, a2} and N == abs(a - b)
        d = 1
        for ls, rec in B["primes"].items():
            ell = int(ls)
            for K in rec["kills"]:
                J = Jacobian(g, p, K["r"])
                gens = [decode(J, G) for G in K["generators"]]
                k = K["level"]
                for G in gens:
                    assert J.is_zero(J.mul(ell ** k, G)) and J.is_zero(phi_apply(J, G, b))
                assert span_size(J, gens, ell ** (4 * k)) == ell ** (4 * k), "kill generators do not span J[l^k]"
            levels = [K["level"] for K in rec["kills"]]
            assert levels == list(range(1, len(levels) + 1))
            s = len(levels)
            if rec.get("witness"):
                W = rec["witness"]
                assert W["level"] == s + 1
                J = Jacobian(g, p, W["r"])
                Wd = decode(J, W["point"])
                assert J.is_zero(J.mul(ell ** W["level"], Wd)) and not J.is_zero(phi_apply(J, Wd, b))
                d *= ell ** (rec["v_N"] - s)
            elif rec.get("closure"):
                assert s == rec["v_N"] - 1 and len(B["primes"]) == 1
                d *= ell
            else:
                assert s == rec["v_N"]
        assert d == B["prime_to_p_degree"]
        degrees.append(d)
    return degrees


def check_map():
    z = sp.symbols("z")
    m = json.load(open("map4.json"))
    u = sp.sympify(m["u"], locals={"z": z})
    Yw = sp.sympify(m["Y_times_w"], locals={"z": z})
    w2 = 11 * (z**2 - 10*z + 13) * (2*z**2 + 2*z - 7) * (7*z**2 + 18*z + 3)
    assert sp.cancel(Yw**2 / w2 - (u - 1) * (u - 4) * (u - 9)) == 0
    num, den = sp.fraction(sp.cancel(u))
    return max(sp.degree(num, z), sp.degree(den, z))


if __name__ == "__main__":
    certs = json.load(open("certificates.json"))
    for c in certs:
        print(c["name"], "p =", c["p"], "degrees", check_certificate(c))
    print("explicit map C' -> E1: u has degree", check_map(), "in z, so the map has degree 4")
    print("all checks passed")
