"""Certify absence of geometric elliptic factors using one good reduction."""

import argparse
import json

import flint
import sympy as sp

from jac import _all_elements, power_sums


EXPONENT = 240
x = sp.Symbol("x")


def good_model(f, p):
    """Require an odd good prime for the stated quintic or sextic model."""
    if p == 2 or not sp.isprime(p) or len(f) not in (6, 7) or f[-1] % p == 0:
        return False
    poly = sp.Poly.from_list(f[::-1], x, modulus=p)
    return poly.is_sqf


def count_points(f, p, r):
    """Count affine points and the rational branches at infinity."""
    field = flint.fq_default_ctx(p, r)
    polynomial = flint.fq_default_poly_ctx(field)([field(c) for c in f])
    count = 1 if len(f) == 6 else 2*int(field(f[-1]).is_square())
    for t in _all_elements(field, p, r):
        rhs = polynomial(t)
        count += 1 if rhs.is_zero() else 2*int(rhs.is_square())
    return count


def weil_polynomial(f, p):
    t1 = p+1-count_points(f, p, 1)
    t2 = p*p+1-count_points(f, p, 2)
    return [1, -t1, (t1*t1-t2)//2, -p*t1, p*p]


def power_polynomial(weil, p, exponent):
    powers = power_sums(weil, p, 4*exponent)
    es = [1]
    for k in range(1, 5):
        es.append(sum((-1)**(i-1)*es[k-i]*powers[exponent*i]
                      for i in range(1, k+1))//k)
    return sp.Poly.from_list([(-1)**i*coefficient for i, coefficient in enumerate(es)], x)


def simplicity_certificate(f, primes=range(3, 200)):
    for p in primes:
        if not good_model(f, p):
            continue
        weil = weil_polynomial(f, p)
        if not sp.Poly.from_list(weil, x).is_irreducible:
            continue
        if power_polynomial(weil, p, EXPONENT).is_irreducible:
            return {"p": p, "weil": weil, "exponent": EXPONENT, "simple_over_Fpbar": True}
    return None


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coefficients", type=json.loads, help="Integer coefficients in ascending order")
    parser.add_argument("--prime", type=int, action="append", help="Candidate prime, repeatable")
    args = parser.parse_args()
    print(json.dumps(simplicity_certificate(args.coefficients, args.prime or range(3, 200)), indent=2))
