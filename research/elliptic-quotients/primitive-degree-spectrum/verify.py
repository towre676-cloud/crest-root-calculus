"""Replay the primitive parity calculation, isogeny sieve and degree counts."""
import json
from math import gcd, isqrt
from pathlib import Path
import sys

import sympy as sp

from degrees import RATIONAL_MODULUS, SIEVE_PRIMES, spectrum


def check(condition, why):
    if not condition:
        raise ArithmeticError(why)


def elliptic_traces(a, p):
    qa = [0]*p
    for y in range(p):
        qa[y*y % p] += 1
    c1 = 1 + sum(qa[x*(x-1)*(x-a) % p] for x in range(p))
    c2 = 2 + sum(qa[(x-1)*(x-a)*(x+1)*(x+1-a) % p] for x in range(p))
    check(len({1 % p, a % p, -1 % p, (a-1) % p}) == 4
          and len({0, 1, a % p}) == 3, "Singular reduction")
    return [p+1-c1, p+1-c2]


def verify():
    z = sp.Symbol("a")
    # Cross ratio of (1,-1,a,a-1) and the two Legendre j-invariants.
    cross = z*(z-1)/((z-2)*(z+1))
    def invariant(lam):
        return sp.factor(256*(1-lam+lam**2)**3/(lam**2*(1-lam)**2))
    j1, j2 = invariant(z), invariant(cross)
    check(j1.subs(z, 4) == sp.Rational(35152, 9), "First j-invariant")
    val = j2.subs(z, 4)
    check(val.q % 5 == 0 and val.p % 5 != 0, "Second j-invariant")
    # For a=4 modulo 15, B1 and B2 each have four branch values.
    # Track the three kernel tests with labels, independently of their values.
    infinity = "infinity"
    B1 = {0, 1, z, infinity}
    B2 = {-1, 0, z-1, infinity}
    killed = [set(), B1, B2, B1 ^ B2]
    for e in (0, 1, z):
        for odd_support in ({e, infinity}, {e-1, infinity}, {e, e-1}):
            check(odd_support not in killed, "A nonzero two-torsion class was killed")
    odd_primes = {3, 5, 7, 11, 13, 17, 19, 37, 43, 67, 163}
    left = [odd_primes.copy(), odd_primes.copy()]
    data = {}
    for p in SIEVE_PRIMES:
        tr = elliptic_traces(4, p)
        excluded = []
        for i, trace in enumerate(tr):
            delta = trace**2 - 4*p
            hits = sorted(ell for ell in left[i] if ell != p and pow(delta, (ell-1)//2, ell) == ell-1)
            left[i].difference_update(hits)
            excluded.append(hits)
        data[str(p)] = {"traces": tr, "newly_excluded_odd_primes": excluded}
    check(not any(left) and RATIONAL_MODULUS == 4849845, "Incomplete odd-isogeny sieve")
    for d in range(1, 501):
        sample = spectrum(4, d, True)
        half = d//2
        reps = [(m, n) for m in range(-isqrt(half), isqrt(half)+1)
                for n in range(-isqrt(half), isqrt(half)+1) if 2*(m*m+n*n) == d]
        check(sample["rational_map"] == bool(reps), "Rational degree witness")
        check(sample["primitive_classes_over_Q_and_Qbar"] ==
              int(d == 2) + sum(gcd(m, n) == 1 for m, n in reps)//2,
              "Primitive class count")
    check([spectrum(4, n)["all_geometric_classes"] for n in (2,3,4)] == [3,0,11],
          "Small-degree isogeny multiplicities")
    return {"j_at_a4": [str(j1.subs(z, 4)), str(j2.subs(z, 4))],
            "primitive_family_modulus": 15, "stronger_rational_modulus": RATIONAL_MODULUS,
            "witnesses": data, "first_examples": [spectrum(4, d, True) for d in (2, 4, 6, 8, 10, 26)]}


if __name__ == "__main__":
    report = verify()
    if sys.argv[1:] == ["--compute"]:
        print(json.dumps(report, indent=2))
    elif not sys.argv[1:] and report == json.loads(Path(__file__).with_name("receipt.json").read_text()):
        print("Verified primitive counts, arithmetic degree sieve and first 500 degrees.")
    else:
        raise SystemExit("Supply --compute to write the result or replay against the stored receipt.")
