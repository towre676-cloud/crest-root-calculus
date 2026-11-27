"""A height-bounded mod-3 Frobenius sieve in the marked 7- and 10-torsion families."""
from __future__ import annotations

import argparse
from collections import defaultdict
from math import gcd

import sympy as sp

from reconstruct import pari


def family(order, height):
    curves = {}
    for denominator in range(1, height+1):
        for numerator in range(-height, height+1):
            if gcd(numerator, denominator) != 1:
                continue
            t = pari(numerator)/denominator
            if order == 7:
                c = t*t-t
                b = t*c
            elif order == 9:
                c = t*t*(t-1)
                b = c*(t*t-t+1)
            elif order == 12:
                if t == 1:
                    continue
                c = (3*t*t-3*t+1)*(t-2*t*t)/(t-1)**3
                b = c*(2*t-2*t*t-1)/(t-1)
            elif order == 10:
                d = t-(t-1)**2
                c = (2*t**3-3*t*t+t)/d
                b = c*t*t/d
            else:
                raise ValueError(f"Unsupported marked torsion order: {order}")
            model = pari.ellinit([1-c, -b, -b, 0, 0])
            if len(model) == 0:
                continue
            assert int(model.ellorder([0, 0])) == order
            minimal = model.ellminimalmodel()[0]
            key = tuple(map(int, minimal[:5]))
            curves.setdefault(key, (minimal, (numerator, denominator)))
    return curves


HOWE_7 = (1, 0, 0, -5774401, 5346023177)
HOWE_10 = (1, 0, 0, 115, 561)
OTHER_10 = (1, 0, 0, -45, 81)


def trace(model, prime):
    return int(model.ellap(prime)) % 3 if int(model.disc()) % prime else None


def signature(model, primes):
    value = mask = 0
    for i, p in enumerate(primes):
        a = trace(model, p)
        if a is not None:
            value |= a << (2*i)
            mask |= 3 << (2*i)
    return value, mask


def sieve(height, prime_bound):
    e7, e10 = (family(order, height) for order in (7, 10))
    print(f"{len(e7)} order-7 curves × {len(e10)} order-10 curves = {len(e7)*len(e10)} pairs",
          flush=True)
    primes = list(sp.primerange(5, prime_bound+1))
    early, late = primes[:12], primes[12:]
    signatures_10 = [(key, *signature(model, early)) for key, (model, _) in e10.items()]
    survivors = []
    for a, (model, _) in e7.items():
        value, mask = signature(model, early)
        survivors.extend((a, b) for b, v, m in signatures_10 if not (value ^ v) & mask & m)
    print(f"{len(survivors)} survive the first {len(early)} primes", flush=True)
    for p in late:
        left = {k: trace(e7[k][0], p) for k in {a for a, _ in survivors}}
        right = {k: trace(e10[k][0], p) for k in {b for _, b in survivors}}
        survivors = [
            (a, b) for a, b in survivors
            if left[a] is None or right[b] is None or left[a] == right[b]
        ]
        if not survivors:
            break
    print(f"{len(survivors)} survivors after good primes through {prime_bound}")
    for a, b in survivors:
        print(a, b, "at parameters", e7[a][1], e10[b][1])
    return set(survivors)


def check_orientation_obstruction():
    """The alternate 10-torsion target is 2-isogenous and flips the pairing."""
    target = pari.ellinit(list(HOWE_10))
    assert int(pari.polgalois(target.elldivpol(3)/3)[0]) == 24
    point = [pari("-17/4"), pari("17/8")]
    assert target.ellisoncurve(point) and int(target.ellorder(point)) == 2
    quotient, _ = target.ellisogeny(point)
    assert tuple(map(int, pari.ellinit(quotient).ellminimalmodel()[0][:5])) == OTHER_10
    print("The second survivor is excluded by pairing orientation under the rational 2-isogeny.")


def projective_fields(curves):
    """Group by the irreducible nonlinear orbit of projective 3-torsion."""
    buckets = defaultdict(list)
    for minimal, (model, parameter) in curves.items():
        factors = list(model.elldivpol(3).factor()[0])
        degrees = sorted(int(f.poldegree()) for f in factors)
        assert degrees in ([1, 3], [1, 1, 2])
        orbit = max(factors, key=lambda f: int(f.poldegree()))
        key = (tuple(degrees), int(orbit.nfdisc()))
        buckets[key].append((minimal, parameter, orbit))
    return buckets


def local_lift_obstruction(first, second, prime_bound):
    """At this prime the first marked point cannot lift and the second can."""
    e9, e12 = (pari.ellinit(list(model)) for model in (first, second))
    for p in sp.primerange(5, prime_bound+1):
        if int(e9.disc()) % p == 0 or int(e12.disc()) % p == 0:
            continue
        n9, n12 = int(e9.ellcard(p)), int(e12.ellcard(p))
        if n9 % 27 == 0 or n12 % 9:
            continue
        assert n9 % 9 == 0
        invariants = list(map(int, e12.ellgroup(p)))
        if len(invariants) == 1 or invariants[1] % 3:
            return p, n9, n12, invariants
    return None


def search_108(height, prime_bound):
    nine, twelve = (family(order, height) for order in (9, 12))
    print(f"{len(nine)} order-9 curves × {len(twelve)} order-12 curves = {len(nine)*len(twelve)} pairs",
          flush=True)
    a, b = map(projective_fields, (nine, twelve))
    matches = []
    for key in sorted(a.keys() & b.keys()):
        for e9, t9, f in a[key]:
            for e12, t12, g in b[key]:
                if f.nfisisom(g) == 0:
                    continue
                witness = local_lift_obstruction(e9, e12, prime_bound)
                matches.append((t9, t12, witness))
                print(t9, t12, f"ruled out by {witness}" if witness else "no local obstruction found")
    print(f"{len(matches)} projective-field matches")
    return matches


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--height", type=int, default=100)
    parser.add_argument("--prime-bound", type=int, default=251)
    parser.add_argument("--order-108", action="store_true", help="Test the marked order-9 and order-12 families.")
    args = parser.parse_args()
    if args.height < 1 or args.prime_bound < 5:
        parser.error("use a positive height and a prime bound of at least 5")
    if args.order_108:
        matches = search_108(args.height, args.prime_bound)
        if (args.height, args.prime_bound) == (100, 251):
            assert len(matches) == 4 and all(witness for _, _, witness in matches)
    else:
        pairs = sieve(args.height, args.prime_bound)
        if args.height in (40, 100) and args.prime_bound == 251:
            assert pairs == {(HOWE_7, HOWE_10), (HOWE_7, OTHER_10)}
            check_orientation_obstruction()
