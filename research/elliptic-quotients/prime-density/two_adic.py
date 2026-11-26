"""Exclude the two extra candidate characters at larger discriminant -8 levels."""
import argparse
import json

from flint import nmod_mat

from certify import Closure, allowed_characters, form_gram, legendre, primes_up_to, residue_field, scalar_orbits


MODULUS = 65521


def certify(ell, primes, weights):
    field = residue_field(-8, ell)
    if not primes or len(primes) != len(set(primes)) or any(p in (2, ell) or p not in primes_up_to(p) for p in primes):
        raise ValueError("Use distinct good auxiliary primes.")
    odd_character = next(d for d in allowed_characters(-8, ell) if d % 2)
    other_characters = [d for d in allowed_characters(-8, ell) if d != odd_character]
    closure = Closure([form_gram(field, a) for a in scalar_orbits(field)], primes)
    n = len(closure.forms)
    base = [0] * (n*n)
    twisted = {d: base.copy() for d in other_characters}
    for p, weight in zip(primes, weights):
        signs = {d: legendre(d, p) for d in other_characters}
        for i, row in enumerate(closure.rows[p]):
            for j in row:
                offset = i*n + j
                base[offset] += weight
                for d in other_characters:
                    twisted[d][offset] += signs[d]*weight
    polynomial = nmod_mat(n, n, base, MODULUS).charpoly()
    degrees = {d: polynomial.gcd(nmod_mat(n, n, entries, MODULUS).charpoly()).degree()
               for d, entries in twisted.items()}
    if any(degrees.values()):
        raise ArithmeticError("A proposed character was not excluded.")
    return {"D":-8,"ell":ell,"modulus":MODULUS,"closure":n,"primes":primes,
            "weights":weights,"excluded_characters":other_characters,"gcd_degrees":degrees}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("levels", nargs="*", help="ell:prime,prime:weight,weight")
    args = parser.parse_args()
    levels = args.levels or ["23:3,31:1,2","47:3,5,7,11,13,17,19,23:1,2,3,4,5,6,7,8"]
    results = []
    for level in levels:
        ell, primes, weights = level.split(":")
        primes = list(map(int,primes.split(",")))
        weights = list(map(int,weights.split(",")))
        if len(primes) != len(weights):
            parser.error("Provide one weight per prime.")
        results.append(certify(int(ell),primes,weights))
    print(json.dumps(results,indent=2))


if __name__ == "__main__":
    main()
