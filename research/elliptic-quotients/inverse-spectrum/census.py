"""Exact square-spectrum census using the closed inert-level form."""
import argparse
from collections import defaultdict
from itertools import product
import json
from math import isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "inert-gluing"))
from verify import ResidueField, degree_counts, humbert_orbit, square_shells  # noqa: E402


CLASS_NUMBER_ONE = (-3, -4, -7, -8, -11, -19, -43, -67, -163)


def residue_field(discriminant, prime):
    trace = discriminant % 2
    norm = (trace*trace-discriminant)//4
    if prime < 3 or prime % 2 == 0 or any(prime % d == 0 for d in range(3, isqrt(prime)+1, 2)):
        raise ValueError("An odd prime is required.")
    field = ResidueField(prime, trace, norm)
    if field.legendre(discriminant) != -1:
        raise ValueError("The prime must be inert.")
    return field


def scalar_orbits(field):
    t, n = field.trace_omega, field.norm_omega
    units = [field.reduce((x, y)) for x, y in product(range(-1, 2), repeat=2)
             if x*x+t*x*y+n*y*y == 1]
    unseen = {a for a in product(range(field.p), repeat=2) if (1+field.norm(a)) % field.p}
    orbits = {}
    for a in sorted(unseen):
        if a in unseen:
            orbit = humbert_orbit(field, units, a)
            unseen.difference_update(orbit)
            orbits[a] = orbit
    return orbits


def separation(discriminant, prime, cutoff):
    """Partition by exact multiplicities until every pair separates or cutoff is hit."""
    field = residue_field(discriminant, prime)
    orbits = scalar_orbits(field)
    genera = defaultdict(list)
    for a in orbits:
        genera[field.legendre(1+field.norm(a))].append(a)
    unresolved = [v for v in genera.values() if len(v) > 1]
    latest = 0
    witness = None
    for d, shell in enumerate(square_shells(field, cutoff), 1):
        if not unresolved:
            break
        finer = []
        for group in unresolved:
            by_count = defaultdict(list)
            for a in group:
                by_count[degree_counts(field, a, [shell])[0]].append(a)
            if len(by_count) > 1:
                first, second = sorted(by_count)[:2]
                witness = {"centers": [by_count[first][0], by_count[second][0]],
                           "degree": d, "counts": [first, second]}
                latest = d
            finer.extend(v for v in by_count.values() if len(v) > 1)
        unresolved = finer
    return {"D": discriminant, "ell": prime, "forms": len(orbits),
            "genus_sizes": {k: len(v) for k, v in sorted(genera.items())},
            "last_first_difference": latest, "latest_witness": witness,
            "unresolved_through_cutoff": unresolved, "cutoff": cutoff}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-limit", type=int, default=109)
    parser.add_argument("--minimum-degree-cutoff", type=int, default=65)
    parser.add_argument("--maximum-degree", type=int, default=512)
    args = parser.parse_args()
    if (args.prime_limit < 3 or args.minimum_degree_cutoff < 1
            or args.maximum_degree < args.minimum_degree_cutoff):
        parser.error("Use prime limit >= 3 and positive, ordered degree cutoffs.")
    results = []
    for discriminant in CLASS_NUMBER_ONE:
        for prime in range(3, args.prime_limit+1, 2):
            try:
                residue_field(discriminant, prime)
            except ValueError:
                continue
            cutoff = min(args.maximum_degree, max(2*prime, args.minimum_degree_cutoff))
            result = separation(discriminant, prime, cutoff)
            while result["unresolved_through_cutoff"] and cutoff < args.maximum_degree:
                cutoff = min(2*cutoff, args.maximum_degree)
                result = separation(discriminant, prime, cutoff)
            results.append(result)
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
