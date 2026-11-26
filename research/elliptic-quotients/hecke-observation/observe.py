"""Certify separation using only quotient degrees that are powers of one auxiliary prime."""
import argparse
from collections import defaultdict
import json
from math import gcd, isqrt
from pathlib import Path
import sys

import numpy as np

from neighbours import Closure

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "inverse-spectrum"))
from census import residue_field, scalar_orbits  # noqa: E402
from verify import degree_counts, form_gram, square_shells  # noqa: E402


MODULUS = 65521


def refine(groups, exponent, h, centers):
    unresolved = []
    witness = None
    for group in groups:
        by_value = defaultdict(list)
        for j in group:
            by_value[int(h[j])].append(j)
        if len(by_value) > 1:
            first, second = sorted(by_value)[:2]
            witness = {"exponent": exponent, "centers": [centers[by_value[first][0]], centers[by_value[second][0]]],
                       "primitive_counts_modulus": [first, second]}
        unresolved.extend(g for g in by_value.values() if len(g) > 1)
    return unresolved, witness


def run(D, ell, p):
    field = residue_field(D, ell)
    if p < 3 or gcd(p, 2*D*ell) != 1 or any(p % j == 0 for j in range(2, isqrt(p) + 1)):
        raise ValueError("Use an odd good auxiliary prime.")
    centers = list(scalar_orbits(field))
    closures = []
    genus_indices = []
    for sign in (1, -1):
        indices = [i for i, a in enumerate(centers) if field.legendre(1 + field.norm(a)) == sign]
        seeds = [form_gram(field, centers[i]) for i in indices]
        closures.append(Closure(seeds, (p,)))
        genus_indices.append(indices)
    operators = [np.array(closure.rows[p], dtype=np.int64) for closure in closures]
    current = [np.array([2 * int(q.qfrep(1)[0]) for q in closure.forms], dtype=np.int64) for closure in closures]
    previous = [np.zeros_like(h) for h in current]
    if any(any(h[:closure.nseeds]) for h, closure in zip(current, closures)):
        raise ArithmeticError("A gluing lattice represents norm one.")
    chi = pow(D % p, (p - 1)//2, p)
    chi = -1 if chi == p - 1 else chi
    within = [[group] for group in genus_indices]
    joint = [list(range(len(centers)))]
    last_within = [None, None]
    last_joint = None
    shells = square_shells(field, p*p)
    for exponent in range(sum(len(c.forms) for c in closures)):
        values = np.zeros(len(centers), dtype=np.int64)
        for gi, closure in enumerate(closures):
            count = current[gi] - previous[gi]
            values[genus_indices[gi]] = count[:closure.nseeds] % MODULUS
        if exponent <= 2:
            degree = p**exponent
            expected = [degree_counts(field, a, [shells[degree-1]])[0] for a in centers]
            if values.tolist() != [v % MODULUS for v in expected]:
                raise ArithmeticError("Neighbour recurrence disagrees with direct primitive enumeration.")
        for gi in range(2):
            within[gi], witness = refine(within[gi], exponent, values, centers)
            last_within[gi] = witness or last_within[gi]
        joint, witness = refine(joint, exponent, values, centers)
        last_joint = witness or last_joint
        if not joint and exponent >= 2:
            break
        for gi, op in enumerate(operators):
            correction = chi * current[gi] if exponent == 0 else p * previous[gi]
            following = (current[gi][op].sum(axis=1) - correction) % MODULUS
            previous[gi], current[gi] = current[gi], following
    return {
        "D": D, "ell": ell, "p": p, "modulus": MODULUS, "tested_through_exponent": exponent,
        "forms": len(centers), "genera": [
            {"sign": sign, "forms": len(indices), "closure": len(closure.forms),
             "norm_one_classes": sum(int(q.qfrep(1)[0]) != 0 for q in closure.forms),
             "last_witness": last_within[gi],
             "unresolved": [[centers[i] for i in group] for group in within[gi]]}
            for gi, (sign, indices, closure) in enumerate(zip((1, -1), genus_indices, closures))
        ],
        "last_joint_witness": last_joint, "joint_unresolved": [[centers[i] for i in group] for group in joint],
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("levels", nargs="*", help="discriminant:gluing prime:auxiliary prime")
    parser.add_argument("--verify", type=Path, help="Recompute and check the levels in this receipt.")
    args = parser.parse_args()
    if bool(args.levels) == bool(args.verify):
        parser.error("Provide levels or --verify RECEIPT.")
    recorded = json.loads(args.verify.read_text()) if args.verify else []
    levels = [f"{r['D']}:{r['ell']}:{r['p']}" for r in recorded] if args.verify else args.levels
    results = []
    for spec in levels:
        result = run(*map(int, spec.split(":")))
        results.append(result)
        print(spec, result["forms"], result["tested_through_exponent"], len(result["joint_unresolved"]),
              file=sys.stderr, flush=True)
    rendered = json.dumps(results, indent=2)
    if args.verify and json.loads(rendered) != recorded:
        raise ArithmeticError("The recomputed receipt differs.")
    print(rendered)


if __name__ == "__main__":
    main()
