"""Certify the opposite-secant twins and uniform bounds on forced relations."""

import argparse
from collections import Counter
from fractions import Fraction
from itertools import combinations
import json
from math import gcd
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "forced-square-kernel"))
from star import (CLASS_NUMBER_ONE, act, ambient_genus, cross, fiber, gram_matrix,
                  matrix_entries, normal_partition, orbits, pari, residue_field, sublattice)
from rank import bound_on_genus
from neighbours import Closure
from observe import scalar_orbits


def secant_orbits(gram, automorphisms, field):
    vertices, _ = normal_partition(gram, field)
    actions = [{x: act(g, x, field.p) for x in vertices} for g in automorphisms]
    labels, sizes = orbits(combinations(sorted(vertices), 2), [
        lambda e, table=table: tuple(sorted(table[x] for x in e)) for table in actions
    ])
    return vertices, labels, sizes


def twin_pairs(gram, automorphisms, field):
    vertices, labels, _ = secant_orbits(gram, automorphisms, field)
    internal = normal_partition(gram, field)[1]
    reps = {}
    opposite = {}
    for edge, index in labels.items():
        complement = tuple(sorted(set(vertices).difference(edge)))
        normal = cross(*edge)
        inside = tuple(v for v in internal if sum(x*y for x, y in zip(v, normal)) % 3 == 0)
        other_normal = cross(*complement)
        other_inside = tuple(v for v in internal if sum(x*y for x, y in zip(v, other_normal)) % 3 == 0)
        if len(inside) != 1 or inside != other_inside:
            raise ArithmeticError("Opposite secants have different internal incidence.")
        target = labels[complement]
        if opposite.setdefault(index, target) != target:
            raise ArithmeticError("Opposite secants did not descend to automorphism orbits.")
        reps.setdefault(index, gram_matrix(sublattice(gram, normal, 3)))

    for i, j in opposite.items():
        if i < j:
            yield reps[i], reps[j]


def first_difference(left, right, cutoff):
    left_reps, right_reps = left.qfrep(cutoff**2), right.qfrep(cutoff**2)
    first = next(m for m in range(1, cutoff+1) if left_reps[m*m-1] != right_reps[m*m-1])
    if first % 3 or any(left_reps[m*m-1] != right_reps[m*m-1]
                        for m in range(1, cutoff+1) if m % 3):
        raise ArithmeticError("The exact square counts disagree with the local proof.")
    return first


def level_three(D, ambient, archive):
    field = residue_field(D, 3)
    genera = {"split": [], "nonsplit": []}
    source_split = source_nonsplit = forced = 0
    witnesses = Counter()
    twins = []
    for index, (gram, automorphisms) in enumerate(ambient):
        block = fiber(gram, automorphisms, field, genera)
        pairs = list(twin_pairs(gram, automorphisms, field))
        if len(pairs) != block["forced_kernel"]:
            raise ArithmeticError("Opposite pairs did not exhaust the forced subspace.")
        if index == 0:
            source_split = block["split_classes"]
            source_nonsplit = block["nonsplit_classes"]
            if pairs:
                raise ArithmeticError("The source ambient class contains an opposite-secant collision.")
        forced += len(pairs)
        for left, right in pairs:
            if left.qfrep(1)[0] or right.qfrep(1)[0] or left.qfisom(right):
                raise ArithmeticError("Expected distinct non-product refined Humbert forms.")
            degree = first_difference(left, right, 64)
            counts = left.qfrep(128**2)
            good_square = next(m for m in range(1, 129) if gcd(m, 6*D) == 1 and counts[m*m-1])
            witnesses[degree] += 1
            twins.append({"parent": index, "left": matrix_entries(left),
                          "right": matrix_entries(right), "first_degree": degree,
                          "good_square_witness": good_square})
    # Finite rank upper bounds and the geometric lower bounds agree, so no
    # additional same-genus relations are possible at any good degree.
    ranks = {}
    signatures = {}
    for genus in genera:
        lower = forced if genus == "split" else 0
        bound = bound_on_genus(genera[genus], D, 3, lower, 128)
        if bound["actual_kernel_upper"] != lower:
            raise ArithmeticError("The finite rank failed to certify the all-good-degree kernel.")
        ranks[genus] = bound["classes"] - lower
        sigs = set()
        for gram, _ in genera[genus]:
            counts = gram_matrix(gram).qfrep(128**2)
            sigs.add(tuple(int(counts[m*m-1]) for m in range(1, 129) if gcd(m, 6*D) == 1))
        signatures[genus] = sigs
    if signatures["split"] & signatures["nonsplit"]:
        raise ArithmeticError("The good square prefix failed to separate the two genera.")
    if archive is not None:
        unmatched = [r for r in archive if r["d"] == D and r["l"] == 3]
        for twin in twins:
            lhs, rhs = [gram_matrix(twin[key]) for key in ("left", "right")]
            for i, row in enumerate(unmatched):
                x, y = gram_matrix(row["q1"]), gram_matrix(row["q2"])
                if (lhs.qfisom(x) and rhs.qfisom(y)) or (lhs.qfisom(y) and rhs.qfisom(x)):
                    if row["first_square_difference"] != twin["first_degree"]:
                        raise ArithmeticError("The first differing degree changed from the archive.")
                    unmatched.pop(i)
                    break
            else:
                raise ArithmeticError("An opposite-secant pair is missing in the archived sweep.")
        if unmatched:
            raise ArithmeticError("Some archived pairs have no opposite-secant witness.")
    return {"D": D, "parent_classes": len(ambient), "split_classes": len(genera["split"]),
            "nonsplit_classes": len(genera["nonsplit"]), "observed_ranks": ranks,
            "source_classes": [source_split, source_nonsplit],
            "gluing_classes": len(scalar_orbits(field)), "twin_pairs": twins,
            "first_degree_counts": dict(sorted(witnesses.items()))}


def forced_bounds(ambients):
    results = []
    for D, ambient in ambients.items():
        group_bound = max(int(pari.qfauto(gram_matrix(gram))[0]) // 2 for gram, _ in ambient)
        # For larger primes every pair energy is at most 4 |Gamma/{+-1}|/(ell+1)
        # and hence below 1/2. Exact diagonals cover the finite remainder.
        prime_limit = 8 * group_bound - 1
        exact_checks = []
        for ell in range(5, prime_limit+1, 2):
            try:
                field = residue_field(D, ell)
            except ValueError:
                continue
            diagonals = [Fraction(fiber(gram, autos, field)["max_diagonal"]) for gram, autos in ambient]
            pair_bound = 2 * max(diagonals)
            if pair_bound >= 1:
                raise ArithmeticError("The diagonal bound permits a forced pair collision.")
            exact_checks.append({"ell": ell, "pair_energy_upper": str(pair_bound)})
        max_bound = max([Fraction(1, 2)] + [
            Fraction(row["pair_energy_upper"]) for row in exact_checks])
        results.append({"D": D, "max_projective_group_order": group_bound,
                        "finite_prime_limit": prime_limit, "exact_checks": exact_checks,
                        "all_ell_gt_3_pair_energy_upper": str(max_bound)})
    return results


def apply_operator(rows, prime, character, values, degree):
    zero = [0] * len(values)
    previous, current = zero, values
    for k in range(degree):
        neighbour_counts = [sum(current[j] for j in row) for row in rows]
        if k == 0:
            following = [n - character*x for n, x in zip(neighbour_counts, current)]
        else:
            following = [n - prime*x for n, x in zip(neighbour_counts, previous)]
        previous, current = current, following
    return current, previous


def walk_replay(D, ell, ambient):
    field = residue_field(D, ell)
    genera = {"split": [], "nonsplit": []}
    for gram, automorphisms in ambient:
        fiber(gram, automorphisms, field, genera)
    primes = tuple(p for p in [3, 5, 7, 11, 13, 17, 19] if gcd(p, 2*D*ell) == 1)[:2]
    p, q = primes
    p_character, q_character = (int(pari.kronecker(D, prime)) for prime in primes)
    exponents = [(0, 0), (1, 0), (2, 0), (0, 1), (0, 2), (1, 1), (2, 1), (1, 2)]
    degrees = [p**a * q**b for a, b in exponents]
    checked = 0
    for records in genera.values():
        closure = Closure([gram for gram, _ in records], primes)
        if len(closure.forms) != len(records):
            raise ArithmeticError("A good-prime neighbour left the full genus.")
        reps = [form.qfrep(max(degrees)**2) for form in closure.forms]
        seed = [2 * int(rep[0]) for rep in reps]
        if any(x not in (0, 2) for x in seed):
            raise ArithmeticError("Norm one should occur on a single signed pair.")
        for (a, b), m in zip(exponents, degrees):
            ordinary_p, old_p = apply_operator(closure.rows[p], p, p_character, seed, a)
            ordinary, _ = apply_operator(closure.rows[q], q, q_character, ordinary_p, b)
            primitive_p = [x-y for x, y in zip(ordinary_p, old_p)]
            primitive_q, old_q = apply_operator(closure.rows[q], q, q_character, primitive_p, b)
            predicted = [x-y for x, y in zip(primitive_q, old_q)]
            divisors = [(1, 1), (p, -1), (q, -1), (p*q, 1)]
            for index, rep in enumerate(reps):
                expected = 2*sum(mu * int(rep[(m//d)**2-1]) for d, mu in divisors if m % d == 0)
                if ordinary[index] != 2 * int(rep[m*m-1]) or predicted[index] != expected:
                    raise ArithmeticError("The Hecke walk differs from ordinary or primitive enumeration.")
                checked += 2
    return {"D": D, "ell": ell, "primes": primes, "degrees": sorted(degrees), "entries": checked}


def explicit_family():
    replays = []
    for degree in (9, 33, 57, 81):
        k = degree**2 - 4
        f = (25*k-1)//4
        parent = ((25, 0, 1), (0, 4, 0), (1, 0, k))
        if f % 300 != 181 or int(pari.qfauto(gram_matrix(parent))[0]) != 4:
            raise ArithmeticError("The parent congruence or automorphism group is incorrect.")
        left, right = [gram_matrix(sublattice(parent, row, 3))
                       for row in ((1, 0, 0), (1, 0, -1))]
        generated = list(twin_pairs(parent, [((1, 0, 0), (0, -1, 0), (0, 0, 1))],
                                    residue_field(-4*f, 3)))
        if len(generated) != 1:
            raise ArithmeticError("The explicit parent should give one opposite-secant pair.")
        a, b = generated[0]
        if not ((left.qfisom(a) and right.qfisom(b)) or (left.qfisom(b) and right.qfisom(a))):
            raise ArithmeticError("The explicit planes are not the opposite-secant pair.")
        if first_difference(left, right, degree) != degree or left.qfisom(right):
            raise ArithmeticError("Unexpected first distinguishing degree.")
        difference = 2*int(left.qfrep(degree**2)[degree**2-1] - right.qfrep(degree**2)[degree**2-1])
        if difference != 4:
            raise ArithmeticError("The first differing square should contribute four signed vectors.")
        replays.append({"degree": degree, "f": f, "left": matrix_entries(left),
                        "right": matrix_entries(right), "signed_difference": difference})
    d, a, b = 481, 964140, 43961
    u, v = a*a + d*b*b, 2*a*b
    x, y = 45*u + 2*d*v, 45*v + 2*u
    if u*u-d*v*v != 1 or (u % 120, v % 120) != (1, 0) or x*x-d*y*y != 101:
        raise ArithmeticError("The fixed-field Pell step failed its norm or congruence check.")
    return {"direct_replay": replays, "pell_norm_one_unit": [u, v],
            "next_pell_degree": x//5, "next_pell_conductor": y//2}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--archive", type=Path, help="Optional collision_counts.json from commit c754b46.")
    args = parser.parse_args()
    archive = json.loads(args.archive.read_text()) if args.archive else None
    pari.allocatemem(64_000_000, 768_000_000, silent=True)
    ambients = {D: ambient_genus(D) for D in CLASS_NUMBER_ONE}
    threes = [level_three(D, ambient, archive) for D, ambient in ambients.items()
              if int(pari.kronecker(D, 3)) == -1]
    result = {"level_three": threes, "forced_bounds": forced_bounds(ambients),
              "walk_replay": [walk_replay(D, ell, ambients[D])
                              for D, ell in [(-4, 3), (-19, 3), (-8, 13)]],
              "explicit_family": explicit_family()}
    output = json.dumps(result, indent=2)
    if args.verify and json.loads(output) != json.loads(args.verify.read_text()):
        raise ArithmeticError("The exact replay differs from the receipt.")
    print(output)


if __name__ == "__main__":
    main()
