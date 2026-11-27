"""Certify a local square exclusion by exact sublattice isometries."""

import argparse
from collections import defaultdict
from itertools import product
from math import lcm
import json
from pathlib import Path

from local_square import case_data, delay_census, pari
from neighbours import neighbour_basis, quadratic


def restriction(gram, basis):
    return basis.mattranspose() * gram * basis


def preimage(columns, modulus):
    identity = [tuple(modulus * int(i == j) for i in range(3)) for j in range(3)]
    return pari.mathnf(pari.Mat([pari.Col(v) for v in identity + list(columns)]))


def split_frame(gram, p):
    """Find coordinates with norm congruent to a unit times xy modulo p²."""
    if p <= 2 or not pari.isprime(p) or gram.matdet() % (p * p):
        raise ValueError("Expected an odd prime whose square divides the determinant.")
    kernel = (pari.Mod(1, p) * gram).matker()
    if len(kernel) != 1:
        raise ValueError("Expected a rank-two reduction.")
    left = right = None
    for vector in product(range(p), repeat=3):
        v = pari.Col(vector)
        pairings = gram * v
        if not any(pairings % p) or quadratic(gram, vector) % p:
            continue
        if left is None:
            left = v
        elif int((left.mattranspose() * pairings)[0]) % p:
            right = v
            break
    if right is None:
        raise ValueError("The rank-two reduction is nonsplit.")

    lifted = []
    for v in (left, right):
        pairings = list(map(int, gram * v))
        j = next(i for i, x in enumerate(pairings) if x % p)
        value = int((v.mattranspose() * gram * v)[0])
        entries = list(map(int, v))
        entries[j] += p * (-value // p * pow(2 * pairings[j], -1, p) % p)
        lifted.append(tuple(entries))
    u, v = lifted
    pair_u, pair_v = (tuple(map(int, gram * pari.Col(w))) for w in (u, v))
    r = delay_census.cross(pair_u, pair_v)
    frame = pari.Mat([pari.Col(w) for w in (u, v, r)])
    local = restriction(gram, frame)
    if (not frame.matdet() % p or any(local[i, i] % (p * p) for i in range(3))
            or local[0, 2] or local[1, 2] or not local[0, 1] % p):
        raise ArithmeticError("Failed to construct a split local frame.")
    return (u, v, r), int(2 * local[0, 1]) % p


def split_terms(gram, p):
    """Three signed lists for the radical, square lines, and p² divisibility."""
    (u, v, r), coefficient = split_frame(gram, p)
    yield "radical", 1, preimage([r], p)
    for t in range(1, p):
        if pow(coefficient * t % p, (p - 1) // 2, p) == 1:
            yield "unit_square_lines", 1, preimage([tuple(a + t * b for a, b in zip(u, v)), r], p)
    for a, b, sign in ((2, 0, 1), (1, 1, 1), (0, 2, 1), (2, 1, -1), (1, 2, -1)):
        yield "norm_divisible_by_p_squared", sign, preimage(
            [tuple(p**a * w for w in u), tuple(p**b * w for w in v), r], p * p
        )


def dyadic_terms(gram):
    """Represent the unit-square coset and norm-divisible-by-four sublattice."""
    residues = tuple(product(range(2), repeat=3))
    norms = {v: int(quadratic(gram, v)) for v in residues}
    if any(n % 4 not in (0, 1) for n in norms.values()):
        raise ValueError("Expected all norms to be zero or one modulo four.")
    even = [v for v, n in norms.items() if n % 4 == 0]
    yield "norm_divisible_by_p_squared", 1, preimage(even, 2)
    units = [v for v, n in norms.items() if n % 8 == 1]
    if not units:
        return
    representative = units[0]
    directions = [tuple((a - b) % 2 for a, b in zip(v, representative)) for v in units]
    low = preimage(directions, 2)
    high = preimage([*directions, representative], 2)
    if abs(low.matdet()) * len(units) != 8 or abs(low.matdet()) != 2 * abs(high.matdet()):
        raise ValueError("The unit-square residues do not form one affine coset.")
    yield "unit_square_coset", 1, high
    yield "unit_square_coset", -1, low


def next_dyadic_terms(gram):
    """Resolve the next valuation layer on an integral Gram form."""
    if gram.denominator() != 1:
        raise ValueError("The second-layer test requires an integral Gram matrix.")
    residues = tuple(product(range(2), repeat=3))
    divisible = [v for v in residues if quadratic(gram, v) % 4 == 0]
    preimage_of_zero_set = preimage(divisible, 2)
    label = "second_norm_divisible_by_four"
    if abs(preimage_of_zero_set.matdet()) * len(divisible) == 8:
        yield label, 1, preimage_of_zero_set / 2
    else:
        yield label, 1, pari.matid(3)
        for v in divisible[1:]:
            yield label, 1, preimage([v], 2) / 2
            yield label, -1, pari.matid(3)

    for residue in residues:
        if quadratic(gram, residue) % 4 != 1:
            continue
        lifts = (tuple(a + 2 * b for a, b in zip(residue, r)) for r in residues)
        units = [v for v in lifts if quadratic(gram, v) % 8 == 1]
        if not units:
            continue
        representative = units[0]
        directions = [tuple((a - b) % 4 for a, b in zip(v, representative)) for v in units]
        low = preimage(directions, 4)
        high = preimage([*directions, representative], 4)
        if abs(low.matdet()) * len(units) != 64 or abs(low.matdet()) != 2 * abs(high.matdet()):
            raise ArithmeticError("Failed to resolve the second unit-square coset.")
        yield "second_unit_square_cosets", 1, high
        yield "second_unit_square_cosets", -1, low


def neighbour_terms(gram, edges, p, layers=1):
    if layers != 1 and (p != 2 or layers != 2):
        raise ValueError("Only the two-adic criterion has a second layer.")
    terms = defaultdict(list)
    for edge, sign in zip(edges, (1, -1), strict=True):
        for point in edge:
            basis = neighbour_basis(gram, 3, point)
            q = restriction(gram, basis)
            projected = dyadic_terms(q) if p == 2 else split_terms(q, p)
            for name, multiplicity, sub in projected:
                terms[name].append((sign * multiplicity, basis * sub))
    if layers == 2:
        first_divisible = terms.pop("norm_divisible_by_p_squared")
        for sign, raw in first_divisible:
            basis = raw / 2
            for name, multiplicity, sub in next_dyadic_terms(restriction(gram, basis)):
                terms[name].append((sign * multiplicity, basis * sub))
    return terms


class UnmatchedSublattices(Exception):
    """The available isometries do not give a cancellation certificate."""


def match_terms(gram, terms, pencil=None):
    """Return the exact isometries pairing the positive and negative terms."""
    matches = {}
    for name, summands in terms.items():
        positive = [(i, b) for i, (sign, b) in enumerate(summands) if sign == 1]
        negative = [(i, b) for i, (sign, b) in enumerate(summands) if sign == -1]
        if len(positive) != len(negative):
            raise UnmatchedSublattices(f"Unbalanced signed terms in {name}.")
        matches[name] = []
        for i, a in positive:
            left = restriction(gram, a)
            left_pencil = restriction(pencil, a) if pencil is not None else 0
            for j, (index, b) in enumerate(negative):
                right = restriction(gram, b)
                right_pencil = restriction(pencil, b) if pencil is not None else 0
                # A large perturbation proposes a simultaneous isometry. Test both matrices.
                search_left = left + 531 * left_pencil
                search_right = right + 531 * right_pencil
                scale = lcm(int(search_left.denominator()), int(search_right.denominator()))
                transform = (scale * search_right).qfisom(scale * search_left)
                if not transform or restriction(left, transform) != right:
                    continue
                if pencil is not None and restriction(left_pencil, transform) != right_pencil:
                    continue
                if abs(transform.matdet()) != 1:
                    raise ArithmeticError("Found a non-unimodular matching matrix.")
                matches[name].append([i, index, delay_census.matrix_entries(transform)])
                del negative[j]
                break
            else:
                raise UnmatchedSublattices(f"No certified match in {name} for term {i}.")
    return matches


def family(gram, edges, terms, step, child_norm, neighbour_norms, expected_coefficients):
    """Check the bounded ingredients of the uniform deformation argument."""
    z_squared = delay_census.gram_matrix(((0, 0, 0), (0, 0, 0), (0, 0, 1)))
    matches = match_terms(gram, terms, z_squared)
    parents = [gram, gram + step * z_squared]
    children = [[delay_census.gram_matrix(delay_census.sublattice(
        delay_census.matrix_entries(g), delay_census.cross(*edge), 3
    )) for edge in edges] for g in parents]
    if any(child.qfrep(1)[0] for child in children[0]):
        raise ArithmeticError("The starting pair contains a product.")
    child_counts = [[2 * int(child.qfrep(child_norm)[child_norm - 1]) for child in pair] for pair in children]
    if child_counts != [[2, 0], [2, 0]]:
        raise ArithmeticError("The children lack the ordinary theta witness.")
    neighbour_counts = []
    for g, n in zip(parents, neighbour_norms, strict=True):
        forms = [restriction(g, neighbour_basis(g, 3, v)) for edge in edges for v in edge]
        values = [int(q.qfrep(n)[n-1]) for q in forms]
        neighbour_counts.append(2 * (sum(values[:2]) - sum(values[2:])))
    if neighbour_counts != expected_coefficients:
        raise ArithmeticError("The signed neighbour sums need nonzero ordinary witnesses.")
    return {"step": step, "child_witness_norm": child_norm, "child_counts": child_counts,
            "neighbour_witness_norms": neighbour_norms, "signed_neighbour_coefficients": neighbour_counts,
            "simultaneous_matches": matches}


def analyse(records):
    receipts = {2: [], 5: []}
    deformations = {}
    for index, entry in enumerate(records):
        tested = next((test for test in entry.get("local_obstructions", [])
                       if test["prime"] == 5 or test["orbit_dimension"] <= 2), None)
        if tested is None:
            continue
        prime = tested["prime"]
        layers = tested["orbit_dimension"]
        parent, _theta, _primes = case_data(entry)
        _, a, b, c, u, v, w = parent[0]
        gram = delay_census.gram_matrix(((a, u, v), (u, b, w), (v, w, c)))
        terms = neighbour_terms(gram, parent[1], prime, layers)
        try:
            matching = match_terms(gram, terms)
        except UnmatchedSublattices:
            if layers != 2:
                raise
            continue
        receipts[prime].append({"case": index, "parent_determinant": parent[0][0],
                                "valuation_layers": layers, "matches": matching})
        if parent[0] == (4000, 8, 13, 41, 2, 2, 3):
            deformations["five_adic"] = family(gram, parent[1], terms, 1800, 8, [5, 8], [-2, -4])
        if parent[0] == (1024, 5, 13, 17, 1, 1, -3):
            deformations["dyadic_one_layer"] = family(gram, parent[1], terms, 72, 5, [5, 5], [-4, -4])
        if parent[0] == (10240, 8, 32, 41, 0, 2, 4):
            deformations["dyadic_two_layers"] = family(gram, parent[1], terms, 288, 8, [8, 8], [-4, -4])
    if [len(receipts[p]) for p in (2, 5)] != [37, 18] or len(deformations) != 3:
        raise ArithmeticError("Incomplete census or missing deformation seed.")
    return {"dyadic_cases": receipts[2], "five_adic_cases": receipts[5], "deformations": deformations}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    records = json.loads(Path(__file__).with_name("local-square.json").read_text())
    result = json.dumps(analyse(records), indent=2)
    if args.verify and json.loads(result) != json.loads(args.verify.read_text()):
        raise ArithmeticError("The isometry replay differs from its receipt.")
    print(result)


if __name__ == "__main__":
    main()
