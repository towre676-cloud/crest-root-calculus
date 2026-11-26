"""Certify the bielliptic gluing and the Humbert form of each D=-19 model.

The exact class-polynomial test in descend.py makes the CM ideal list
exhaustive. Arb discards impossible classes and torsion labelings. The
surviving labels determine a rational lattice, which PARI checks exactly.
"""

import argparse
import itertools
import json
from fractions import Fraction as Fr
from math import lcm
from pathlib import Path
import sys

from cypari2 import Pari
from flint import acb, acb_poly, arb, ctx

from descend import L_POLY, compute as check_descent, poly

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "inverse-spectrum"))
from gluing import E_of, K, gram, q_from_ns, reduce_form
from polarised import ideals, ns_general

PARI = Pari()
D = -19
PRECISION = 1200
FIELD = K(D)
TWINS = ("q1", "q2")
CERTIFICATE = HERE / "certified_link_-19.json"


def require(condition, reason):
    if not condition:
        raise ArithmeticError(reason)


def read(path):
    return json.loads(path.read_text())


def integer_matrix(matrix, n):
    return [[int(matrix[i, j]) for j in range(n)] for i in range(n)]


def checked_roots(coefficients):
    degree = len(coefficients) - 1
    roots = acb_poly(coefficients).roots(
        tol=arb(2)**-(PRECISION - 100), maxprec=4 * PRECISION
    )
    require(
        len(roots) == degree
        and all(root.is_finite() for root in roots)
        and all(not a.overlaps(b) for a, b in itertools.combinations(roots, 2)),
        "Root isolation did not give the full set of disjoint balls.",
    )
    return roots


def eval_L(expression, eta):
    total = acb(0)
    for value in poly(expression).all_coeffs():
        total = total * eta + acb(arb(int(value.p)) / int(value.q))
    return total


def ideal_candidates(class_polynomial):
    # Every admissible j is already proved to satisfy H_-171.
    candidates = []
    w = (acb(D) + acb(D).sqrt()) / 2

    def as_complex(element):
        p, q = element
        return acb(arb(p.numerator)/p.denominator) + acb(arb(q.numerator)/q.denominator)*w

    for basis in ideals(D, 3):
        alpha, beta = map(as_complex, basis)
        tau = beta / alpha
        require(tau.imag > 0, "The ideal basis must be positively oriented.")
        j = tau.modular_j()
        require(j.is_finite(), "The modular j enclosure is not finite.")
        candidates.append((basis, tau, j))
    require(
        len(candidates) == len(class_polynomial) - 1
        and all(not a[2].overlaps(b[2]) for a, b in itertools.combinations(candidates, 2)),
        "The class-polynomial roots are not exhausted by distinct ideal classes.",
    )
    return candidates


def classify(g2, g3, candidates):
    j = 1728*g2**3/(g2**3 - 27*g3**2)
    require(j.is_finite(), "The curve j enclosure is not finite.")
    hits = [(basis, tau) for basis, tau, tested_j in candidates if j.overlaps(tested_j)]
    require(len(hits) == 1, f"The CM ideal class is ambiguous ({len(hits)} candidates).")
    return hits[0]


def weierstrass(a, b, c):
    """Put the cubic in 4X^3-g2*X-g3 form by X=U+a/3."""
    return -4*(b - a*a/3), -4*(2*a**3/27 - a*b/3 + c)


def half_period_labels(roots, a, g2, g3, tau):
    hp = [acb(1)/2, tau/2, (1 + tau)/2]
    g2_tau, g3_tau = acb.elliptic_invariants(tau)
    # Both invariants are nonzero for the conductor-three classes.
    scale = (g3/g3_tau)/(g2/g2_tau)
    values = [scale*acb.elliptic_p(h, tau) for h in hp]
    require(all(v.is_finite() for v in values), "The scaled half periods are not finite.")
    hits = [
        perm for perm in itertools.permutations(range(3))
        if all((root + a/3).overlaps(values[j]) for root, j in zip(roots, perm))
    ]
    require(len(hits) == 1, f"The half-period labeling is ambiguous ({len(hits)} candidates).")
    return hits[0]


def label_elliptic_quotient(s1, s2, candidates):
    a, b, c = -s1, s2, acb(-1)
    roots = checked_roots([c, b, a, 1])
    g2, g3 = weierstrass(a, b, c)
    basis, tau = classify(g2, g3, candidates)
    labels = half_period_labels(roots, a, g2, g3, tau)
    return roots, basis, labels


def reciprocal_permutation(roots1, roots2, labels1, labels2):
    perm = [None] * 3
    for r, h in zip(roots1, labels1):
        require(not r.contains(0), "A cubic root enclosure contains zero.")
        hits = [label for s, label in zip(roots2, labels2) if (1/r).overlaps(s)]
        require(len(hits) == 1, "The reciprocal-root match is ambiguous.")
        perm[h] = hits[0]
    require(len(set(perm)) == 3, "The reciprocal-root map must be bijective.")
    return perm


def encode_basis(basis):
    return [[str(x) for x in row] for row in basis]


def half_period(basis, index):
    alpha, beta = basis
    values = [alpha, beta, tuple(a+b for a, b in zip(alpha, beta))]
    return tuple(x/2 for x in values[index])


def glued_lattice(a1, a2, perm):
    zero = (Fr(0), Fr(0))
    gens = [u + zero for u in a1] + [zero + v for v in a2]
    gens.extend(half_period(a1, h) + half_period(a2, j) for h, j in enumerate(perm))
    den = lcm(*(x.denominator for row in gens for x in row))
    matrix = PARI.matrix(4, len(gens), [int(row[i]*den) for i in range(4) for row in gens])
    hnf = PARI.mathnf(matrix)
    return [tuple(Fr(int(hnf[i, j]), den) for i in range(4)) for j in range(4)]


def humbert_from_gluing(a1, a2, perm):
    basis = glued_lattice(a1, a2, perm)

    def principal(a):
        alpha, beta = a
        return 1/abs(FIELD.mul(alpha, FIELD.conj(beta))[1])

    theta = (2*principal(a1), 2*principal(a2), Fr(0), Fr(0))
    et = [[E_of(FIELD, theta, x, y) for y in basis] for x in basis]
    require(all(e.denominator == 1 for row in et for e in row), "The polarization is not integral.")
    em = PARI.matrix(4, 4, [int(e) for row in et for e in row])
    require(em.matdet() == 1, "The polarization is not principal.")
    ns, G = ns_general(FIELD, basis)
    columns = PARI.matrix(4, 4, [str(row[i]) for i in range(4) for row in ns])
    classes = PARI.matsolve(columns, PARI.matrix(4, 1, [str(x) for x in theta]))
    coords = [Fr(str(classes[i, 0])) for i in range(4)]
    require(all(x.denominator == 1 for x in coords), "The polarization is outside NS.")
    q = reduce_form(q_from_ns(G, [int(x) for x in coords]))
    return {
        "glued_basis": encode_basis(basis),
        "polarization_matrix": integer_matrix(em, 4),
        "humbert_form": q,
    }


def isometry(source, target):
    change = PARI.qfisom(source, target)
    require(
        change != 0
        and abs(change.matdet()) == 1
        and change.mattranspose() * target * change == source,
        "The computed form does not match the prescribed twin by an integral isometry.",
    )
    return integer_matrix(change, 3)


def incidence_targets():
    archive = read(HERE / "twin_lattices_-19.json")
    evidence = read(HERE.parent / "opposite-secant-spectra" / "certificate.json")
    level = next(row for row in evidence["level_three"] if row["D"] == -19)
    secants, = level["twin_pairs"]
    targets = {name: gram(archive[name]["q"]) for name in TWINS}
    changes = {
        name: isometry(targets[name], gram(secants[side]))
        for name, side in zip(TWINS, ("left", "right"))
    }
    counts = []
    for target in targets.values():
        half_counts = PARI.qfrep(target, 9)
        require(half_counts[0] == 0, "A twin unexpectedly represents one.")
        counts.append(2*int(half_counts[8]))
    require(counts == [0, 4], "The prescribed targets must have different degree-three counts.")
    return targets, changes, dict(zip(TWINS, counts))


def certify_equation(model, eta, candidates, target):
    s1, s2 = (eval_L(model[f"s{i}"], eta) for i in (1, 2))
    roots1, a1, labels1 = label_elliptic_quotient(s1, s2, candidates)
    roots2, a2, labels2 = label_elliptic_quotient(s2, s1, candidates)
    perm = reciprocal_permutation(roots1, roots2, labels1, labels2)
    data = humbert_from_gluing(a1, a2, perm)
    return {
        "E1_ideal_basis": encode_basis(a1),
        "E2_ideal_basis": encode_basis(a2),
        "half_period_permutation": perm,
        **data,
        "isometry_to_archived_twin": isometry(gram(data["humbert_form"]), target),
    }


def run():
    descent = check_descent()
    class_polynomial = descent["H_minus_171_descending"]
    require(
        list(map(int, PARI.Vec(PARI.polclass(-171)))) == class_polynomial,
        "The verified class polynomial differs from H_-171.",
    )
    L = poly(L_POLY)
    targets, changes, counts = incidence_targets()
    equations = read(HERE / "equations_-19_reduced.json")
    with ctx.workprec(PRECISION):
        candidates = ideal_candidates(class_polynomial)
        roots = checked_roots(list(reversed([int(c) for c in L.all_coeffs()])))
        embeddings = []
        for eta in roots:
            embeddings.append({
                name: certify_equation(equations[name], eta, candidates, targets[name])
                for name in TWINS
            })
    # Sort by the full records rather than by Arb root order.
    embeddings.sort(key=lambda row: json.dumps(row, sort_keys=True))
    return {
        "precision_bits": PRECISION,
        "conductor": 3,
        "embedding_count": len(embeddings),
        "archived_to_opposite_secant_isometries": changes,
        "degree_three_counts": counts,
        "embeddings": embeddings,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--verify", action="store_true", help="Compare to the committed certificate.")
    group.add_argument("--write-certificate", action="store_true", help="Replace the certificate with a fresh replay.")
    args = parser.parse_args()
    result = run()
    if args.verify:
        require(result == read(CERTIFICATE), "The replay differs from the committed certificate.")
    if args.write_certificate:
        CERTIFICATE.write_text(json.dumps(result, indent=2) + "\n")
    print(f"Checked both equations in all {result['embedding_count']} embeddings.")
    print("Exact geometric quotient counts at degree 3:", result["degree_three_counts"])


if __name__ == "__main__":
    main()
