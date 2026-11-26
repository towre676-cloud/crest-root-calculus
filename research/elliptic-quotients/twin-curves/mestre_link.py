"""Recover exact geometric endomorphism orders of the quartic Mestre models."""

import argparse
import itertools
import json
import sys
from fractions import Fraction
from pathlib import Path

from cypari2 import Pari
from flint import acb, acb_mat, acb_poly, arb, arb_mat, ctx

from certify_link import HERE, read, require
from mestre_obstruction_d43 import FIELD_POLYNOMIAL

sys.path.insert(0, str(HERE.parent / "rigorous-search"))
sys.path.insert(0, str(HERE.parent / "inverse-spectrum"))

import periods_arb as PA
from enumeration import short_vectors
from gluing import q_from_ns, reduce_form, isometric
from verify_d43 import compute as verify_p3, pair_spectrum
from search import endomorphism_form, frame, integral_lower_bound

pari = Pari(256_000_000)


def arb_q(c):
    c = pari(c)
    return arb(int(pari.numerator(c))) / int(pari.denominator(c))


def eval_poly(expr, var, x):
    P = pari(expr)
    s = acb(0)
    for k in range(int(pari.poldegree(P, var)), -1, -1):
        s = s * x + acb(arb_q(pari.polcoef(P, k, var)))
    return s


def normalised_quintic(model, field, emb_index, hi_bits):
    ctx.prec = hi_bits
    Kp = pari(field)
    var = "w"
    kc = [acb(arb_q(pari.polcoef(Kp, k, var))) for k in range(int(pari.poldegree(Kp, var)) + 1)]
    roots = acb_poly(kc).roots(tol=arb(2) ** -(hi_bits - 200), maxprec=4 * hi_bits)
    require(len(roots) == Kp.poldegree() and all(not a.overlaps(b)
                for i, a in enumerate(roots) for b in roots[:i]), "The field embeddings are not isolated.")
    w0 = roots[emb_index]
    co = [eval_poly(c.replace("theta", "w"), var, w0) for c in model]
    br = acb_poly(co).roots(tol=arb(2) ** -(hi_bits // 2), maxprec=4 * hi_bits)
    require(len(br) == 6, "The sextic must have six isolated branch points.")
    best = None
    for a, b, c in itertools.permutations(range(6), 3):
        r1, r2, r3 = br[a], br[b], br[c]
        T = lambda x: (x - r1) * (r2 - r3) / ((x - r3) * (r2 - r1))
        rest = [T(br[k]) for k in range(6) if k not in (a, b, c)]
        pts = [complex(0), complex(1)] + [complex(float(z.real.mid()), float(z.imag.mid())) for z in rest]
        sep = min(abs(p - q) for i, p in enumerate(pts) for q in pts[:i])
        span = max(abs(p) for p in pts)
        score = sep / (1 + span)
        if best is None or score > best[0]:
            best = (score, rest)
    lam = best[1]
    # y^2 = x (x - 1) (x - l1)(x - l2)(x - l3): coefficients, low to high
    P = acb_poly([acb(0), acb(1)]) * acb_poly([acb(-1), acb(1)])
    for l in lam:
        P = P * acb_poly([-l, acb(1)])
    f = [P[k] for k in range(6)]
    return f


def period_matrix_acb(f, bits):
    PA.set_prec(bits)
    best = None
    for z0 in PA.Z0_CHOICES:
        g, z = PA.moved_poly(f, z0)
        if g[6].contains(0):
            continue
        em = [PA.mid(r) for r in acb_poly(g).roots(tol=arb(2) ** -60, maxprec=4000)]
        span = max(abs(a - b) for a in em for b in em)
        cl = max(PA.clearance(em, i) for i in range(6)) / span
        if best is None or cl > best[0]:
            best = (cl, z0)
    require(best is not None, "No Möbius shift away from the branch locus.")
    g, z = PA.moved_poly(f, best[1])
    S = PA.star_periods_arb(g)
    P = [[-w1, -(z * w1 + w0)] for w0, w1 in S["periods"]]
    S["Pi"] = [[P[c][r] for c in range(4)] for r in range(2)]
    S["E"] = [row[:4] for row in S["intersection"][:4]]
    Pm, Em = acb_mat(S["Pi"]), acb_mat(S["E"])
    rie = Pm * Em.inv() * Pm.transpose()
    require(all(rie[i, j].contains(0) for i in range(2) for j in range(2)), "Riemann relation failed.")
    return S


def flat(matrix):
    return [int(matrix[i, j]) for i in range(4) for j in range(4)]


def hnf(matrices):
    values = list(matrices)
    return pari.mathnf(pari.matrix(16, len(values), [
        values[j][i] for i in range(16) for j in range(len(values))
    ]))


def columns_to_lists(matrix):
    rows, cols = map(int, matrix.matsize())
    return [[int(matrix[i, j]) for i in range(rows)] for j in range(cols)]


def integral_coordinates(lattice, vector):
    coefficients = pari.matinverseimage(lattice, pari.Col(vector))
    require(len(coefficients) == len(lattice) and all(c.denominator() == 1 for c in coefficients),
            "The candidate order is not closed under the requested operation.")
    return coefficients


def endomorphism_shortlist(data, trace_bound):
    """Return every integral matrix passing the bounded enclosure test."""
    J, G = frame(data)
    scale = 2**32
    form = integral_lower_bound(endomorphism_form(G, J, arb(2)**20), scale)
    exact = pari.matrix(16, 16, [entry for row in form for entry in row])
    transform = exact.qflllgram()
    require(abs(transform.matdet()) == 1, "The enumeration basis change is not unimodular.")
    reduced = transform.mattranspose() * exact * transform
    reduced_rows = columns_to_lists(reduced)  # The form is symmetric.
    candidates = 0
    survivors = []
    for vector in short_vectors(reduced_rows, trace_bound * scale):
        values = list(map(int, transform * pari.Col(vector)))
        R = arb_mat([values[i:i+4] for i in range(0, 16, 4)])
        commutator = J * R - R * J
        candidates += 1
        if all(commutator[i, j].contains(0) for i in range(4) for j in range(4)):
            survivors.append(values)
    return hnf(survivors), candidates, len(survivors)


def intersection_form(data):
    E = pari.matrix(4, 4, [e for row in data["E"] for e in row])
    require(E.mattranspose() == -E and abs(E.matdet()) == 1,
            "The intersection form is not alternating and unimodular.")
    return E


def humbert_from_order(data, lattice):
    E = intersection_form(data)
    inverse = E**-1
    rank = len(lattice)
    matrices = [pari.matrix(4, 4, v) for v in columns_to_lists(lattice)]
    adjoints = [inverse * R.mattranspose() * E for R in matrices]
    dag = [integral_coordinates(lattice, flat(A)) for A in adjoints]
    involution = pari.matrix(rank, rank, [dag[j][i] for i in range(rank) for j in range(rank)])
    sym = lattice * pari.matkerint(involution - pari.matid(rank))
    rho = len(sym)
    require(rho == 4, "Expected four symmetric classes.")
    ns_matrices = [pari.matrix(4, 4, v) for v in columns_to_lists(sym)]
    theta = list(map(int, integral_coordinates(sym, flat(pari.matid(4)))))
    ns = [[Fraction(int(a.trace()*b.trace()), 4) - Fraction(int((a*b).trace()), 2)
           for b in ns_matrices] for a in ns_matrices]
    require(all(entry.denominator == 1 for row in ns for entry in row), "Nonintegral intersection form.")
    ns = [[int(entry) for entry in row] for row in ns]
    return reduce_form(q_from_ns(ns, theta))


def certify_cm_order(data, multiplier, discriminant=-43):
    """Isolate the integral scalar-CM matrix, with integrality from the isogeny chain."""
    J, _ = frame(data)
    enclosure = J * (multiplier * arb(-discriminant).sqrt())
    entries = [enclosure[i, j].unique_fmpz() for i in range(4) for j in range(4)]
    require(all(entry is not None for entry in entries), "The scalar-CM enclosure is inconclusive.")
    C = pari.matrix(4, 4, list(map(int, entries)))
    require(C*C == discriminant * multiplier**2 * pari.matid(4), "The scalar-CM matrix has the wrong square.")
    units = [pari.matrix(4, 4, [int(i == j) for i in range(16)]) for j in range(16)]
    commutators = [flat(C*unit - unit*C) for unit in units]
    linear_map = pari.matrix(16, 16, [commutators[j][i] for i in range(16) for j in range(16)])
    order = pari.mathnf(pari.matkerint(linear_map))
    require(len(order) == 8, "Expected rank eight for a CM square.")
    return C, order


def bounded_comparison(data, certified, bound=128):
    lattice, candidates, survivors = endomorphism_shortlist(data, bound)
    require(lattice == certified, "The bounded shortlist does not span the certified order.")
    return {"trace_bound": bound, "enumerated_up_to_sign": candidates, "survivors_up_to_sign": survivors}


def replay(embeddings, names, compare_shortlist=False):
    from verify_richelot_d43 import verify_paths

    require(verify_p3() == read(HERE / "certified_link_-43_p3.json"), "The p3 anchor certificate changed.")
    distances = verify_paths()
    models = read(HERE / "mestre_quartic_models_-43.json")
    archive = read(HERE / "twin_lattices_-43.json")
    results = {}
    for name in names or sorted(models):
        for embedding in embeddings:
            f = normalised_quintic(models[name]["model_constant_first"], FIELD_POLYNOMIAL, embedding, 1800)
            data = period_matrix_acb(f, 600)
            multiplier = 6 * 2**distances[name]
            matrix, order = certify_cm_order(data, multiplier)
            result = {"richelot_steps": distances[name], "cm_multiplier": multiplier,
                      "scalar_cm_matrix": columns_to_lists(matrix.mattranspose()),
                      "rank_endomorphisms": len(order), "humbert_form": humbert_from_order(data, order)}
            if compare_shortlist:
                result.update(bounded_comparison(data, order))
            matches = [label for label, entry in sorted(archive.items())
                       if isometric(result["humbert_form"], entry["q"])]
            require(matches == [name], f"{name}, embedding {embedding}: matches {matches}.")
            results.setdefault(name, {})[str(embedding)] = result
            print(f"{name}, embedding {embedding}: certified form isometric to {name}.", flush=True)
    spectra = {}
    for pair in ("p0", "p1", "p2"):
        names = tuple(f"{pair}_q{index}" for index in (1, 2))
        if all(name in results for name in names):
            targets = {name: pari.matrix(3, 3, [entry for row in next(iter(results[name].values()))["humbert_form"]
                                                 for entry in row]) for name in names}
            spectra[pair] = pair_spectrum(targets, names)
            first = spectra[pair]["first_degree"]
            counts = [spectra[pair]["primitive_counts"][name][-1] for name in names]
            print(f"{pair}: first distinguishing degree {first}, counts {counts}.")
    return {"models": results, "pair_spectra": spectra}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--embedding", type=int, choices=range(4), action="append")
    parser.add_argument("--model", action="append")
    parser.add_argument("--output", type=Path)
    parser.add_argument("--compare-shortlist", action="store_true")
    args = parser.parse_args()
    result = replay(args.embedding or range(4), args.model, args.compare_shortlist)
    if args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
