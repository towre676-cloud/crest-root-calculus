"""Replay the exact discriminant -43 bielliptic lattice identification."""

import argparse
from itertools import combinations
import json

from flint import acb, ctx

from certify_link import (
    HERE, PARI, certify_gluing, checked_roots, eval_L, gram, ideal_candidates,
    isometry, read, require,
)
from descend import elliptic_j_pair, from_dihedral, multiplication_polynomial, poly
from igusa_exact import absolute_invariants

DISCRIMINANT = -43
FIELD_POLYNOMIAL = "x^4-x^3+6*x^2-x+7"
TWINS = ("p3_q1", "p3_q2")
PRECISION = 3000
PARI.allocatemem(64_000_000)


def check_models(moduli, models, F, class_polynomial):
    result = {}
    for name in TWINS:
        model, recognized = models[name], moduli[name]
        require(poly(model["F"]) == F == poly(recognized["F"]), "The model fields differ.")
        u, v = (poly(model[key]).rem(F) for key in ("u", "v"))
        require(all(poly(recognized[key]).rem(F) == value for key, value in (("u", u), ("v", v))),
                "The dihedral invariants differ from the input record.")
        coefficients = from_dihedral(u, v, F)
        require(coefficients == [poly(c).rem(F) for c in model["model_constant_first"]],
                "The stored quartic-field equation fails the exact descent formula.")
        disc = (v*v-4*u+18*v-27).rem(F)
        require(bool(disc), "The even sextic must be smooth.")
        js = elliptic_j_pair(u, v, disc, F, class_polynomial)
        mod = PARI(str(F.as_expr()).replace("**", "^")).subst("x", "w")
        abs_invariants = absolute_invariants(PARI, [
            PARI.Mod(PARI(str(c.as_expr()).replace("**", "^")).subst("x", "w"), mod)
            for c in coefficients
        ])
        for key, value in zip(("J2^5/J10", "J2^3J4/J10", "J2^2J6/J10"), abs_invariants, strict=True):
            target = PARI.Mod(PARI(recognized[key]).subst("x", "w"), mod)
            require(value == target, f"The exact sextic has a different invariant {key}.")
        z = poly(recognized["J2^5/J10"]).rem(F)
        minimal = multiplication_polynomial(z, F)
        require(minimal.degree() == 4 and minimal.is_irreducible,
                "The absolute invariant does not generate the quartic field.")
        result[name] = {"quotient_j_sum_and_product": js, "moduli_generator_polynomial": str(minimal.as_expr())}
    return result


def pair_spectrum(targets):
    data = read(HERE.parent / "opposite-secant-spectra" / "certificate.json")
    level = next(row for row in data["level_three"] if row["D"] == DISCRIMINANT)
    pairs = []
    for pair in level["twin_pairs"]:
        left, right = (gram(pair[key]) for key in ("left", "right"))
        for a, b in ((left, right), (right, left)):
            if PARI.qfisom(targets[TWINS[0]], a) and PARI.qfisom(targets[TWINS[1]], b):
                pairs.append((pair, a, b))
    require(len(pairs) == 1, "The archived forms are not a unique opposite-secant pair.")
    pair, a, b = pairs[0]
    first = pair["first_degree"]
    transitions = [isometry(targets[name], form) for name, form in zip(TWINS, (a, b), strict=True)]
    counts = {}
    for name in TWINS:
        reps = PARI.qfrep(targets[name], first**2)
        primitive = []
        for m in range(1, first+1):
            primitive.append(2*int(reps[m*m-1]) - sum(primitive[j-1] for j in range(1, m) if m % j == 0))
        counts[name] = primitive
    left, right = (counts[name] for name in TWINS)
    require(left[:-1] == right[:-1] and left[-1] != right[-1], "The first distinguishing degree changed.")
    return {"first_degree": first, "primitive_counts": counts, "isometries_to_opposite_secants": transitions}


def compute():
    F = poly(FIELD_POLYNOMIAL)
    require(F.is_irreducible, "The quartic field polynomial must be irreducible.")
    class_polynomial = list(map(int, PARI.Vec(PARI.polclass(9*DISCRIMINANT))))
    models = read(HERE / "models_certified_-43.json")
    moduli = read(HERE / "invariants_F_-43.json")
    records = check_models(moduli, models, F, class_polynomial)
    archive = read(HERE / "twin_lattices_-43.json")
    targets = {name: gram(archive[name]["q"]) for name in models}
    for name in TWINS:
        records[name]["embeddings"] = []

    with ctx.workprec(PRECISION):
        candidates = ideal_candidates(class_polynomial, DISCRIMINANT)
        for th in checked_roots([int(F.nth(i)) for i in range(5)]):
            for name in TWINS:
                u, v = (eval_L(models[name][key], th) for key in ("u", "v"))
                # The exact parameters satisfy s1^3+s2^3=u and s1*s2=v
                # by construction. Any choice gives a geometrically isomorphic sextic.
                s1 = ((u+(u*u-4*v**3).sqrt())/2)**(acb(1)/3)
                require(s1.is_finite() and not s1.contains(0), "Could not enclose a nonzero dihedral parameter.")
                glued = certify_gluing(s1, v/s1, candidates, targets[name], DISCRIMINANT)
                form = gram(glued["humbert_form"])
                matches = sorted(other for other in archive if PARI.qfisom(form, gram(archive[other]["q"])))
                require(matches == [name], "The model must match exactly its recorded twin.")
                records[name]["embeddings"].append(glued)
    for name in TWINS:
        records[name]["embeddings"].sort(key=lambda row: json.dumps(row, sort_keys=True))
    return {"field": FIELD_POLYNOMIAL, "ring_class_polynomial": class_polynomial,
            "models": records, "pair": pair_spectrum(targets)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args()
    result = compute()
    path = HERE / "certified_link_-43_p3.json"
    if args.write_certificate:
        path.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    else:
        require(result == read(path), "The p3 replay differs from the committed certificate.")
    print("Certified both p3 equations in all four quartic embeddings.")
    print("First distinguishing degree:", result["pair"]["first_degree"])
    print("Geometric quotient counts:", result["pair"]["primitive_counts"])


if __name__ == "__main__":
    main()
