"""Verify explicit level-three twins at discriminant -163 via CM torsion."""

import argparse
from itertools import permutations
import json

from cypari2 import Pari
from flint import acb, ctx

from certify_link import (
    HERE, checked_roots, classify, eval_L, gram, half_period_labels,
    humbert_from_gluing, ideal_candidates, isometry, read, require,
)
from verify_d43 import pair_spectrum

D = -163
PARI = Pari(128_000_000)
X = PARI("x")
A = PARI(-2174420)
B = PARI(1234136692) + PARI(1)/4
F = X**3 + A*X + B
PSI3 = 3*X**4 + 6*A*X**2 + 12*B*X - A**2


def targets():
    data = read(HERE.parent / "walk-to-products" / "collision_counts.json")
    return {f"p{i}_{side}": gram(record[side])
            for i, record in enumerate(row for row in data if row["d"] == D and row["l"] == 3)
            for side in ("q1", "q2")}


def torsion_data(record):
    field = PARI(record["field"])
    require(field.poldegree() == 24 and field.polisirreducible(), "The equation field must have degree 24.")
    points, kernels = ([PARI.Mod(PARI(value), field) for value in record[key]]
                       for key in ("torsion_2", "kernel_3"))
    require(len(set(map(str, points))) == len(points) == 3 and all(F.subst("x", t) == 0 for t in points),
            "The three points must be all the 2-torsion abscissae.")
    require(len(set(map(str, kernels))) == len(kernels) == 4 and all(PSI3.subst("x", t) == 0 for t in kernels),
            "The four points must be all the cyclic 3-kernel abscissae.")
    return field, points, kernels


def elliptic_data(torsion, kernels):
    polynomial = PARI.polclass(9*D)
    data = []
    for t in kernels:
        tangent = 2*(3*t**2 + A)
        a = A-5*tangent
        b = B-7*(4*F.subst("x", t)+t*tangent)
        roots = [u + tangent/(u-t) + 4*F.subst("x", t)/(u-t)**2 for u in torsion]
        require(len(set(map(str, roots))) == 3 and all(r**3+a*r+b == 0 for r in roots),
                "Vélu abscissae do not give three distinct roots of the new cubic.")
        require(4*a**3+27*b**2 != 0, "Singular elliptic factor.")
        j = 1728*4*a**3/(4*a**3+27*b**2)
        require(polynomial.subst("x", j) == 0, "The quotient does not have the asserted CM.")
        data.append((roots, a, b))
    return list(map(int, PARI.Vec(polynomial))), data


def even_sextic(left, right):
    """Return y²=∏(x²-u_i), with ordered u_i matched to the right-hand roots."""
    T = PARI.matrix(3, 4, [v for a, b in zip(left, right, strict=True) for v in (a, 1, -a*b, -b)])
    basis = T.matker()
    require(len(basis) == 1, "The projective root match is not unique.")
    a, b, c, d = basis[0]
    require(c != 0 and a*d-b*c != 0, "The match must have a finite pole and nonzero determinant.")
    u = [v+d/c for v in left]
    require(all(u) and all(w == a/c + (b*c-a*d)/(c*c*v) for w, v in zip(right, u, strict=True)),
            "The reciprocal root match fails in the equation field.")
    return [(X**2-v) for v in u]


def embeddings(field):
    return checked_roots([int(field.polcoef(i, "w")) for i in range(int(field.poldegree())+1)])


def elliptic_labels(data, embedding, candidates):
    def evaluate(t):
        return eval_L(str(PARI.lift(t)).replace("w", "x"), embedding)
    roots, a, b = data
    g2, g3 = -4*evaluate(a), -4*evaluate(b)
    ideal, tau = classify(g2, g3, candidates)
    labels = half_period_labels([evaluate(root) for root in roots], acb(0), g2, g3, tau)
    return ideal, labels


def form_for_match(left, right, pairing):
    ideal1, labels1 = left
    ideal2, labels2 = right
    perm = [None]*3
    for i in range(3):
        perm[labels1[i]] = labels2[pairing[i]]
    return humbert_from_gluing(ideal1, ideal2, perm, D)["humbert_form"]


def compute():
    field, torsion, kernels = torsion_data(read(HERE / "d163_torsion_field.json"))
    polynomial, elliptic = elliptic_data(torsion, kernels)
    archive = targets()
    names = ("p16_q1", "p16_q2")
    choices = {}
    with ctx.workprec(2600):
        roots = embeddings(field)
        candidates = ideal_candidates(polynomial, D)
        for i, left in enumerate(elliptic):
            for j, right in enumerate(elliptic[i+1:], i+1):
                L = elliptic_labels(left, roots[0], candidates)
                R = elliptic_labels(right, roots[0], candidates)
                for permutation in permutations(range(3)):
                    form = gram(form_for_match(L, R, permutation))
                    matches = [name for name in names if PARI.qfisom(form, archive[name])]
                    if not matches:
                        continue
                    name, = matches
                    if name in choices:
                        continue
                    right_roots = [right[0][k] for k in permutation]
                    factors = even_sextic(left[0], right_roots)
                    model = factors[0]*factors[1]*factors[2]
                    for root in roots:
                        l = elliptic_labels(left, root, candidates)
                        r = elliptic_labels(right, root, candidates)
                        q = gram(form_for_match(l, r, permutation))
                        isometry(q, archive[name])
                    choices[name] = {
                        "elliptic_kernel_indices": [i, j],
                        "torsion_permutation": list(permutation),
                        "model_constant_first": [str(PARI.lift(model.polcoef(k, "x"))) for k in range(7)],
                        "verified_embeddings": len(roots),
                        "humbert_form": str(form),
                    }
    require(set(choices) == set(names), "The search did not find both bielliptic twins.")
    return {"field": str(field), "ring_class_polynomial": polynomial, "models": choices,
            "spectrum": pair_spectrum({n: archive[n] for n in names}, names, discriminant=D)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args()
    data = compute()
    path = HERE / "models_certified_-163.json"
    if args.write_certificate:
        path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    else:
        require(data == read(path), "The stored certificate differs from the replay.")
    print("Verified D=-163 models:", ", ".join(sorted(data["models"])))
    print("First distinguishing degree:", data["spectrum"]["first_degree"])
    print("Primitive counts:", data["spectrum"]["primitive_counts"])
