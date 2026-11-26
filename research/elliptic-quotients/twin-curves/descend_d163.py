"""Descend the two bielliptic D=-163 models to their common quartic field."""

import argparse
import json

from descend import from_dihedral, poly
from igusa_exact import absolute_invariants
from verify_d163 import HERE, PARI, read, require


def compute():
    initial = read(HERE / "models_certified_-163.json")
    field = PARI(initial["field"])
    result = {}
    for name, row in initial["models"].items():
        old_coeffs = [PARI.Mod(PARI(c), field) for c in row["model_constant_first"]]
        e1, e2, product = -old_coeffs[4], old_coeffs[2], -old_coeffs[0]
        u = e1**3/product + e2**3/product**2
        v = e1*e2/product
        # These are the dihedral invariants after scaling the cubic roots
        # by a cube root of their product.
        modulus = PARI.minpoly(u+v, "x")
        require(modulus.poldegree() == 4 and modulus.polisirreducible(),
                "The sum of the dihedral invariants must generate a quartic field.")
        reduced, embedding = PARI.polredbest(modulus, 1)
        th = PARI.lift(PARI.modreverse(embedding)).subst("x", u+v)
        require(reduced.subst("x", th) == 0, "The reduced generator does not satisfy its polynomial.")

        basis = PARI.matrix(24, 4, [PARI.lift(th**j).polcoef(i, "w") for i in range(24) for j in range(4)])
        desc = []
        for old in (u, v):
            coord = PARI.matinverseimage(basis, PARI.Col([PARI.lift(old).polcoef(i, "w") for i in range(24)]))
            require(len(coord) == 4, "The dihedral invariant does not descend to the quartic field.")
            exact = sum(c*PARI("x")**j for j, c in enumerate(coord))
            require(exact.subst("x", th) == old, "The dihedral-invariant embedding failed.")
            desc.append(poly(str(exact)))

        coeffs = from_dihedral(*desc, poly(str(reduced)))
        evals = [PARI(str(c.as_expr()).replace("**", "^")).subst("x", th) for c in coeffs]
        require(absolute_invariants(PARI, evals) == absolute_invariants(PARI, old_coeffs),
                "The original and quartic models have different geometric invariants.")
        result[name] = {
            "field": str(reduced),
            "embedding_in_degree_24_field": str(PARI.lift(th)),
            "u": str(desc[0].as_expr()),
            "v": str(desc[1].as_expr()),
            "model_constant_first": [str(c.as_expr()).replace("**", "^") for c in coeffs],
        }
    require(len({row["field"] for row in result.values()}) == 1, "The quartic fields differ.")
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args()
    data = compute()
    path = HERE / "quartic_models_-163.json"
    if args.write_certificate:
        path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")
    else:
        require(data == read(path), "The quartic-model certificate differs from the replay.")
    print("Verified both explicit quartic descents:", next(iter(data.values()))["field"])
