"""Verify quartic Mestre points and integral sextics."""

import argparse

from igusa_exact import absolute_invariants
from mestre_obstruction_d43 import INVARIANTS, PARI, conic_matrix
from verify_mestre_d43 import HERE, read, require


def replay(discriminant=-43):
    recognized = read(HERE / f"invariants_F_{discriminant}.json")
    models = read(HERE / f"mestre_quartic_models_{discriminant}.json")
    formulas = read(HERE / "mestre_formulas.json")
    field = PARI(next(iter(recognized.values()))["F"]).subst("x", "w")
    require(field.poldegree() == 4 and field.polisirreducible(), "Expected a quartic field.")
    nf = PARI.nfinit(field)
    require(len(models) == {-43: 6, -67: 12}[discriminant] and set(models) == {name for name in recognized if not recognized[name].get("u")},
            "The receipt must contain all the non-dihedral cases.")

    for name, model in sorted(models.items()):
        values = recognized[name]
        require(PARI(values["F"]).subst("x", "w") == field, f"{name}: invariant field mismatch.")
        point = PARI.Col([PARI.Mod(PARI(entry), field) for entry in model["conic_point"]])
        gram = conic_matrix(values, formulas, field)
        require(len(point) == 3 and any(point) and point.mattranspose()*gram*point == 0,
                f"{name}: the stored point is not on the Mestre conic.")

        coeffs = [PARI.Mod(PARI(entry), field) for entry in model["model_constant_first"]]
        require(len(coeffs) == 7 and coeffs[-1] != 0, f"{name}: expected a sextic.")
        require(all(c.denominator() == 1 for value in coeffs for c in PARI.nfalgtobasis(nf, value)),
                f"{name}: the sextic coefficients are not integral.")
        desired = [PARI.Mod(PARI(values[key]).subst("x", "w"), field) for key in INVARIANTS]
        require(list(absolute_invariants(PARI, coeffs)) == desired,
                f"{name}: the sextic does not have the target Igusa invariants.")
        print(f"{name}: quartic point and integral sextic verified.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discriminant", type=int, choices=(-43, -67), default=-43)
    replay(parser.parse_args().discriminant)
