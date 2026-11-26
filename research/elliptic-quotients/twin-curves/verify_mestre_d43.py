"""Check the six archived equations against the exact algebraic moduli triples.

This validates the field and sextic data directly. The original comparison of
those triples to the CM lattices comes from numerical recognition.
"""

from cypari2 import Pari

from certify_link import HERE, read, require
from igusa_exact import absolute_invariants

KEYS = ("J2^5/J10", "J2^3J4/J10", "J2^2J6/J10")


def replay():
    pari = Pari(128_000_000)
    models = read(HERE / "mestre_models_-43.json")
    recognized = read(HERE / "invariants_F_-43.json")
    require(len(models) == 6 and set(models) == {name for name in recognized if not recognized[name].get("u")},
            "The certificate must contain each of the six non-dihedral models.")
    for name, record in sorted(models.items()):
        value = recognized[name]
        F = pari(record["F"]).subst("x", "w")
        L = pari(record["model_field"])
        require(F == pari(value["F"]).subst("x", "w") and F.poldegree() == 4 and F.polisirreducible(),
                "The recognized invariant field is not a quartic field.")
        require(L.poldegree() == 8 and L.polisirreducible(), "The equation field is not a degree-eight field.")
        base = [pari.Mod(pari(value[key]).subst("x", "w"), F) for key in KEYS]
        exact = absolute_invariants(pari, [
            pari.Mod(pari(c.replace("theta", "w")), L) for c in record["model_constant_first"]
        ])
        # Rebuild the field identification from the first invariant. The
        # power basis is invertible because its minimal polynomial is quartic.
        first = pari.minpoly(base[0], "t")
        require(first.poldegree() == 4 and first == pari.minpoly(exact[0], "t"),
                "The equation has a different absolute invariant field generator.")
        basis = pari.matrix(4, 4, [
            pari.polcoef(pari.lift(base[0]**j), i, "w") for i in range(4) for j in range(4)
        ])
        for wanted, actual in zip(base[1:], exact[1:], strict=True):
            coordinates = basis.matsolve(pari.Col([pari.polcoef(pari.lift(wanted), i, "w") for i in range(4)]))
            require(actual == sum(c*exact[0]**i for i, c in enumerate(coordinates)),
                    f"The equation {name} fails the joint invariant comparison.")
        print(f"{name}: verified degree-eight equation and exact quartic invariant triple.", flush=True)


if __name__ == "__main__":
    replay()
