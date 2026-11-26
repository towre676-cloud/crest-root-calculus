"""Verify Richelot paths from the two certified bielliptic models."""

from igusa_exact import absolute_invariants
from mestre_obstruction_d43 import FIELD_POLYNOMIAL, PARI
from certify_link import HERE, read, require

X = PARI("x")


def invariants(rhs):
    return absolute_invariants(PARI, [rhs.polcoef(i, "x") for i in range(7)])


def quartic_signatures(models):
    """Encode each triple as a quartic root and two polynomials in that root."""
    field = PARI(FIELD_POLYNOMIAL)
    require(field.polisirreducible(), "The common coefficient field is reducible.")
    result = {}
    for name, model in models.items():
        values = [PARI.Mod(PARI(c).subst("x", "w"), field) for c in model["model_constant_first"]]
        first, second, third = absolute_invariants(PARI, values)
        minimal = PARI.minpoly(first, "t")
        require(minimal.poldegree() == 4, "The first invariant must generate the coefficient field.")
        basis = PARI.matrix(4, 4, [
            PARI.lift(first**j).polcoef(i, "w") for i in range(4) for j in range(4)
        ])
        expressions = [basis.matsolve(PARI.Col([PARI.lift(value).polcoef(i, "w") for i in range(4)]))
                       for value in (second, third)]
        result[name] = minimal, expressions
    return result


def matches(signature, triple):
    minimal, expressions = signature
    first = triple[0]
    return minimal.subst("t", first) == 0 and all(
        sum(coefficient*first**i for i, coefficient in enumerate(expression)) == value
        for expression, value in zip(expressions, triple[1:], strict=True)
    )


def verify_edge(record, source, target):
    field = PARI(record["field"])
    require(field.polisirreducible(), "The Richelot witness requires an irreducible field.")
    roots = [PARI.Mod(PARI(c), field) for c in record["roots"]]
    require(len(roots) == 6 and all(a != b for i, a in enumerate(roots) for b in roots[:i]),
            "The witness does not have six distinct roots.")
    pairing = record["pairs"]
    require(len(pairing) == 3 and all(len(pair) == 2 for pair in pairing) and
            sorted(i for pair in pairing for i in pair) == list(range(6)),
            "The witness does not partition the branch points.")
    factors = [(X - roots[i])*(X - roots[j]) for i, j in pairing]
    require(matches(source, invariants(factors[0]*factors[1]*factors[2])),
            "The Richelot source does not match the quartic model.")
    coefficient_matrix = PARI.matrix(3, 3, [
        q.polcoef(i, "x") for q in factors for i in range(3)
    ])
    require(coefficient_matrix.matdet() != 0, "The Richelot system is degenerate.")
    brackets = [a.deriv()*b - a*b.deriv() for a, b in
                zip(factors, factors[1:] + factors[:1], strict=True)]
    require(matches(target, invariants(brackets[0]*brackets[1]*brackets[2])),
            "The Richelot target does not match the quartic model.")


def verify_paths():
    sources = read(HERE / "models_certified_-43.json")
    models = read(HERE / "mestre_quartic_models_-43.json")
    signatures = quartic_signatures(sources | models)
    distances = {name: 0 for name in sources}
    records = read(HERE / "richelot_paths_-43.json")
    require(len(records) == len(models), "Expected one incoming Richelot edge for each quartic sextic.")
    for record in records:
        source, target = record["source"], record["target"]
        require(source in distances and target in models and target not in distances,
                "The Richelot witnesses do not form paths rooted at the bielliptic models.")
        verify_edge(record, signatures[source], signatures[target])
        distances[target] = distances[source] + 1
    return {name: distances[name] for name in models}


if __name__ == "__main__":
    for model, steps in verify_paths().items():
        print(f"{model}: {steps} exact Richelot steps from a certified CM model.")
