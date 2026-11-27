"""Descend the determinant-1440 sextics by an exact projective norm witness."""

import argparse
import json

from exact_models import HERE, MODULI_POLYNOMIAL, absolute_invariants, field, pari, require


def coefficients(poly, variable, size):
    return [pari.polcoef(poly, k, variable) for k in range(size)]


def serialize(value):
    return str(pari.lift(value))


def transform_sextic(f, matrix, x):
    a, b, c, d = (matrix[i, j] for i in range(2) for j in range(2))
    return sum(coef * (a * x + b)**i * (c * x + d)**(6-i)
               for i, coef in enumerate(coefficients(f, x, 7)))


def normalized(f, x):
    return f / pari.polcoef(f, pari.poldegree(f, x), x)


def descend(model, field_data, witness):
    x, y, w = (pari(v) for v in ("x", "y", "w"))
    modulus = field_data["Pabs"]

    def parse(value):
        return pari.Mod(pari(value.replace("theta", "y")), modulus)

    def polynomial(values):
        return pari.Polrev([parse(value) for value in values], x)

    sigma = parse(witness["conjugation"])
    require(sigma != pari.Mod(y, modulus) and not pari.subst(modulus, y, sigma),
            "The proposed conjugation is not a nonidentity automorphism.")

    def conjugate(value):
        return pari.subst(pari.lift(value), y, sigma)

    require(conjugate(sigma) == pari.Mod(y, modulus), "Conjugation must have order two.")
    g = pari(MODULI_POLYNOMIAL).subst(x, w)
    generator = parse(witness["octic_generator_in_ring_class_field"])
    require(pari.polisirreducible(g) and not pari.subst(g, w, generator)
            and conjugate(generator) == generator, "The octic field must be the fixed field.")
    rad = field_data["s"]
    require(conjugate(rad) == -rad, "The extension is not generated over the fixed field by sqrt(-10).")

    f = polynomial(model["sextic_constant_first"])
    inv = absolute_invariants(pari, coefficients(f, x, 7))
    require(pari.minpoly(inv[0], x).poldegree() == 8 and all(conjugate(z) == z for z in inv),
            "The fixed field does not equal the field of moduli.")
    j1, j2 = (parse(model[k]["j"]) for k in ("phi1", "phi2"))
    require(j1 != j2 and conjugate(j1) == j2, "Conjugation must exchange the two labels.")

    # The unmarked critical point r of P2/R maps to infinity in the conjugate chart.
    # The other point s in that fibre maps to zero.
    R = polynomial(model["R"])
    P2 = polynomial(model["phi2"]["P2"])
    S2 = polynomial(model["phi2"]["S"])
    wronskian = pari.deriv(P2, x) * R - P2 * pari.deriv(R, x)
    line, remainder = pari.divrem(wronskian, S2, x)
    require(not remainder and pari.poldegree(line, x) == 1, "The residual critical point is not determined.")
    r = -pari.polcoef(line, 0, x)/pari.polcoef(line, 1, x)
    numerator = pari.subst(R, x, r) * P2 - pari.subst(P2, x, r) * R
    line, remainder = pari.divrem(numerator, (x-r)**2, x)
    require(not remainder and pari.poldegree(line, x) == 1, "The residual fibre has the wrong shape.")
    s = -pari.polcoef(line, 0, x)/pari.polcoef(line, 1, x)
    scale = -1/(3-(r-s)*pari.subst(pari.deriv(R, x), x, r)/pari.subst(R, x, r))
    entries = (scale, -scale*s, pari.Mod(1, modulus), -r)
    matrix = pari.matrix(2, 2, entries)
    conj_matrix = pari.matrix(2, 2, [conjugate(z) for z in entries])
    product = conj_matrix * matrix
    alpha = product[0, 0]
    require(matrix.matdet() and product == alpha * pari.matid(2), "The projectivity is not a cocycle.")
    conjugate_f = pari.Polrev([conjugate(z) for z in coefficients(f, x, 7)], x)
    require(normalized(transform_sextic(conjugate_f, matrix, x), x) == normalized(f, x),
            "The projectivity does not identify the branch divisors.")

    norm_witness = parse(witness["norm_witness"])
    require(norm_witness * conjugate(norm_witness) == alpha, "The relative norm identity failed.")
    linear_cocycle = matrix / norm_witness
    inverse_cocycle = linear_cocycle**-1
    for shift in (pari.Mod(1, modulus), rad, 1 + rad):
        change = shift * pari.matid(2) + conjugate(shift) * inverse_cocycle
        if change.matdet():
            break
    else:
        raise ArithmeticError("No invertible Hilbert 90 matrix in the prescribed choices.")
    conjugate_change = pari.matrix(2, 2, [conjugate(change[i,j]) for i in range(2) for j in range(2)])
    require(linear_cocycle * change == conjugate_change, "Hilbert 90 matrix identity failed.")
    new = normalized(transform_sextic(f, change, x), x)
    require(pari.poldegree(new, x) == 6 and pari.poldegree(pari.gcd(new, pari.deriv(new, x)), x) == 0,
            "The descended polynomial is not a squarefree sextic.")

    def vector(z):
        return pari.Col(coefficients(pari.lift(z), y, 16))

    basis = pari.Mat([vector(generator**j) for j in range(8)])
    output = []
    for v in coefficients(new, x, 7):
        coordinates = basis.matinverseimage(vector(v))
        require(len(coordinates) == 8, "A descended coefficient is outside the fixed field.")
        output.append(serialize(pari.Polrev(list(coordinates), w)))
    descended = [pari.Mod(pari(c), g) for c in output]
    require([pari.subst(pari.lift(c), w, generator) for c in descended] == coefficients(new, x, 7),
            "Embedding the octic model failed to recover the transformed sextic.")
    return {"field_polynomial_w": str(g), "sextic_constant_first": output}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true")
    args = parser.parse_args()
    witnesses = json.loads((HERE / "octic-descent-witnesses.json").read_text())
    models = json.loads((HERE / "models_1440.json").read_text())
    data = field()
    result = {name: descend(models[name], data, witnesses[name]) for name in ("A", "B")}
    archive = HERE / "models_1440_octic.json"
    if args.verify:
        require(result == json.loads(archive.read_text()), "The descended sextics differ from the receipt.")
        print("Both octic models and their exact isomorphisms verified.")
    else:
        archive.write_text(json.dumps(result, indent=2) + "\n")


if __name__ == "__main__":
    main()
