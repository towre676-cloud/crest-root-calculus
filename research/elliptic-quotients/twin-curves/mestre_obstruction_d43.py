"""Decide the six D=-43 Mestre conics by exact local arithmetic."""

import argparse
import json
from math import prod

from cypari2 import Pari
from flint import fmpz

from certify_link import HERE, read, require

PARI = Pari(64_000_000)
INVARIANTS = ("J2^5/J10", "J2^3J4/J10", "J2^2J6/J10")
FIELD_POLYNOMIAL = "w^4-w^3+6*w^2-w+7"


def conic_matrix(invariants, formulas, field):
    values = [PARI.Mod(PARI(invariants[key]).subst("x", "w"), field) for key in INVARIANTS]
    z1, z2, z3 = values

    def evaluate(key):
        return PARI.substvec(PARI(formulas[key]), ["J2", "J4", "J6", "J10"],
                             [PARI.Mod(1, field), z2/z1, z3/z1, 1/z1])

    I2 = evaluate("I2")
    require(I2 != 0, "The Mestre normalization requires I2 nonzero.")
    return PARI.matrix(3, 3, [
        evaluate(f"A{min(i, j)}{max(i, j)}") / I2**(i+j+1)
        for i in range(1, 4) for j in range(1, 4)
    ])


def denominator(values):
    return PARI.lcm([PARI.denominator(c) for value in values for c in PARI.Vec(PARI.lift(value))])


def exact_prime_support(integer):
    """Factor by trial division and perfect powers before invoking PARI."""
    require(integer.denominator() == 1 and integer > 0, "The norm must be a positive integer.")
    norms = {}
    for base, power in fmpz(int(integer)).factor_smooth(23, proved=-1):
        if PARI.isprime(int(base)):
            factors = [(int(base), 1)]
        else:
            factors = [(int(q), int(e)) for q, e in zip(*PARI.factorint(int(base)), strict=True)]
        for p, e in factors:
            require(PARI.isprime(p), "A determinant factor is not a certified rational prime.")
            norms[p] = norms.get(p, 0) + int(power)*e
    require(prod(PARI(p)**e for p, e in norms.items()) == integer,
            "Incomplete factorization of the integral determinant norm.")
    return sorted(norms.items())


def local_representative(nf, value, prime, power):
    """Find a small integral representative modulo prime^power, via integer HNF."""
    coeffs = PARI.nfalgtobasis(nf, value)
    degree = len(coeffs)
    scale = PARI.lcm([PARI.denominator(c) for c in coeffs])
    ideal = PARI.idealpow(nf, prime, power+PARI.idealval(nf, scale, prime))
    augmented = PARI.matrix(degree, 2*degree, [
        (scale if i == j else 0) if j < degree else ideal[i, j-degree]
        for i in range(degree) for j in range(2*degree)
    ])
    lattice, transform = PARI.mathnf(augmented, 1)
    weights = lattice.matsolve(scale*coeffs)
    require(all(PARI.denominator(c) == 1 for c in weights), "The input is not locally integral.")
    bound = prime[0]**power
    lifted = PARI.Col([
        sum(transform[i, j+degree]*weights[j] for j in range(degree)) % bound
        for i in range(degree)
    ])
    representative = PARI.nfbasistoalg(nf, lifted)
    require(value == representative or PARI.idealval(nf, value-representative, prime) >= power,
            "The local unit representative does not have the required precision.")
    return representative


def local_symbol(nf, a, b, prime):
    p, e, f = (int(prime[i]) for i in (0, 2, 3))
    v, w = (int(PARI.idealval(nf, x, prime)) for x in (a, b))
    if p != 2 and v % 2 == w % 2 == 0:
        return v, w, 1

    uniformizer = PARI.nfbasistoalg(nf, prime[1])
    require(PARI.idealval(nf, uniformizer, prime) == 1, "The chosen element is not a uniformizer.")
    u, z = a/uniformizer**v, b/uniformizer**w
    if p == 2:
        # Every unit congruent to 1 modulo prime^(2e+1) is a local square.
        aa, bb = (local_representative(nf, x, prime, 2*e+1) for x in (u, z))
        symbol = PARI.nfhilbert(nf, aa*uniformizer**(v % 2), bb*uniformizer**(w % 2), prime)
        return v, w, int(symbol)

    reduction = PARI.nfmodprinit(nf, prime)
    res = PARI.nfmodpr(nf, u, reduction)**(w % 2)*PARI.nfmodpr(nf, z, reduction)**(v % 2)
    q = p**f
    require(res != 0, "A locally normalized unit has zero residue.")
    sign = -1 if v % 2 and w % 2 and q % 4 == 3 else 1
    return v, w, sign if res**((q-1)//2) == 1 else -sign


def compute():
    fields = read(HERE / "mestre_models_-43.json")
    invariants = read(HERE / "invariants_F_-43.json")
    formulas = read(HERE / "mestre_formulas.json")
    F = PARI(FIELD_POLYNOMIAL)
    require(F.polisirreducible(), "The field polynomial is reducible.")
    nf = PARI.nfinit(F)
    require(list(map(int, nf[1])) == [0, 2], "The quartic field is expected to be totally imaginary.")
    records = {}
    for name, archived in sorted(fields.items()):
        matrix = conic_matrix(invariants[name], formulas, F)
        d0, d1, d2 = (PARI.Mod(PARI(value), F) for value in archived["conic_diagonal"])
        det = matrix.matdet()
        require(d0*d1*d2 != 0 and matrix[0, 0] == d0
                and matrix[0, 0]*matrix[1, 1]-matrix[0, 1]**2 == d0*d1 and det == d0*d1*d2,
                "The diagonal coefficients do not agree with the conic's principal minors.")
        scale = denominator(value for column in matrix for value in column)
        support = exact_prime_support(abs(PARI.norm(scale**3*det)))
        locals_ = []
        for rational_prime in sorted({2, *(p for p, _ in support)}):
            for prime in PARI.idealprimedec(nf, rational_prime):
                v, w, symbol = local_symbol(nf, -d1/d0, -d2/d0, prime)
                locals_.append([rational_prime, [int(c) for c in prime[1]],
                                int(prime[2]), int(prime[3]), v, w, symbol])
        records[name] = {"integral_scale": str(scale), "determinant_norm_factors": support,
                         "local_symbols": locals_, "has_quartic_point": all(row[-1] == 1 for row in locals_)}
    return json.loads(json.dumps(records))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args()
    result = compute()
    path = HERE / "mestre_obstruction_-43.json"
    if args.write_certificate:
        # Keep each prime-ideal record on one line.
        rows = ["{"]
        for name, record in result.items():
            rows.extend([
                f'  "{name}": {{',
                '    "integral_scale": ' + json.dumps(record["integral_scale"]) + ",",
                '    "determinant_norm_factors": ' + json.dumps(record["determinant_norm_factors"]) + ",",
                '    "local_symbols": [',
                ",\n".join("      " + json.dumps(place) for place in record["local_symbols"]),
                "    ],",
                '    "has_quartic_point": ' + json.dumps(record["has_quartic_point"]),
                "  },",
            ])
        rows[-1] = rows[-1].rstrip(",")
        rows.append("}")
        path.write_text("\n".join(rows) + "\n", encoding="utf-8")
    else:
        require(result == read(path), "The local conic replay differs from the receipt.")
    for name, row in result.items():
        print(f"{name}: {len(row['local_symbols'])} finite places; quartic point exists: {row['has_quartic_point']}.")


if __name__ == "__main__":
    main()
