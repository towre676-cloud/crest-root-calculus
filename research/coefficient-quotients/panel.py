"""Generate isolated coefficient inputs and a separate provenance sidecar."""

import argparse
from collections import Counter
import json
from pathlib import Path
import random
import sys

import sympy as sp

from discover import X, coefficients, normalization

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "periodic-spectra"))
from spectral_curve import cell_discriminant, trace_and_quotient


def poly(expression):
    return sp.Poly(expression, X, domain=sp.QQ)


def dense_random(rng, degree, bound):
    return poly(X**degree + sum(rng.randint(-bound, bound) * X**i for i in range(degree)))


def affine(polynomial, scale, shift):
    return polynomial.compose(poly(scale * X + shift))


def composite(rng, d, m, bound, family):
    for _ in range(50):
        if family == "dense":
            inner = dense_random(rng, d, bound)
        else:
            parameter = rng.randint(1, bound)
            inner = poly(2 * sp.chebyshevt(d, X / 2)).compose(poly(X + parameter))
        outer = dense_random(rng, m, bound)
        value = outer.compose(inner)
        if value.gcd(value.diff()).degree() == 0:
            return value, {"construction": family, "inner": coefficients(inner),
                           "outer": coefficients(outer), "coefficient_bound": bound}
    raise ArithmeticError("Could not generate a smooth composite within the frozen attempts.")


def variable_fixtures(rng):
    u = poly(X**3 - 3 * X)
    used = set()
    fixtures = []
    for _ in range(3):
        a = rng.choice([n for n in range(-12, 13) if n not in used and abs(n) > 2])
        used.add(a)
        fixtures.append((u, poly((X**2 - 4) * (X - a)), poly(X**2 - 1), "cubic-saturated"))
        fixtures.append((u, poly((X + 2) * (X**2 + abs(a))), poly(X - 1), "cubic-proper"))
    for a, b in ((1, -2), (2, -3), (1, 4)):
        u = poly(sp.integrate(5 * X**2 * (X - a) * (X - b), X))
        values = [u.eval(c) for c in (0, a, b)]
        if len(set(values)) != len(values):
            raise ArithmeticError("Coincident branch values in a proper complementary fixture.")
        fixtures.append((u, poly(sp.prod(X - c for c in values)), poly(X * (X - a) * (X - b)),
                         "quintic-proper"))
    u = poly(2 * sp.chebyshevt(5, X / 2))
    fixtures.append((u, poly((X**2 - 4) * (X - 3)), poly(u.diff().as_expr() / 5),
                     "quintic-saturated"))
    return fixtures


def generate(seed):
    rng = random.Random(seed)
    cases = []

    def add(source, study, cohort, detail):
        if not source:
            raise ArithmeticError("Invalid generated radicand.")
        cases.append(({"coefficients": coefficients(source)}, {"study": study, "cohort": cohort, **detail}))

    plans = [(2, 3, 3), (3, 3, 10), (2, 6, 100), (4, 5, 4), (6, 6, 16), (9, 5, 7),
             (5, 7, 2**12), (7, 8, 2**18), (12, 6, 2), (16, 3, 2**20)]
    for index, (d, m, bound) in enumerate(plans):
        family = "dense" if index % 3 else "chebyshev-inner"
        source, detail = composite(rng, d, m, bound, family)
        mul = sp.Rational(rng.choice((-3, 2, 5)), 7)
        scale, shift = sp.Rational(2, rng.choice((3, 5, 7))), rng.choice((-3, 1, 2))
        transformed = affine(source, scale, shift).mul_ground(mul) if index % 2 else source
        base = {**detail, "affine": [str(scale), str(shift), str(mul)] if index % 2 else None}
        add(transformed, "constant", "algebraic-cover", base)
        square = poly((X - rng.randint(-8, 8))**(1 + index % 3) * (X**2 + 1))
        if index % 2 == 0:
            add(transformed * square**2, "constant", "inflated-cover",
                {**base, "inflator": coefficients(square)})
        for degree in range(1, transformed.degree()):
            altered = transformed + poly(X**degree)
            if altered.gcd(altered.diff()).degree() == 0:
                add(altered, "constant", "one-coefficient-perturbation",
                    {**base, "perturbation": [degree, "1"]})
                break

    target = poly(X**3 - X + 3)
    for bound in (2, 4):
        sixth = poly(2 * sp.chebyshevt(6, X / 2)).shift(rng.randint(-bound, bound))
        add(target.compose(sixth), "constant", "multiple-degrees",
            {"construction": "commuting Chebyshev steps", "inner": coefficients(sixth),
             "outer": coefficients(target)})

    for length, repeats in ((3, 4), (4, 1), (4, 5), (6, 3), (3, 5), (1, 6)):
        stiffness = rng.sample(range(1, 21), length)
        uniform = length == 1
        if uniform:
            stiffness *= 3
        center = 2 * max(stiffness) + 3
        diagonal = ([center] * len(stiffness) if length % 2 == 0 or uniform else
                    [stiffness[i - 1] + stiffness[i] + rng.randint(0, 4)
                     for i in range(len(stiffness))])
        bonds = [-v for v in stiffness]
        discriminant, product = cell_discriminant(diagonal, bonds)
        trace, _ = trace_and_quotient(discriminant, product, repeats)
        add(poly(trace**2 - 4 * product**(2 * repeats)), "constant", "jacobi-spring",
            {"diagonal": diagonal, "bonds": bonds, "repetitions": repeats,
             "pinning": [a + bonds[i - 1] + bonds[i] for i, a in enumerate(diagonal)],
             "declared_nonminimal": uniform})

    add(poly(-sp.Rational(7, 5) * (X + 1)**4 * (X**2 + 1)**2),
        "constant", "genus-zero-control", {"construction": "constant times square"})
    add(poly(X**7 + X + 1), "constant", "prime-degree-control", {"construction": "degree seven"})

    regression_u, regression_g = poly(X**3 + 2 * X**2), poly(X**3 + X + 1)
    add(regression_g.compose(regression_u), "constant", "stock-regression-control",
        {"construction": "nonzero consecutive right-component coefficients",
         "inner": coefficients(regression_u), "outer": coefficients(regression_g)})

    for index, (u, g, v, family) in enumerate(variable_fixtures(rng)):
        f = g.compose(u).exquo(v**2)
        if f.gcd(f.diff()).degree() != 0:
            raise ArithmeticError("A ramified-cover fixture is singular.")
        shift = rng.randint(-3, 3)
        transformed = f.shift(shift)
        details = {"construction": family, "inner": coefficients(u),
                   "outer": coefficients(g), "vertical": coefficients(v), "shift": shift,
                   "degree_stratum": [u.degree(), g.degree(), v.degree()]}
        add(transformed, "nonconstant", family, details)
        altered = transformed + (index + 1)
        if altered.gcd(altered.diff()).degree() != 0:
            raise ArithmeticError("A perturbed fixture is singular.")
        add(altered, "nonconstant", "one-coefficient-perturbation", {**details, "perturbation": [0, str(index + 1)]})

    rng.shuffle(cases)
    public, sidecar = {}, {}
    for i, (input_case, provenance) in enumerate(cases, 1):
        identifier = f"q{i:03}"
        public[identifier] = input_case
        sidecar[identifier] = provenance
    return public, {"seed": seed, "cohorts": dict(sorted(Counter(v["cohort"] for v in sidecar.values()).items())),
                    "cases": sidecar}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("seed", type=int)
    parser.add_argument("destination", type=Path)
    args = parser.parse_args()
    public, sidecar = generate(args.seed)
    args.destination.mkdir(parents=True, exist_ok=True)
    for name, data in (("inputs", public), ("provenance", sidecar)):
        (args.destination / f"{name}.json").write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


if __name__ == "__main__":
    main()
