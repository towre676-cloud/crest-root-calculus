"""Compute the non-CM elliptic spectrum of a genus-two cyclic cubic equation."""

import argparse
import json
from pathlib import Path
import sys

import sympy as sp

from reconstruct import invariant, leading_constant, non_cm_witness

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "cyclic-cover-reconstruction"))
from cyclic import coefficients, polynomial, power_reduce, source_model


def hyperelliptic_chart(values):
    """Write xi^3 = A/C and take the quadratic discriminant in x."""
    removed, branch, genus, _ = source_model(values, 3)
    if genus != 2:
        raise ValueError("The normalized cyclic cubic equation must have genus two.")

    # The four branch points have weights 1,1,2,2. Make infinity a
    # weight-one point, leaving deg(A)=1 and deg(C)=2 at finite places.
    character = 1 if branch.degree() == 5 else 2
    additional, normalized = power_reduce(branch, 3, character)
    unit, components = normalized.sqf_list()
    factors = {order: poly for poly, order in components}
    A, C = factors[1].mul_ground(unit), factors[2]
    if A.degree() != 1 or C.degree() != 2:
        raise ArithmeticError("The normalized genus-two branch profile is inconsistent.")
    denominator = removed**character*additional*C
    if polynomial(values)**character*C != denominator**3*A:
        raise ArithmeticError("The exact cubic change of coordinates failed.")

    a0, a1 = (A.nth(i) for i in range(2))
    c0, c1, c2 = (C.nth(i) for i in range(3))
    P, Q, R = c1*c1-4*c2*c0, 4*c2*a0-2*c1*a1, a1*a1
    if not P or not R or Q*Q == 4*P*R:
        raise ArithmeticError("The hyperelliptic sextic is singular.")
    return {
        "source_y_power": character,
        "xi_denominator": coefficients(denominator),
        "A": coefficients(A), "C": coefficients(C),
        "eta_constant": str(-a1),
        "eta_xi_cubed_coefficient": [str(c1), str(2*c2)],
        "sextic_coefficients": [str(R), "0", "0", str(Q), "0", "0", str(P)],
    }


def quotient_j(P, Q, R):
    t = sp.Symbol("t")
    scale = sp.sqrt(P*R)
    lam = Q/scale
    field = sp.QQ if scale.is_Rational else sp.QQ.algebraic_field(scale)
    j = []
    for parameter in (lam, -lam):
        cubic = sp.Poly(t**3+(-6*parameter-15)*t+parameter**2+14*parameter+22,
                        t, domain=field)
        j.append(invariant(cubic))
    return j, field


def degree_counts(bound):
    """Sieve the Eisenstein norm multiplicities through a degree bound."""
    limit = bound//2
    counts = [6]*(limit+1)
    for p in sp.primerange(2, limit+1):
        factor = 2 if p % 3 == 1 else 0
        for m in range(p, limit+1, p):
            counts[m] *= factor
    return {2*m: counts[m] for m in range(1, limit+1) if counts[m]}


def compute(values, bound=150):
    if type(bound) is not int or bound < 0:
        raise ValueError("Use a nonnegative integer spectrum bound.")
    chart = hyperelliptic_chart(values)
    P, Q, R = (sp.Rational(chart["sextic_coefficients"][i]) for i in (6, 3, 0))
    j, field = quotient_j(P, Q, R)
    result = {"chart": chart, "quotient_j": list(map(str, j)), "spectrum_bound": bound}
    witness = non_cm_witness(j[0], field) or non_cm_witness(j[1], field)
    if witness is None:
        result["status"] = "non-CM not certified"
        return result
    # Both reciprocal lifts descend exactly when R/P is a rational sixth power.
    ratio = R/P
    sixth_power = (ratio > 0 and sp.integer_nthroot(int(ratio.p), 6)[1]
                   and sp.integer_nthroot(int(ratio.q), 6)[1])
    result.update(status="certified", non_cm_witness=witness, minimal_isogeny_degree=3,
                  fixed_nonzero_two_torsion=3, equivariant_elliptic_maps=False,
                  geometric_optimal_counts=degree_counts(bound),
                  rational_optimal_counts={2: 2} if sixth_power and bound >= 2 else {},
                  rational_bielliptic_lifts=bool(sixth_power),
                  leading_numerator=str(leading_constant(3, 3)))
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON with exponent 3 and ascending coefficients")
    parser.add_argument("--bound", type=int, default=150)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    if data["exponent"] != 3:
        parser.error("the exponent must be 3")
    print(json.dumps(compute(data["coefficients"], args.bound), indent=2))


if __name__ == "__main__":
    main()
