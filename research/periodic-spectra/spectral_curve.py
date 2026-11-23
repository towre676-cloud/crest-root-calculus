"""Verify the classical normalization identity for repeated Jacobi cells.

This is an exact companion certificate, not a timed spectral solver. The
general proof is in spectral-curve.md; finite symbolic checks exercise it.
"""
import json

import sympy as sp


X = sp.Symbol("x")


def trace_and_quotient(discriminant, product, repetitions):
    """Return R_N(D,B) and S_(N-1)(D,B), with no division by B."""
    if type(repetitions) is not int or repetitions < 1:
        raise ValueError("The cell repetition count must be a positive integer.")
    previous, current = sp.Integer(2), sp.sympify(discriminant)
    q_previous, quotient = sp.Integer(0), sp.Integer(1)
    for _ in range(2, repetitions + 1):
        previous, current = current, sp.expand(
            discriminant * current - product**2 * previous)
        q_previous, quotient = quotient, sp.expand(
            discriminant * quotient - product**2 * q_previous)
    return current, quotient


def cell_discriminant(diagonal, bonds):
    """Independent symbolic multiplication of the two-by-two cell matrices."""
    if (not diagonal or len(diagonal) != len(bonds)
            or any(type(v) is not int for v in diagonal + bonds)
            or any(v == 0 for v in bonds)):
        raise ValueError("Use equal nonempty integer arrays and nonzero bonds.")
    matrix = sp.eye(2)
    for i, a in enumerate(diagonal):
        matrix = sp.Matrix([[X - a, -bonds[i - 1]], [bonds[i], 0]]) * matrix
    return sp.expand(sp.trace(matrix)), sp.prod(bonds)


def remove_square_factors(expression):
    """Return monic H,F with expression=H^2 F and F squarefree.

    Even-multiplicity factors disappear from F. The radical retaining every
    distinct factor would generally give the WRONG normalization curve.
    """
    polynomial = sp.Poly(expression, X, domain=sp.QQ)
    if polynomial.degree() < 1 or polynomial.LC() != 1:
        raise ValueError("This certificate handles nonconstant monic radicands.")
    content, factors = polynomial.sqf_list()
    square = sp.Poly(1, X, domain=sp.QQ)
    reduced = sp.Poly(content, X, domain=sp.QQ)
    for factor, multiplicity in factors:
        square *= factor ** (multiplicity // 2)
        if multiplicity % 2:
            reduced *= factor
    if square**2 * reduced != polynomial or reduced.gcd(reduced.diff()).degree() != 0:
        raise ArithmeticError("Square removal failed its exact identity.")
    return square, reduced


def ascending(polynomial):
    polynomial = sp.Poly(polynomial, X, domain=sp.QQ)
    return [str(polynomial.nth(i)) for i in range(polynomial.degree() + 1)]


def verify_example(name, diagonal, bonds, repetitions, expected_genus):
    from model import coefficient_period, discriminant

    d, product = cell_discriminant(diagonal, bonds)
    coefficients, implementation_product = discriminant(diagonal, bonds)
    implementation_d = sum(c * X**i for i, c in enumerate(coefficients))
    if d != implementation_d or product != implementation_product:
        raise ArithmeticError("Independent cell matrix disagrees with the solver.")
    trace, quotient = trace_and_quotient(d, product, repetitions)
    base = sp.Poly(d**2 - 4 * product**2, X, domain=sp.QQ)
    repeated = sp.Poly(trace**2 - 4 * product**(2 * repetitions), X, domain=sp.QQ)
    q = sp.Poly(quotient, X, domain=sp.QQ)
    if repeated != q**2 * base:
        raise ArithmeticError("Repeated-cell radicand identity failed.")
    h, reduced = remove_square_factors(base)
    repeated_h, repeated_reduced = remove_square_factors(repeated)
    if repeated_reduced != reduced or repeated_h != q * h:
        raise ArithmeticError("Repeated and original curves have different normalizations.")
    degree = int(reduced.degree())
    if degree < 2 or degree % 2:
        raise ArithmeticError("Unexpected normalized spectral-curve degree.")
    genus = (degree - 2) // 2
    if genus != expected_genus:
        raise ArithmeticError("An example's exact normalization genus changed.")
    real_branch_points = int(reduced.count_roots(-sp.oo, sp.oo))
    if real_branch_points != degree:
        raise ArithmeticError("A real Jacobi spectral branch point became nonreal.")
    return {
        "example": name, "diagonal": diagonal, "bonds": bonds,
        "declared_cell_length": len(diagonal),
        "least_coefficient_period": coefficient_period(diagonal, bonds),
        "repetitions": repetitions, "bond_product": int(product),
        "discriminant_coefficients": ascending(d),
        "base_radicand_coefficients": ascending(base),
        "normalized_radicand_coefficients": ascending(reduced),
        "base_square_multiplier_coefficients": ascending(h),
        "repetition_multiplier_degree": int(q.degree()),
        "repeated_radicand_degree": int(repeated.degree()),
        "normalized_radicand_degree": degree, "normalization_genus": genus,
        "real_branch_points": real_branch_points,
        "declared_cell_closed_gaps": len(diagonal) - 1 - genus,
    }


def verify_all():
    d, b = sp.symbols("d b")
    identity_checks = 0
    for n in range(1, 13):
        trace, quotient = trace_and_quotient(d, b, n)
        residual = sp.expand(trace**2 - 4*b**(2*n) - quotient**2*(d**2 - 4*b**2))
        if residual != 0:
            raise ArithmeticError("A universal polynomial recurrence identity failed.")
        identity_checks += 1
    examples = [
        ("spring-two-cell", [5, 6], [-1, -3], 4, 1),
        ("spring-three-cell", [3, 7, 7], [-1, -4, -2], 5, 2),
        ("spring-four-cell", [6, 7, 9, 6], [-2, -5, -1, -3], 3, 3),
        ("signed-three-cell", [3, -2, 5], [1, -2, 3], 3, 2),
        ("nonminimal-uniform-three-cell", [2, 2, 2], [-1, -1, -1], 4, 0),
    ]
    reports = [verify_example(*example) for example in examples]
    return {"universal_identity_checks": identity_checks,
            "examples": reports, "checks_passed": True}


if __name__ == "__main__":
    print(json.dumps(verify_all(), indent=2))
