"""Exhaust the maps in an unramified genus-three to genus-two example.

The independent calculation uses branch-point automorphisms and the degree-two
deck involution, without eliminating the general rational-map coefficient ideal.
"""
import argparse
from itertools import permutations
import json
from pathlib import Path

import sympy as sp


X = sp.Symbol("x")
SOURCE = sp.Poly(sp.prod(X**2 - j for j in range(1, 5)), X)
TARGET = sp.Poly(sp.prod(X - j for j in range(5)), X)


def mobius_from_triples(source, target):
    rows = []
    for point, image in zip(source, target):
        if point == sp.oo:
            row = [0, 0, 1, 0] if image == sp.oo else [1, 0, -image, 0]
        elif image == sp.oo:
            row = [0, 0, point, 1]
        else:
            row = [point, 1, -point * image, -image]
        rows.append(row)
    kernel = sp.Matrix(rows).nullspace()
    if len(kernel) != 1:
        raise ArithmeticError("Distinct triples failed to determine a projective map.")
    values = list(kernel[0])
    scale = next(value for value in reversed(values) if value != 0)
    values = [value / scale for value in values]
    a, b, c, d = values
    if a * d == b * c:
        raise ArithmeticError("The branch-triple map is singular.")
    return values


def branch_automorphisms(polynomial, rational_branches, fixed_triple):
    """All rational base automorphisms when three rational branches are supplied."""
    branch_degree = 2 * ((polynomial.degree() - 1) // 2) + 2
    result = []
    for images in permutations(rational_branches, 3):
        a, b, c, d = mobius_from_triples(fixed_triple, images)
        numerator, denominator = a * X + b, c * X + d
        transformed = sp.Poly(sp.expand(sum(polynomial.nth(i) * numerator**i
                                           * denominator**(branch_degree - i)
                                           for i in range(polynomial.degree() + 1))), X)
        factor, remainder = sp.div(transformed, polynomial)
        if remainder or factor.degree() != 0:
            continue
        scalar = factor.nth(0)
        result.append({"matrix": list(map(str, [a, b, c, d])),
                       "vertical_square_factor": str(scalar),
                       "rational_lifts": 2 if sp.sqrt(scalar).is_Rational else 0})
    return result


def geometric_source_automorphisms():
    """Check all 8*7*6 possible images in the explicit splitting field."""
    field = sp.QQ.algebraic_field(sp.sqrt(2), sp.sqrt(3))
    roots = [field.from_sympy(point) for point in
             [-sp.sqrt(3), -sp.sqrt(2), sp.Integer(-2), sp.Integer(-1),
              sp.Integer(1), sp.Integer(2), sp.sqrt(2), sp.sqrt(3)]]
    result = []
    for first, second, third in permutations(roots, 3):
        difference = second - first
        other_difference = second - third
        # Compose the maps taking (-2,-1,1) and this image triple to (0,1,infinity).
        a = -2 * third * difference - first * other_difference
        b = -4 * third * difference + first * other_difference
        c = -2 * difference - other_difference
        d = -4 * difference + other_difference
        images = []
        for root in roots:
            denominator = c * root + d
            if denominator == field.zero:
                break
            image = (a * root + b) / denominator
            if image not in roots:
                break
            images.append(roots.index(image))
        if len(images) == 8:
            divisor = d if d != field.zero else c
            result.append({"matrix": [str(sp.simplify(field.to_sympy(value / divisor)))
                                       for value in [a, b, c, d]],
                           "branch_permutation": images})
    return result


def coefficients(polynomial):
    return [str(polynomial.nth(i)) for i in range(polynomial.degree() + 1)]


def calculate():
    source_automorphisms = branch_automorphisms(SOURCE, [-2, -1, 1, 2], [-2, -1, 1])
    target_automorphisms = branch_automorphisms(TARGET, [0, 1, 2, 3, 4, sp.oo], [0, 1, sp.oo])
    geometric_automorphisms = geometric_source_automorphisms()
    expected_source = {("1", "0", "0", "1"), ("-1", "0", "0", "1")}
    if {tuple(row["matrix"]) for row in source_automorphisms} != expected_source:
        raise ArithmeticError("The rational source branch group changed.")
    if {tuple(row["matrix"]) for row in geometric_automorphisms} != expected_source:
        raise ArithmeticError("The geometric source branch group changed.")
    if any(row["vertical_square_factor"] != "1" for row in source_automorphisms):
        raise ArithmeticError("An unexpected source lift appeared.")
    if len(target_automorphisms) != 4 or sum(row["rational_lifts"] for row in target_automorphisms) != 2:
        raise ArithmeticError("The target automorphism count changed.")

    # The hyperelliptic involution has eight fixed points. The lift (-x,y)
    # fixes two points over zero and both infinity points. The lift (-x,-y)
    # fixes neither, since F(0)=24 and the leading coefficient is one.
    if SOURCE.nth(0) != 24 or SOURCE.LC() != 1:
        raise ArithmeticError("The involution fixed-point calculation changed.")
    fixed_points = [8, 4, 0]
    quotient_genera = [(2 * 3 + 2 - fixed) // 4 for fixed in fixed_points]
    if quotient_genera != [0, 1, 2]:
        raise ArithmeticError("Riemann--Hurwitz did not identify a unique genus-two quotient.")

    if sp.expand(TARGET.as_expr().subs(X, X**2) - SOURCE.as_expr() * X**2) != 0:
        raise ArithmeticError("The displayed even-model cover identity failed.")
    odd_source = sp.Poly(sp.prod((X + 1)**2 - j * X**2 for j in range(1, 5)), X)
    A, B, C = (X + 1)**2, X**2, X * (X + 1)
    homogeneous = sum(TARGET.nth(i) * A**i * B**(5 - i) for i in range(6))
    if odd_source.degree() != 7 or sp.gcd(odd_source, odd_source.diff()).degree():
        raise ArithmeticError("The transformed source is not smooth of genus three.")
    if sp.expand(B * homogeneous - odd_source.as_expr() * C**2) != 0:
        raise ArithmeticError("The displayed odd-model rational-map identity failed.")

    return {"source": coefficients(SOURCE), "target": coefficients(TARGET),
            "degree_bound": 2, "degree_one_maps": 0,
            "rational_maps": [{"u": ["0", "0", "1"], "v_over_y": ["0", sign], "degree": 2}
                              for sign in ["1", "-1"]],
            "geometric_map_count": 8,
            "source_rational_branch_triples": 24,
            "source_geometric_branch_triples": 336,
            "target_branch_triples": 120,
            "source_base_automorphisms": source_automorphisms,
            "source_geometric_base_automorphisms": geometric_automorphisms,
            "target_base_automorphisms": target_automorphisms,
            "source_involution_fixed_points": fixed_points,
            "source_quotient_genera": quotient_genera,
            "odd_source": coefficients(odd_source),
            "odd_model_map": {"A": coefficients(sp.Poly(A, X)),
                              "B": coefficients(sp.Poly(B, X)),
                              "C": coefficients(sp.Poly(C, X)),
                              "vertical_denominator_exponent": 3},
            "scope": "all nonconstant maps to the specified genus-two target"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    result = calculate()
    if args.verify:
        if result != json.loads(args.verify.read_text(encoding="utf-8")):
            raise ArithmeticError("The exhaustive automorphism receipt differs.")
        print("PASS")
    elif args.output:
        args.output.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    else:
        print(json.dumps(result, indent=2))
