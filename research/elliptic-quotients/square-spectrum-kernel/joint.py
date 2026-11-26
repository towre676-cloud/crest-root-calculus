"""Certify the square and full-theta kernels across both genera at each level."""
import argparse
from fractions import Fraction
import json
from math import gcd, lcm
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory

from sympy import Matrix, factorint

from replay import HERE, center, coefficient_rows, image, pivot_indices, residue_field


def primitive(column):
    values = [Fraction(x) for x in column]
    denominator = lcm(*(x.denominator for x in values))
    integral = [int(x * denominator) for x in values]
    common = gcd(*integral)
    return [x // common for x in integral]


def integer_kernel(rows):
    """Return primitive integer vectors forming a basis over the rationals."""
    independent = [rows[i] for i in pivot_indices(rows)]
    return [primitive(column) for column in Matrix(independent).nullspace()]


def certify_level(binary, level):
    D, ell = level["D"], level["ell"]
    genera = level["genera"]
    forms = [form for genus in genera for form in genus["forms"]]
    width = len(forms)
    modulus = -4 * D * ell * ell
    if len(genera) != 2 or any(g["level"] != modulus for g in genera):
        raise ArithmeticError("Unexpected common modular level.")
    index = modulus
    for p in factorint(modulus):
        index = index // p * (p + 1)
    square_bound = index // 6
    theta_bound = index // 8
    squares = [m*m for m in range(1, square_bound + 1)]
    indices = sorted(set(squares) | set(range(1, theta_bound + 1)))
    coefficients = coefficient_rows(binary, D, ell, forms, indices)

    cut = len(genera[0]["forms"])
    genus_indicators = [[int(i < cut) for i in range(width)], [int(i >= cut) for i in range(width)]]
    short = genus_indicators + [coefficients[n] for n in squares[:300]]
    kernel = integer_kernel(short)
    if kernel:
        restrictions = [image(coefficients[n], kernel) for n in squares[300:]]
        if any(map(any, restrictions)):
            kernel = [primitive(image(c, list(zip(*kernel)))) for c in integer_kernel(restrictions)]
    if any(any(image(coefficients[n], kernel)) for n in squares):
        raise ArithmeticError("Square-kernel certificate failed.")
    theta_rows = [image(coefficients[n], kernel) for n in range(1, theta_bound + 1)]
    theta_indices = [i + 1 for i in pivot_indices(theta_rows)] if kernel else []
    theta_rank = len(theta_indices)
    # Relations supported in either genus are the kernels of the two projections.
    coordinate_rows = list(zip(*kernel))
    within = 2 * len(kernel) - sum(
        len(list(pivot_indices(part))) for part in (coordinate_rows[:cut], coordinate_rows[cut:])
    )
    result = {
        "D": D, "ell": ell, "forms": width, "sturm": square_bound, "theta_cutoff": theta_bound,
        "joint_square_nullity": len(kernel), "within_genus_nullity": within,
        "joint_full_theta_nullity": len(kernel) - theta_rank,
        "first_full_theta_indices": theta_indices,
    }
    if (D, ell) == (-3, 5):
        if kernel != [[1, -1, -2, 2]] or theta_indices != [12]:
            raise ArithmeticError("Four-form example changed.")
        field = residue_field(D, ell)
        result["example"] = {
            "centers_in_order": [center(field, D, f["psi"]) for f in forms],
            "relation": kernel[0],
            "first_full_theta_index": 12,
            "ordinary_counts_at_first_difference": coefficients[12],
            "difference": image(coefficients[12], kernel)[0],
        }
    print((D, ell), len(kernel), len(kernel) - theta_rank, file=sys.stderr, flush=True)
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    with TemporaryDirectory() as td:
        binary = Path(td) / "count"
        run([args.compiler, "-std=c++17", "-O3", str(HERE / "count.cpp"), "-o", str(binary)], check=True)
        source = json.loads((HERE.parent / "inverse-spectrum" / "imported-census.json").read_text())
        results = [certify_level(binary, level) for level in source["levels"]]
        assert sum(r["forms"] for r in results) == 463
        assert sum(r["joint_square_nullity"] for r in results) == 202
        assert sum(r["within_genus_nullity"] for r in results) == 6
        assert all(r["joint_full_theta_nullity"] == 0 for r in results)
        print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
