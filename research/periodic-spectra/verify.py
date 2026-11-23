"""Independently verify full cyclic spectra and adversarial structural controls.

The verifier expands the entire input and uses open-path determinants. It
does not use the period detector, cell transfer matrices, or Chebyshev
recurrence of the measured solver. Coefficients are ascending.
"""

import argparse
import copy
import json
import math
from pathlib import Path
import random
import re


HERE = Path(__file__).resolve().parent


def _integer(value, label):
    if type(value) is not int:
        raise ValueError(f"{label} must be an integer, not a float or boolean")
    return value


def expanded_input(case):
    """Independent input expansion; reject ambiguous matrix conventions."""
    diagonal = list(case["diagonal"])
    bonds = list(case["bonds"])
    repeats = _integer(case["repeats"], "repeats")
    if not diagonal or len(diagonal) != len(bonds) or repeats < 1:
        raise ValueError("nonempty equal-length cells and positive repeats required")
    diagonal = [_integer(a, "diagonal entry") for a in diagonal]
    bonds = [_integer(b, "bond") for b in bonds]
    if any(b == 0 for b in bonds):
        raise ValueError("every cycle bond must be nonzero")
    a, b = diagonal * repeats, bonds * repeats
    if len(a) < 3:
        raise ValueError("the expanded cycle must have at least three vertices")
    if "defect" in case:
        defect = case["defect"]
        index = _integer(defect["index"], "defect index")
        delta = _integer(defect["diagonal_delta"], "diagonal defect")
        if index < 0 or index >= len(a):
            raise ValueError("defect lies outside the cycle")
        a[index] += delta
    return a, b


def independent_matrix(case):
    a, b = expanded_input(case)
    size = len(a)
    matrix = [[0] * size for _ in range(size)]
    for i in range(size):
        matrix[i][i] = a[i]
        j = (i + 1) % size
        matrix[i][j] += b[i]
        matrix[j][i] += b[i]
    return matrix


def _open_path(a, b):
    """det(x I - A) for an open path, by expansion at its last vertex."""
    previous_previous = [0]
    previous = [1]
    for i, diagonal in enumerate(a):
        current = [0] * (len(previous) + 1)
        for degree, value in enumerate(previous):
            current[degree] -= diagonal * value
            current[degree + 1] += value
        if i:
            edge_square = b[i - 1] ** 2
            for degree, value in enumerate(previous_previous):
                current[degree] -= edge_square * value
        previous_previous, previous = previous, current
    return previous


def independent_coefficients(case):
    """Full exact cyclic determinant; deliberately ignore all periodicity."""
    a, b = expanded_input(case)
    coefficients = _open_path(a, b[:-1])
    interior = _open_path(a[1:-1], b[1:-2])
    for degree, value in enumerate(interior):
        coefficients[degree] -= b[-1] ** 2 * value
    coefficients[0] -= 2 * math.prod(b)
    if len(coefficients) != len(a) + 1 or coefficients[-1] != 1:
        raise ArithmeticError("independent determinant has the wrong degree")
    return coefficients


def _output_integer(value):
    if type(value) is int:
        return value
    if isinstance(value, str) and re.fullmatch(r"-?(0|[1-9][0-9]*)", value):
        return int(value)
    raise ArithmeticError("an output coefficient is not an exact integer")


def verify(case, actual):
    """Check every coefficient, not just samples, traces or a hash."""
    expected = independent_coefficients(case)
    if not isinstance(actual, (list, tuple)):
        raise ArithmeticError("the output must be a complete coefficient sequence")
    actual = [_output_integer(value) for value in actual]
    if actual != expected:
        raise ArithmeticError(f"incorrect characteristic polynomial: {case.get('id', '?')}")
    return {
        "id": case.get("id"),
        "degree": len(expected) - 1,
        "all_coefficients_exact": True,
        "verification": "full-ring open-path determinants",
    }


def _period(a, b):
    length = len(a)
    return next(p for p in range(1, length + 1)
                if length % p == 0 and all(a[i] == a[i % p]
                                           and b[i] == b[i % p]
                                           for i in range(length)))


def _small_cases():
    cases = [
        {"id": "triangle-uniform", "diagonal": [0], "bonds": [1], "repeats": 3},
        {"id": "triangle-negative", "diagonal": [0], "bonds": [-1], "repeats": 3},
        {"id": "triangle-mixed", "diagonal": [2, -3, 1], "bonds": [2, -1, 3], "repeats": 1},
        {"id": "cell-two", "diagonal": [2, -1], "bonds": [-2, 3], "repeats": 2},
        {"id": "cell-two-odd", "diagonal": [2, -1], "bonds": [-2, 3], "repeats": 3},
        {"id": "hidden-period", "diagonal": [2, -1, 2, -1],
         "bonds": [-2, 3, -2, 3], "repeats": 2},
        {"id": "nonperiodic", "diagonal": [2, 1, -4, 3, 0],
         "bonds": [1, -2, 3, -4, 5], "repeats": 1},
        {"id": "single-defect", "diagonal": [2, -1], "bonds": [-2, 3],
         "repeats": 4, "defect": {"index": 3, "diagonal_delta": 1}},
        {"id": "zero-defect", "diagonal": [2, -1], "bonds": [-2, 3],
         "repeats": 3, "defect": {"index": 3, "diagonal_delta": 0}},
    ]
    rng = random.Random(736210)
    for index in range(24):
        size = rng.randrange(1, 5)
        repeats = rng.randrange(1, 5)
        if size * repeats < 3:
            repeats = 3
        cases.append({
            "id": f"seeded-{index}",
            "diagonal": [rng.randrange(-4, 5) for _ in range(size)],
            "bonds": [rng.choice([-3, -2, -1, 1, 2, 3]) for _ in range(size)],
            "repeats": repeats,
        })
    return cases


def _spring_cases():
    for repeats in (1, 2, 3):
        for pinned in (False, True):
            springs = [2, 3, 5]
            pinning = [0, 1, 0] if pinned else [0, 0, 0]
            yield {
                "id": f"springs-{repeats}-{'pinned' if pinned else 'free'}",
                "diagonal": [springs[i - 1] + springs[i] + pinning[i]
                             for i in range(3)],
                "bonds": [-k for k in springs],
                "repeats": repeats,
            }, pinning * repeats


def verify_controls():
    """Exercise actual solver entry points and independent matrix identities."""
    import sympy as sp
    import model

    x = sp.Symbol("x")
    counters = {"exact_model_cases": 0, "direct_symbolic_determinants": 0,
                "symbolic_characteristic_polynomials": 0,
                "altered_outputs_rejected": 0, "malformed_matrices_rejected": 0,
                "spring_energy_controls": 0}
    for case in _small_cases():
        expected_matrix = independent_matrix(case)
        if model.matrix_from_case(case) != expected_matrix:
            raise ArithmeticError("matrix constructor changed the supplied cycle")
        a, b = expanded_input(case)
        p = _period(a, b)
        actual_a, actual_b, actual_n = model.detect_period(expected_matrix)
        if (list(actual_a), list(actual_b), actual_n) != (a[:p], b[:p], len(a) // p):
            raise ArithmeticError("period recognition did not return the true least period")
        coefficients = model.solve(case)
        verify(case, coefficients)
        counters["exact_model_cases"] += 1
        matrix = sp.Matrix(expected_matrix)
        if len(a) <= 6:
            polynomial = sp.Poly((x * sp.eye(len(a)) - matrix).det(method="domain-ge"), x)
            counters["direct_symbolic_determinants"] += 1
        else:
            polynomial = matrix.charpoly(x).as_poly()
            counters["symbolic_characteristic_polynomials"] += 1
        if [int(polynomial.nth(i)) for i in range(len(a) + 1)] != coefficients:
            raise ArithmeticError("independent symbolic matrix determinant disagrees")
        altered = list(coefficients)
        altered[len(a) // 2] += 1
        try:
            verify(case, altered)
        except ArithmeticError:
            counters["altered_outputs_rejected"] += 1
        else:
            raise ArithmeticError("an altered integer coefficient was accepted")

    for case, pinning in _spring_cases():
        coefficients = model.solve(case)
        verify(case, coefficients)
        matrix = independent_matrix(case)
        length = len(matrix)
        if [sum(row) for row in matrix] != pinning:
            raise ArithmeticError("spring energy row sums differ from pinning")
        if not any(pinning):
            if coefficients[0] != 0 or coefficients[1] == 0:
                raise ArithmeticError("an unpinned connected ring lost its simple zero mode")
        elif (-1) ** length * coefficients[0] <= 0:
            raise ArithmeticError("a pinned spring ring has a nonpositive determinant")
        vector = [i % 4 - 1 for i in range(length)]
        a, b = expanded_input(case)
        quadratic = sum(vector[i] * matrix[i][j] * vector[j]
                        for i in range(length) for j in range(length))
        energy = sum(-b[i] * (vector[(i + 1) % length] - vector[i]) ** 2
                     + pinning[i] * vector[i] ** 2 for i in range(length))
        if quadratic != energy:
            raise ArithmeticError("the matrix does not represent the claimed spring energy")
        counters["spring_energy_controls"] += 1

    base = independent_matrix({"id": "invalid-base", "diagonal": [2, -1],
                               "bonds": [-2, 3], "repeats": 3})
    invalid = [[], [[1]], [[1, 2], [2, 1]], [row[:-1] for row in base]]
    for changes in [((0, 2, 1), (2, 0, 1)), ((0, 1, 7),),
                    ((0, 1, 0), (1, 0, 0)), ((0, 0, 2.0),),
                    ((0, 0, True),)]:
        altered = copy.deepcopy(base)
        for i, j, value in changes:
            altered[i][j] = value
        invalid.append(altered)
    for matrix in invalid:
        try:
            model.detect_period(matrix)
        except (model.StructureError, TypeError, ValueError):
            counters["malformed_matrices_rejected"] += 1
        else:
            raise ArithmeticError("an off-contract matrix was treated as a Jacobi ring")

    return counters


def main():
    from spectral_curve import verify_all as verify_spectral_curves

    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel", type=Path,
                        help="also recompute and independently certify every case in a JSON panel")
    args = parser.parse_args()
    result = verify_controls()
    result["spectral_curve"] = verify_spectral_curves()
    if args.panel:
        import model
        panel = json.loads(args.panel.read_text())
        cases = panel["cases"] if isinstance(panel, dict) else panel
        for case in cases:
            verify(case, model.solve(case))
        result["panel_cases_verified"] = len(cases)
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
