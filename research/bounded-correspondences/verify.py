"""Replay exact map certificates and the finite-local reconstruction controls."""
import argparse
from copy import deepcopy
import json
from pathlib import Path

import sympy as sp

import composition
import finite
import maps
import rational_maps
import frobenius_filter
import unramified_example


HERE = Path(__file__).resolve().parent


def check_record(row, rebuild=False):
    result = row["result"]
    if row["kind"] == "constant":
        composition.verify_constant(result)
        fresh = composition.recover_constant(result["source"], result["target"]) if rebuild else result
        counts = [fresh["geometric_maps"], len(fresh["rational_maps"])]
    elif row["kind"] == "polynomial":
        maps.verify_maps(result)
        fresh = maps.recover(result["source"], result["target"]) if rebuild else result
        counts = [sum(r["geometric_maps"] for r in fresh["degrees"]),
                  sum(len(r["rational_maps"]) for r in fresh["degrees"])]
    elif row["kind"] == "rational":
        rational_maps.verify_maps(result)
        fresh = rational_maps.recover(result["source"], result["target"]) if rebuild else result
        counts = [sum(r["geometric_maps"] for r in fresh["charts"]),
                  sum(len(r["rational_maps"]) for r in fresh["charts"])]
    else:
        finite.verify(result)
        if rebuild:
            variables = sp.symbols(result["variables"])
            context = dict(zip(result["variables"], variables))
            equations = [sp.sympify(e, locals=context) for e in result["equations"]]
            fresh = finite.rational_points(equations, variables)
        else:
            fresh = result
        counts = [fresh["geometric_points"], len(fresh["rational_points"])]
    if counts != row["counts"]:
        raise ArithmeticError(row["id"] + " has different geometric or rational counts.")
    return counts


def must_reject(operation):
    try:
        operation()
    except (ArithmeticError, ValueError):
        return
    raise ArithmeticError("An invalid input or altered certificate was accepted.")


def controls():
    x, y = sp.symbols("x y")
    must_reject(lambda: finite.rational_points([x * y], [x, y]))
    must_reject(lambda: maps.curve([0, 0, 0, 1]))
    rational = finite.rational_points([2 * x - 3, 3 * y + 5], [x, y])
    changed = deepcopy(rational)
    changed["algebra_dimension"] += 1
    must_reject(lambda: finite.verify(changed))
    changed = deepcopy(rational)
    changed["coordinates"].pop()
    must_reject(lambda: finite.verify(changed))
    changed = deepcopy(rational)
    changed["lifting"]["branches"][0]["lifted_root"] = "0"
    must_reject(lambda: finite.verify(changed))
    changed = deepcopy(rational)
    changed["multiplication_matrices"][0][0][0] = "7"
    must_reject(lambda: finite.verify(changed))
    must_reject(lambda: finite.rational_points([2 * x - 3], [x], prime=2))
    sample = maps.recover([4, 12, -1, -7, 0, 1], [4, -4, -1, 1])
    changed = deepcopy(sample)
    changed["degrees"][0]["geometric_maps"] += 1
    must_reject(lambda: maps.verify_maps(changed))
    changed = deepcopy(sample)
    changed["degree_bound"] += 2
    must_reject(lambda: maps.verify_maps(changed))
    changed = deepcopy(rational)
    changed["minimal_polynomial"] = [str(3 * sp.Rational(value) / 2)
                                      for value in changed["minimal_polynomial"]]
    must_reject(lambda: finite.verify(changed))
    must_reject(lambda: finite.local_rational_roots(finite.polynomial(["-3/2", "3"]), []))
    must_reject(lambda: finite.local_rational_roots(finite.polynomial([-2, 4]), []))
    must_reject(lambda: finite.local_rational_roots(finite.polynomial([1, -2]), []))
    # Exercise literal ideal membership independently of an installed native solver.
    cofactor = deepcopy(rational)
    cofactor["original_cofactors"] = [["1/2", "0"], ["0", "1/3"]]
    finite.verify(cofactor)
    changed = deepcopy(cofactor)
    changed["original_cofactors"][0][0] = "0"
    must_reject(lambda: finite.verify(changed))
    changed = deepcopy(cofactor)
    changed["original_cofactors"][0].pop()
    must_reject(lambda: finite.verify(changed))
    changed = deepcopy(cofactor)
    changed["groebner_basis"][0] = "x - 5/2"
    must_reject(lambda: finite.verify(changed))
    archived = json.loads((HERE / "receipt.json").read_text())
    complete = next(row["result"] for row in archived["cases"] if row["kind"] == "rational")
    changed = deepcopy(complete)
    changed["charts"].pop()
    must_reject(lambda: rational_maps.verify_maps(changed))
    nonlinear = finite.rational_points([x**2 - y, x*y - 1], [x, y])
    nonlinear["original_cofactors"] = [["1", "0"], ["0", "1"], ["-y", "x"]]
    finite.verify(nonlinear)
    changed = deepcopy(nonlinear)
    changed["groebner_basis"].pop()
    changed["original_cofactors"].pop()
    must_reject(lambda: finite.verify(changed))
    return 18


def supplementary_controls():
    expected = json.loads((HERE / "unramified-example.json").read_text())
    if unramified_example.calculate() != expected:
        raise ArithmeticError("The exhaustive unramified-map calculation differs.")
    screens = json.loads((HERE / "frobenius-controls.json").read_text())
    for row in screens:
        frobenius_filter.verify_screen(row["source"], row["target"], row["result"])
    altered = deepcopy(screens[-1]["result"])
    altered["reductions"][-1]["source"]["counts"][0] += 1
    must_reject(lambda: frobenius_filter.verify_screen(screens[-1]["source"], screens[-1]["target"], altered))
    return {"unramified_maps": {"geometric": expected["geometric_map_count"],
                                "rational": len(expected["rational_maps"])},
            "frobenius_controls": len(screens), "altered_frobenius_receipts_rejected": 1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--rebuild", action="store_true",
                        help="Repeat discovery as well as checking every archived certificate.")
    args = parser.parse_args()
    receipt = json.loads((HERE / "receipt.json").read_text())
    rows = []
    for row in receipt["cases"]:
        counts = check_record(row, args.rebuild)
        rows.append({"id": row["id"], "geometric": counts[0], "rational": counts[1]})
    rejected = controls()
    extra = supplementary_controls()
    print(json.dumps({"status": "PASS", "cases": rows, "supplementary": extra,
                      "invalid_inputs_or_certificates_rejected": rejected}, indent=2))


if __name__ == "__main__":
    main()
