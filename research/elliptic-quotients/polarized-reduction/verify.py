"""Replay the two-prime proof and controls for repeated and singular fibres."""
import json
from screen import HERE, combined, curve_profile, local_class_bound, read_terms, run


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def main():
    output = run(HERE / "quartic.json", [5, 7], 4)
    expected = json.loads((HERE / "receipt.json").read_text(encoding="utf-8"))
    require(output == expected, "The checked receipt differs from the recomputation.")
    first, second = output["reductions"]
    require(first["point_counts"] == [8, 40, 140] and second["point_counts"] == [7, 43, 352],
            "Unexpected point counts.")
    require(output["bounds"]["possible_optimal_degrees"] == [], "The all-degree exclusion failed.")
    require(first["elliptic_blocks"][0]["projector_multiplier"] == -6, "Wrong multiplier at 5.")
    require(second["elliptic_blocks"][0]["possible_special_optimal_degrees"] == [41],
            "The special degree at 7 is not uniquely forced.")
    require(local_class_bound(first, 2*5**3) == 1 and local_class_bound(first, 4*5) == 0,
            "The residue-characteristic exception at 5 was lost.")
    require(local_class_bound(second, 41*7**2) == 1 and local_class_bound(second, 7) == 0,
            "The exact prime-to-7 degree restriction was lost.")
    require(combined([first])["possible_optimal_degrees"] is None,
            "One prime alone cannot bound the exceptional prime powers.")
    row = read_terms(HERE.parent / "row-36518" / "quartic.json")
    lifted = read_terms(HERE / "quartic.json")
    require(all((row.get(k, 0)-lifted.get(k, 0)) % 5 == 0 for k in row.keys() | lifted.keys()),
            "Reduction at 5 disagrees with the known double cover.")
    positive = read_terms(HERE.parent / "primitive-degree-four" / "ciani.json")
    repeated = curve_profile(positive, 7)
    require([b["geometric_multiplicity"] for b in repeated["elliptic_blocks"]] == [3] and
            combined([repeated])["geometric_class_upper_bound"] is None,
            "The repeated-block control produced an unsupported degree-four bound.")
    try:
        curve_profile({(4, 0, 0): 1}, 5)
    except ValueError as err:
        require("Smooth reduction" in str(err), "Unexpected singular-reduction failure.")
    else:
        raise ArithmeticError("A singular quartic was accepted.")
    print("Verified both nonsimple reductions, degree 41 at 7, all-degree exclusion and controls.")


if __name__ == "__main__":
    main()
