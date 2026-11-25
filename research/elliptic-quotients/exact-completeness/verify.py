"""Replay every stored finite kernel and the absence controls."""

import argparse
import json

from absence import good_model, simplicity_certificate
from certify_spectrum import HERE, certify, data


CASES = {
    "C": ((7, 2, 13), (11, 2, 13)),
    "Cprime": ((13, 2, 13), (17, 2, 13)),
    "T1": ((17, 2, None), (29, 1, None), (61, 1, None)),
    "T2": ((17, 2, None), (29, 1, None), (61, 1, None), (89, 1, None)),
}

ABSENCE_SOURCE = [-2, 3, 1, 0, -1, 0, 1]
EXTRA_INVOLUTION_SOURCE = [1, -1, 0, 0, 0, 1, 1]
SPECTRUM = {str(d): (4 if d == 4 else 8) for d in (4, 20, 52, 68, 100, 116, 148)}


def checked_data():
    records = dict(zip(CASES, data().values(), strict=True))
    report = {"curves": {}, "absence": []}
    for name, candidates in CASES.items():
        f = records[name]["f"]
        maps = records[name]["maps"]
        results = []
        for p, r, witness in candidates:
            result = certify(f, maps, p, r, witness)
            if name.startswith("T"):
                assert result["leading_times_pi"] == "1"
                assert result["spectrum"] == SPECTRUM
                assert len(result["residue_K_v"]) == 12
                for key, value in result["residue_K_v"].items():
                    a, b = map(int, key.split(","))
                    assert value == (2 if a % 2 and b % 2 else 1)
            else:
                degree = 2 if name == "C" else 4
                assert result["kernel_order"] == degree**2
                assert result["axis_kernel_orders"] == [1, 1]
                assert result["optimal_degrees"] == [degree, degree]
            results.append(result)
        report["curves"][name] = {"f": f, "certificates": results}
        print(f"{name}: {len(results)} complete finite-kernel certificates passed.")

    f = records["C"]["f"]
    assert not good_model(f, 5)
    for f in (f, records["T2"]["f"], EXTRA_INVOLUTION_SOURCE):
        assert simplicity_certificate(f, (5, 7, 11, 13, 17, 19, 23, 29)) is None

    # At p=7 the first target has affine points beyond two-torsion, but
    # the other one does not. The implementation must return, not sample forever.
    try:
        certify(records["C"]["f"], records["C"]["maps"], 7, 1, 13)
    except ValueError as error:
        assert "No suitable fibre translate" in str(error)
    else:
        raise AssertionError("Expected the base field to lack suitable translates.")

    for p, expected in ((5, [1, 2, 5, 10, 25]), (7, [1, 1, -1, 7, 49])):
        result = simplicity_certificate(ABSENCE_SOURCE, (p,))
        assert result["weil"] == expected and result["exponent"] == 240
        report["absence"].append({"f": ABSENCE_SOURCE, **result})
    assert all(sum(c*x**i for i, c in enumerate(ABSENCE_SOURCE)) % 5 for x in range(5))
    print("Absence: certified at 5 and 7, including the model without a rational branch point at 5.")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--regenerate", action="store_true", help="Rewrite the tracked receipt")
    args = parser.parse_args()
    path = HERE / "certificates.json"
    actual = checked_data()
    if args.regenerate:
        path.write_text(json.dumps(actual, indent=2)+"\n", encoding="utf-8")
    else:
        assert actual == json.loads(path.read_text(encoding="utf-8")), "Certificate receipt differs."
