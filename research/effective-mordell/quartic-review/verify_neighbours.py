"""Recompute the 36 neighbour records. Pass the original solver directory as an argument."""

import argparse
import json
import sys
from pathlib import Path

from verify_review import check_maps, check_orbits, poly_in_x, require, x as pari_x


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("solver_directory", type=Path)
    args = parser.parse_args()
    directory = args.solver_directory.resolve()
    sys.path.insert(0, str(directory))
    from verify_number_fields import polynomial, recover, summary, pari, sp, x
    sources = json.loads((directory / "receipt.json").read_text())["sweep"]
    records = json.loads((directory / "number-field-receipt.json").read_text())["sweep"]
    require(len(sources) == len(records) == 36, "panel length")
    for i,(row,expected) in enumerate(zip(sources,records)):
        data = summary(*recover(polynomial(row["coefficients"])))
        require(row["source"]==expected["source"] and row["pairing"]==expected["pairing"],
                "source identification")
        require(data["geometric_maps"]==expected["geometric_maps"], "map count")
        require([m["point_field"] for m in data["maps"]]==expected["point_fields"], "map fields")
        rational = sorted(m["half_poles"] for m in data["maps"] if m["point_field"]=="Q")
        original = []
        for m in row["maps"]:
            p = sp.Poly.from_list([sp.Rational(c) for c in reversed(m["half_poles"])], x)
            original.append(str(pari(str(p.as_expr()).replace("**","^"))))
        require(rational==sorted(original), "rational half-poles")
        a = data["reciprocal_chart_center"]
        f = poly_in_x(row["coefficients"])
        if a is not None:
            f = sum(c*(a*pari_x+1)**j*pari_x**(6-j) for j,c in enumerate(row["coefficients"]))
        check_orbits(data,f)
        check_maps(data,f)
        print("Verified neighbour",i,"with",data["geometric_maps"],"maps.",flush=True)


if __name__ == "__main__":
    main()
