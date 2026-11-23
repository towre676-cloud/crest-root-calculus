"""Replay the 12 blind quartics, two classical controls and the primitive example."""
import argparse
import json
from pathlib import Path

from classify import classify


HERE = Path(__file__).resolve().parent


def run():
    root = HERE.parents[2]
    catalogue = json.loads((root / "applications/atlas/original/catalogue.json").read_text(encoding="utf-8"))
    if len(catalogue) != 12 or len({c["id"] for c in catalogue}) != 12:
        raise ArithmeticError("Original atlas population mismatch")
    controls = {"Fermat": [[4, 0, 0, 1], [0, 4, 0, 1], [0, 0, 4, 1]],
                "Klein": [[3, 1, 0, 1], [0, 3, 1, 1], [1, 0, 3, 1]]}
    cases = [(row["id"], row["terms"]) for row in catalogue]
    cases.extend(controls.items())
    primitive = json.loads((HERE.parent / "primitive-degree-four/ciani.json").read_text(encoding="utf-8"))
    cases.append(("primitive quartic cover", primitive))
    report = {name: classify(data) for name, data in cases}
    for name, expected in [("Fermat", (15, 9, 12)), ("Klein", (21, 0, 0)),
                           ("primitive quartic cover", (5, 3, 0))]:
        record = report[name]
        if tuple(record[k] for k in ("geometric_involutions", "rational_involutions",
                                     "hyperflexes_removed")) != expected:
            raise ArithmeticError(f"Control {name} disagrees")
    # Any involution fixing the degree-four map would fix its differential
    # U+V. The five axes are represented over Q or Q(sqrt(-2)).
    # The comparison is made in each orbit's certified number field.
    import sympy as sp
    T = sp.Symbol("T")
    for orbit in report["primitive quartic cover"]["orbits"]:
        # Strings here are internal output of classify(), not caller expressions.
        modulus = sp.Poly(orbit["field_polynomial"], T, domain=sp.QQ)
        l, m, n = [sp.Poly(expr, T, domain=sp.QQ) for expr in orbit["axis"]]
        if n.rem(modulus).is_zero and (l-m).rem(modulus).is_zero:
            raise ArithmeticError("Unexpected involution preserving the degree-four differential")
    known = {"Sutherland row 23640": "13997521/225", "Sutherland row 36518": "0"}
    for row in catalogue:
        name = row["id"]
        r = report[name]
        if r["geometric_involutions"] != int(name in known):
            raise ArithmeticError("Blind benchmark count mismatch")
        if name in known and r["orbits"][0]["quotient"]["j"] != known[name]:
            raise ArithmeticError("Blind quotient invariant mismatch")
    try:
        classify([[4, 0, 0, 1], [0, 4, 0, 1]])
    except ValueError as exc:
        if "singular" not in str(exc):
            raise
    else:
        raise ArithmeticError("Singular input accepted")
    return report


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute", action="store_true")
    args = parser.parse_args()
    result = run()
    if args.compute:
        print(json.dumps(result, indent=2))
    else:
        saved = json.loads((HERE / "receipt.json").read_text(encoding="utf-8"))
        if saved != result:
            raise ArithmeticError("Exact atlas replay differs from the receipt")
        print("Verified all geometric degree-two quotient classes on twelve blind rows and three controls.")


if __name__ == "__main__":
    main()
