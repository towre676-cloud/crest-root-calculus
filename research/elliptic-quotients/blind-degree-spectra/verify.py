"""Recompute the geometric and rational elliptic degree spectra on the fixed blind slice."""
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import sympy as sp


HERE = Path(__file__).resolve().parent
REPO = HERE.parents[2]
ATLAS = REPO / "applications" / "atlas"
sys.path.insert(0, str(ATLAS))
from algebra import characteristic, parse_form, point_count, smoothness


# Seven earlier absolute exclusions use 2,3,5. The three new ones use 7.
SIMPLE_PRIMES = {42243: 2, 80491: 5, 55578: 2, 16968: 3, 78269: 2,
                 54964: 7, 23904: 7, 72418: 3, 52173: 7, 32703: 2}
ELLIPTIC_ROWS = {36518: "2*n^2 or 6*n^2", 23640: "2*n^2 or 4*n^2"}


def require(test, message):
    if not test:
        raise ArithmeticError(message)


def certify_simple(terms, p):
    require(smoothness(terms, p), f"The source is not certified smooth at {p}")
    counts = [point_count(terms, p, n) for n in (1, 2, 3)]
    T, U = sp.symbols("T U")
    ascending = characteristic(p, counts)
    h = sum(c*T**i for i, c in enumerate(ascending))
    require(sp.Poly(h, T).is_irreducible, "The Frobenius sextic factors over Q")
    resultant = sp.Poly(sp.resultant(h, h.subs(T, U*T), T), U, domain=sp.ZZ)
    unit, factors = resultant.factor_list()
    require(unit != 0, "Zero root-ratio resultant")
    nontrivial = []
    diagonal = False
    for factor, exponent in factors:
        if factor.as_expr() == U-1:
            require(exponent == 6, "Frobenius has repeated roots")
            diagonal = True
            continue
        require(factor.content() == 1 and factor.LC() > 1,
                "Need another proof to exclude a cyclotomic factor")
        nontrivial.append({"coefficients_ascending": list(map(int, reversed(factor.all_coeffs()))),
                           "multiplicity": exponent})
    require(diagonal, "Missing trivial root ratio")
    return {"p": p, "point_counts": counts, "frobenius_ascending": ascending,
            "other_irreducible_root_ratio_factors": nontrivial}


def replay_bielliptic(row, terms):
    path = HERE.parent/f"row-{row}"
    model = json.loads((path/"quartic.json").read_text(encoding="utf-8"))
    require({tuple(t[:3]): t[3] for t in model["terms"]} == terms,
            "Bielliptic certificate refers to a different quartic")
    raw = subprocess.check_output([sys.executable, "-B", str(path/"verify.py")])
    report = json.loads(raw)
    expected = json.loads((path/"receipt.json").read_text(encoding="utf-8"))
    require(report == expected, "Bielliptic receipt did not reproduce")
    geometry = report.get("geometric_degrees_all_elliptic_targets", report.get("geometric_map_degrees"))
    arithmetic = report.get("rational_degrees_all_elliptic_targets", report.get("rational_map_degrees"))
    require(report["row"] == row and geometry.split(" (")[0] == "2*n"
            and arithmetic.split(" (")[0] == ELLIPTIC_ROWS[row], "Unexpected degree spectrum")
    return {"geometric_degrees": "2*n", "rational_degrees": ELLIPTIC_ROWS[row],
            "local_certificate": f"../row-{row}/receipt.json"}


def certify_local_collision(sources, spectra):
    """Prove local agreement at 2 and geometric Hom vanishing at 3."""
    left, right = 78269, 32703
    a2 = spectra[left]["frobenius_ascending"]
    b2 = spectra[right]["frobenius_ascending"]
    require(a2 == b2, "The full Frobenius polynomials differ at 2")
    T, U = sp.symbols("T U")
    counts3 = {}
    polynomials = []
    for row in (left, right):
        terms = sources[f"Sutherland row {row}"]
        require(smoothness(terms, 3), "Need simultaneous good reduction")
        counts = [point_count(terms, 3, n) for n in (1, 2, 3)]
        counts3[row] = counts
        polynomials.append(sum(c*T**i for i, c in enumerate(characteristic(3, counts))))
    require(counts3[left][0] == counts3[right][0], "Prime-field point counts at 3 differ")
    cross = sp.Poly(sp.resultant(polynomials[0], polynomials[1].subs(T, U*T), T), U)
    primitive = cross.primitive()[1]
    require(primitive.degree() == 36 and primitive.LC() == 19683 and primitive.is_irreducible,
            "Geometric Hom exclusion via the cross-root resultant failed")
    return {"rows": [left, right], "common_frobenius_at_2_ascending": a2,
            "counts_at_3": counts3, "primitive_cross_resultant_at_3_ascending":
            list(map(int, reversed(primitive.all_coeffs()))), "geometric_Hom_rank": 0}


def main():
    data = json.loads((ATLAS/"example.json").read_text(encoding="utf-8"))
    selection = data["selection"]
    require(selection["sha256"] == "d3f06d9a92932210777499c9d995cfd618906a7c547b4c62fafa6fde4bc0a92d",
            "Unexpected Sutherland source dataset")
    require(selection["row_count"] == 82240 and selection["seed"] == "Crest-Root atlas v1",
            "Unexpected blind selection")
    cases = [e for e in data["curves"] if "row" in e]
    seed = selection["seed"]
    expected = sorted(range(1, selection["row_count"]+1),
                      key=lambda n: hashlib.sha256(f"{seed}|{n}".encode("ascii")).digest())[:12]
    require([c["row"] for c in cases] == expected and set(expected) == SIMPLE_PRIMES.keys() | ELLIPTIC_ROWS.keys(),
            "The blind rows have changed")
    # Anchor the exact equations to the original independently frozen catalogue.
    catalogue = json.loads((ATLAS/"original"/"catalogue.json").read_text(encoding="utf-8"))
    source = {e["id"]: {tuple(t[:3]): t[3] for t in e["terms"]} for e in catalogue}
    rows = []
    for c in cases:
        row = c["row"]
        terms = parse_form(c["quartic"])
        require(terms == source[c["id"]], "Catalogue and application inputs disagree")
        if row in SIMPLE_PRIMES:
            proof = certify_simple(terms, SIMPLE_PRIMES[row])
            rows.append({"row": row, "geometric_degrees": [], "rational_degrees": [], "proof": proof})
        else:
            rows.append({"row": row, **replay_bielliptic(row, terms)})
    spectra = {item["row"]: item["proof"] for item in rows if "proof" in item}
    collision = certify_local_collision(source, spectra)
    print(json.dumps({"parameters": "positive integers n", "rows": rows,
                      "local_spectrum_collision": collision}, indent=2))


if __name__ == "__main__":
    main()
