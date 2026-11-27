"""Replay the number-field reconstruction of optimal degree-four elliptic maps."""
from __future__ import annotations

import argparse
import json
import sys
from pathlib import Path

import sympy as sp

from number_field_fibres import X, differential, pairing_orbits, pairwise_distinct, pari, recover, summary, x
from etale_fibres import recover as recover_rational

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "exact-completeness"))
from absence import simplicity_certificate  # noqa: E402

CURVES = {
    # Repeated-factor curves; the finite-kernel spectrum gives four degree-four classes.
    "T1": [-1792, 7200, 5268, -900, -1605, -90, 154],
    "T2": [-9856, 1440, 6420, 900, -1317, -450, 28],
    # The Richelot control, complete with two classes by geometric non-isogeny.
    "Cprime": [-3003, -14850, 11968, 6600, -3575, -990, 154],
    # Only degree-two optimal classes.
    "C": [-36, 0, 49, 0, -14, 0, 1],
    # Irreducible sextic with cyclically permuted pairs.
    "PERMUTED": [1, -12, 12, 14, -11, -6, 1],
    # Galois group S6.
    "S6": [-1, -1, 0, 0, 0, 0, 1],
    # PERMUTED after a rational change of chart; one map has a simple pole at a finite point,
    # so the univariate pole chart deg S = 1 carries a genuine solution.
    "PERMUTED_CHART": [1, 0, -26, -70, -57, -14, -1],
    # prod_{j=0}^{5} (x - j/2): rational non-integral coefficients for the resolvent.
    "HALVES": ["0", "-15/4", "137/8", "-225/8", "85/4", "-15/2", "1"],
    # Supplied by C: maps whose pairing is defined over Q(sqrt(-3)) and not over Q.
    "X6PLUS1": [1, 0, 0, 0, 0, 0, 1],
}
EXPECTED_GEOMETRIC = {"T1": 4, "T2": 4, "Cprime": 2, "C": 0, "PERMUTED": 2, "S6": 0,
                      "PERMUTED_CHART": 2, "HALVES": 0, "X6PLUS1": 12}


def polynomial(coeffs):
    return sp.Poly.from_list([sp.Rational(c) for c in coeffs[::-1]], x, domain=sp.QQ)


def rational_signature(cover):
    """Half-poles and target j-invariant of a rational-solver cover."""
    g = cover.cubic
    a4, a6 = g.nth(1), g.nth(0)
    return str(pari(str(cover.half_poles.as_expr()).replace("**", "^"))), str(1728*4*a4**3/(4*a4**3+27*a6**2))


def check_curve(name, coeffs):
    f = polynomial(coeffs)
    result = recover(f)
    data = summary(*result)
    covers = result[4]
    assert data["geometric_maps"] == EXPECTED_GEOMETRIC[name], (name, data["geometric_maps"])
    assert sum(p["orbit_size"] for p in data["pairing_orbits"]) == 15
    # The rational part agrees with the rational solver of etale_fibres.py.
    rational = sorted((m["half_poles"], m["j"]) for m in data["maps"] if m["point_field"] == "Q")
    expected = sorted(rational_signature(c) for c, _ in recover_rational(f)[2])
    assert rational == expected, (name, rational, expected)
    # The pulled-back differentials span the holomorphic forms when there are two directions.
    forms = [differential(c) for c in covers]
    record = {"name": name, "coefficients": coeffs, **data}
    if name in ("T1", "T2"):
        # Two rational maps with differentials proportional to 1 and x, and one conjugate pair
        # over Q(i). Its differential is a + b x with (a/b)^2 rational and negative: the diagonal
        # directions use an isomorphism between the two quadratic twists defined over Q(i).
        fields = sorted(m["point_field"] for m in data["maps"])
        assert fields == ["Q", "Q", "w^2 + 1"]
        conjugate = next(c for c in covers if c.field.degree == 2)
        omega = differential(conjugate)
        ratio = pari.polcoef(omega, 0, X)/pari.polcoef(omega, 1, X)
        square = pari.simplify(pari.lift(ratio**2))
        assert pari.poldegree(square) <= 0 and square < 0
        rational_forms = sorted(int(pari.poldegree(w, X)) for c, w in zip(covers, forms) if c.field.degree == 1)
        assert rational_forms == [0, 1]
        record["diagonal_ratio_squared"] = str(square)
    if name == "PERMUTED_CHART":
        assert sorted(m["j"] for m in data["maps"]) == ["1792", "48384"]
        assert sorted(m["half_poles"] for m in data["maps"]) == ["x + 1/2", "x^2 - 2*x - 1"]
    if name == "HALVES":
        assert [p["orbit_size"] for p in data["pairing_orbits"]] == [1]*15
    if name == "X6PLUS1":
        assert all(m["point_field"] == "w^2 - w + 1" and m["j"] == "54000" for m in data["maps"])
        unstable = [m for m in data["maps"] if m["pairing_field"] != "Q"]
        assert len(data["maps"]) == 6 and len(unstable) == 4
        assert all(m["pairing_field"] == "w^2 - w + 1" for m in unstable)
    if name == "S6":
        certificate = simplicity_certificate(coeffs, primes=range(3, 30))
        assert certificate is not None and certificate["p"] == 13
        record["absolute_simplicity"] = certificate
    return record


def check_disjointness():
    """The cross-orbit test accepts the certified pairings and rejects a duplicated orbit."""
    _, orbits = pairing_orbits(polynomial(CURVES["T1"]))
    weights = pairwise_distinct(orbits)
    assert weights is not None and sum(p.orbit for p in orbits) == 15
    assert pairwise_distinct(orbits + orbits[1:2]) is None
    assert pairwise_distinct(orbits[:1] + orbits[:1]) is None
    return {"accepted_weights": weights, "orbit_sizes": [p.orbit for p in orbits]}


def check_sweep():
    """All 36 degree-four Richelot neighbours: rational maps unchanged, geometric counts."""
    rows = json.loads((HERE / "receipt.json").read_text())["sweep"]
    result = []
    for row in rows:
        f = polynomial(row["coefficients"])
        data = summary(*recover(f))
        rational = sorted(m["half_poles"] for m in data["maps"] if m["point_field"] == "Q")
        archived = sorted(str(pari(str(sp.Poly.from_list([sp.Rational(c) for c in m["half_poles"][::-1]], x)
                                       .as_expr()).replace("**", "^"))) for m in row["maps"])
        assert rational == archived, (row["source"], rational, archived)
        expected = 4 if row["source"] == [1, 2, 4] else 2
        assert data["geometric_maps"] == expected, (row["source"], data["geometric_maps"])
        result.append({"source": row["source"], "pairing": row["pairing"],
                       "geometric_maps": data["geometric_maps"],
                       "point_fields": [m["point_field"] for m in data["maps"]]})
        print(f"{row['source']} {row['pairing']}: {data['geometric_maps']} geometric maps", flush=True)
    assert len(result) == 36
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", action="store_true", help="Also replay the 36 degree-four neighbours.")
    parser.add_argument("--write", action="store_true", help="Regenerate the stored receipt.")
    args = parser.parse_args()
    result = {"curves": [check_curve(name, coeffs) for name, coeffs in CURVES.items()],
              "disjointness": check_disjointness()}
    for record in result["curves"]:
        print(f"{record['name']}: {record['geometric_maps']} geometric maps", flush=True)
    if args.sweep:
        result["sweep"] = check_sweep()
    path = HERE / "number-field-receipt.json"
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        receipt = json.loads(path.read_text())
        assert all(receipt[k] == v for k, v in result.items()), "Receipt mismatch."
    print("number-field reconstruction verified")
