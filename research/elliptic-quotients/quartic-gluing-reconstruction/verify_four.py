"""Replay both maps and the full signed kernel on the Richelot control and sweep."""
from __future__ import annotations

import argparse
from collections import Counter
from importlib.util import module_from_spec, spec_from_file_location
import json
from pathlib import Path
import sys

import sympy as sp

from reconstruct_four import odd_chart, pari, solve, t, x

HERE = Path(__file__).resolve().parent
TORSION = HERE.parents[1] / "torsion-denominators"
SOURCE = 11*(x*x-10*x+13)*(2*x*x+2*x-7)*(7*x*x+18*x+3)


def read_polynomial(values, variable):
    return sp.Poly.from_list([sp.Rational(c) for c in reversed(values)], variable, domain=sp.QQ)


def short_curve(cubic):
    return pari.ellinit([pari(str(cubic.nth(k))) for k in (1, 0)])


def targets(record):
    return [read_polynomial(record["cubic"], t) for record in record["maps"]]


def check_doubling(result):
    """Check all order-four blocks double to their recorded order-two labels."""
    g1, g2 = targets(result)
    labels = set(tuple(map(sp.Rational, row)) for row in result["order_two_x_pairs"])
    hit = Counter()
    for component in result["order_four"]:
        factor = read_polynomial(component["factor"], t)
        x1 = sp.Poly(t, t, domain=sp.QQ)
        x2 = read_polynomial(component["x"], t)
        multiplier = read_polynomial(component["y_multiplier"], t)
        assert (g1*multiplier**2-g2.compose(x2)).rem(factor).is_zero

        def double(cubic, coordinate):
            square_y = cubic.compose(coordinate).rem(factor)
            numerator = (coordinate**2*3 + cubic.nth(1))**2
            return (numerator*(square_y*4).invert(factor)-coordinate*2).rem(factor)

        roots = double(g1, x1), double(g2, x2)
        assert all(root.degree() <= 0 for root in roots)
        pair = tuple(root.nth(0) for root in roots)
        assert pair in labels
        hit[pair] += factor.degree()
    assert hit == {pair: 2 for pair in labels}


def projective_relation(left, right):
    """Return a rational projectivity relating two abscissa functions."""
    a, b, c, d = sp.symbols("a b c d")
    numerator_left, denominator_left = sp.fraction(sp.cancel(left))
    numerator_right, denominator_right = sp.fraction(sp.cancel(right))
    relation = (a*numerator_left+b*denominator_left)*denominator_right - (
        c*numerator_left+d*denominator_left
    )*numerator_right
    matrix, _ = sp.linear_eq_to_matrix(sp.Poly(sp.expand(relation), x).all_coeffs(), a, b, c, d)
    kernel = matrix.nullspace()
    if len(kernel) != 1:
        return None
    coefficients = kernel[0]
    assert coefficients[0]*coefficients[3] != coefficients[1]*coefficients[2]
    return list(coefficients)


def abscissa(record):
    p, d = (read_polynomial(record[key], x) for key in ("numerator", "half_poles"))
    return p.as_expr()/d.as_expr()**2


def local_replay(source, record):
    sys.path.insert(0, str(TORSION))
    from jac import geometrically_distinct

    spec = spec_from_file_location("torsion_verifier", TORSION / "verify.py")
    verifier = module_from_spec(spec)
    spec.loader.exec_module(verifier)
    existing = json.loads((TORSION / "certificates.json").read_text())
    targets_for_count = [short_curve(g) for g in targets(record)]
    coefficients = list(map(int, reversed(source.all_coeffs())))
    for p in (37, 61):
        certificate = next(c for c in existing if c["curve"] == coefficients and c["p"] == p)
        assert verifier.check_certificate(certificate) == [4, 4]
        traces = [int(e.ellap(p)) for e in targets_for_count]
        assert sorted(traces) == sorted(block["trace"] for block in certificate["blocks"])
        assert geometrically_distinct(*traces, p)
        print(f"Independent local certificate at {p}: traces {traces}, degrees [4, 4]")
    return [37, 61]


def check_control():
    source = sp.Poly(SOURCE, x, domain=sp.QQ)
    result = solve(source)
    check_doubling(result)
    assert [m["half_poles"] for m in result["maps"]] == [["-3", "0", "1"], ["-5", "0", "1"]]
    archive = json.loads((TORSION / "map4.json").read_text())
    old_x = sp.sympify(archive["u"], locals={"z": x})
    old_y_times_y = sp.sympify(archive["Y_times_w"], locals={"z": x})
    assert sp.cancel(old_y_times_y**2-SOURCE*(old_x-1)*(old_x-4)*(old_x-9)) == 0
    new = result["maps"][1]
    change = projective_relation(abscissa(new), old_x)
    assert change is not None
    h = read_polynomial(new["lift"], x).as_expr()
    s = read_polynomial(new["half_poles"], x).as_expr()
    ratio = sp.cancel(sp.diff(old_x, x)*h*SOURCE/(old_y_times_y*sp.diff(abscissa(new), x)*s**3))
    assert ratio.is_Rational and ratio != 0
    result["richelot_abscissa_change"] = list(map(str, change))
    result["richelot_differential_ratio"] = str(ratio)
    result["independent_local_primes"] = local_replay(source, result)
    # Move x=2 to infinity. The search receives only the transformed equation.
    transformed = odd_chart(source, sp.QQ.convert(2), 6, sp.QQ)
    moved = solve(transformed)
    check_doubling(moved)
    assert sorted(str(short_curve(g).j()) for g in targets(result)) == sorted(
        str(short_curve(g).j()) for g in targets(moved)
    )
    moved_functions = [abscissa(m).subs(x, 1/(x-2)) for m in moved["maps"]]
    assert all(sum(projective_relation(abscissa(m), n) is not None for n in moved_functions) == 1
               for m in result["maps"])
    print("C': both maps and 16 kernel points verified in the original and reciprocal charts.")
    return result


def certificate_of_distinction(source, result):
    from jac import geometrically_distinct

    curves = [short_curve(g) for g in targets(result)]
    if curves[0].j() == curves[1].j():
        return None
    for prime in map(int, sp.primerange(5, 150)):
        if source.LC() % prime == 0 or source.discriminant() % prime == 0:
            continue
        if any(e.disc() % prime == 0 for e in curves):
            continue
        traces = [int(e.ellap(prime)) for e in curves]
        if geometrically_distinct(*traces, prime):
            return {"prime": prime, "traces": traces}
    raise ArithmeticError("No prime in range separates the two elliptic factors.")


def check_sweep():
    inputs = json.loads((TORSION / "sweep.json").read_text())
    minima = json.loads((TORSION / "exact_sweep.json").read_text())["rows"]

    def key(row):
        return tuple(row["source"]), tuple(map(tuple, row["pairing"]))

    degree_four = {key(row) for row in minima if row["minimum"] == 4}
    result = []
    for row in inputs:
        if key(row) not in degree_four:
            continue
        f = read_polynomial(row["coeffs"], x)
        reconstructed = solve(f)
        check_doubling(reconstructed)
        distinction = certificate_of_distinction(f, reconstructed)
        assert (distinction is None) == (row["source"] == [1, 2, 4])
        result.append({"source": row["source"], "pairing": row["pairing"],
                       "coefficients": row["coeffs"], "geometric_distinction": distinction,
                       **reconstructed})
        print(f"{row['source']} {row['pairing']}: both maps and the full kernel", flush=True)
    assert len(result) == len(degree_four) == 36
    assert sum(r["geometric_distinction"] is not None for r in result) == 32
    # Send a simple pole of a known paired fibre to infinity, exercising deg(S)=1.
    row = next(r for r in result if r["source"] == [2, 3, 5]
               and r["pairing"] == [[2, 3], [-2, -5], [-3, 5]])
    f = read_polynomial(row["coefficients"], x)
    changed = solve(odd_chart(f, sp.QQ.convert(1), 6, sp.QQ))
    assert sorted(len(m["half_poles"]) for m in changed["maps"]) == [2, 3]
    check_doubling(changed)
    return result


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--sweep", action="store_true", help="Replay the 36 degree-four neighbours.")
    parser.add_argument("--write", action="store_true", help="Regenerate the stored receipt (with --sweep).")
    args = parser.parse_args()
    if args.write and not args.sweep:
        parser.error("--write requires --sweep")
    result = {"control": check_control()}
    path = HERE / "receipt.json"
    if args.sweep:
        result["sweep"] = check_sweep()
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    else:
        receipt = json.loads(path.read_text())
        assert all(receipt[k] == data for k, data in result.items())
