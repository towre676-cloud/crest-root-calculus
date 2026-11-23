"""Exclude rational curve maps by exact good-prime Frobenius divisibility."""
import argparse
from functools import lru_cache
import json
from pathlib import Path
import sys

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "cubic-tower-factors"))
from field import Field
from tower import weil_polynomial

from finite import polynomial

T = sp.Symbol("T")


def curve(values):
    result = polynomial(values, sp.Symbol("x"))
    if result.degree() < 3 or result.gcd(result.diff()).degree():
        raise ValueError("Use a squarefree branch polynomial of degree at least three.")
    return result


def reduction(polynomial, prime):
    """Return the same-degree squarefree model, or None at an excluded prime."""
    values = []
    for coefficient in reversed(polynomial.all_coeffs()):
        if coefficient.q % prime == 0:
            return None
        values.append(int(coefficient.p) * pow(int(coefficient.q), -1, prime) % prime)
    reduced = sp.Poly.from_list(values[::-1], polynomial.gen, modulus=prime)
    if reduced.degree() != polynomial.degree() or reduced.gcd(reduced.diff()).degree():
        return None
    return values


@lru_cache(maxsize=8)
def finite_field(prime, extension, budget):
    return Field(prime, extension, max_field_order=budget)


def point_count(values, prime, extension, budget):
    """Count the affine quadratic-character sum and the full infinity fibre."""
    field = finite_field(prime, extension, budget)
    elements = field.log[np.arange(field.q, dtype=np.int64)]
    evaluated = np.full(field.q, -1, dtype=np.int64)
    for coefficient in reversed(values):
        multiplied = np.where((elements < 0) | (evaluated < 0),
                              -1, (elements + evaluated) % field.order)
        evaluated = field.add(multiplied, field.log[coefficient])
    characters = np.where(evaluated < 0, 0, np.where(evaluated % 2 == 0, 1, -1))
    infinity_points = 1
    if (len(values) - 1) % 2 == 0:
        infinity_points += 1 if field.log[values[-1]] % 2 == 0 else -1
    return int(field.q + infinity_points + characters.sum())


def local_data(polynomial, prime, budget):
    values = reduction(polynomial, prime)
    if values is None:
        return None
    genus = (polynomial.degree() - 1) // 2
    counts = [point_count(values, prime, extension, budget)
              for extension in range(1, genus + 1)]
    traces = [prime**extension + 1 - count for extension, count in enumerate(counts, 1)]
    characteristic = weil_polynomial(traces, prime, genus)
    return {"counts": counts, "zeta_numerator": characteristic}


def screen(source_values, target_values, primes=(5, 7, 11), field_budget=500000):
    """A failed division excludes every nonconstant rational map in every degree."""
    source, target = curve(source_values), curve(target_values)
    largest_genus = (max(source.degree(), target.degree()) - 1) // 2
    reductions = []
    skipped = []
    for prime in primes:
        if type(prime) is not int or prime < 5 or not sp.isprime(prime):
            raise ValueError("Use primes at least five with the existing field-table implementation.")
        if prime**largest_genus > field_budget:
            skipped.append({"prime": prime, "reason": "field table exceeds the requested budget"})
            continue
        if reduction(source, prime) is None or reduction(target, prime) is None:
            skipped.append({"prime": prime, "reason": "the supplied model has bad reduction"})
            continue
        left, right = local_data(source, prime, field_budget), local_data(target, prime, field_budget)
        source_poly = sp.Poly.from_list(left["zeta_numerator"][::-1], T, domain=sp.QQ)
        target_poly = sp.Poly.from_list(right["zeta_numerator"][::-1], T, domain=sp.QQ)
        remainder = source_poly.rem(target_poly)
        row = {"prime": prime, "source": left, "target": right,
               "division_remainder": [str(remainder.nth(i)) for i in range(max(0, remainder.degree()) + 1)],
               "target_divides_source": remainder.is_zero}
        reductions.append(row)
        if not remainder.is_zero:
            return {"status": "rational_maps_excluded", "scope": "all nonconstant rational maps in all degrees",
                    "geometric_status": "undetermined", "witness_prime": prime,
                    "reductions": reductions, "skipped": skipped}
    return {"status": "inconclusive", "geometric_status": "undetermined",
            "reductions": reductions, "skipped": skipped}


def verify_screen(source_values, target_values, result, field_budget=500000):
    primes = [row["prime"] for row in result["reductions"]]
    rebuilt = screen(source_values, target_values, primes, field_budget)
    for key in ("status", "geometric_status", "reductions"):
        if rebuilt[key] != result[key]:
            raise ArithmeticError("The Frobenius exclusion differs from exact recounting.")
    if rebuilt.get("witness_prime") != result.get("witness_prime"):
        raise ArithmeticError("The recorded witness prime differs.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--primes", nargs="+", type=int, default=[5, 7, 11])
    parser.add_argument("--field-budget", type=int, default=500000)
    options = parser.parse_args()
    data = json.loads(options.input.read_text(encoding="utf-8"))
    print(json.dumps(screen(data["source"], data["target"], options.primes, options.field_budget), indent=2))
