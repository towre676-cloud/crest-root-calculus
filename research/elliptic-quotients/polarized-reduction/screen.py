"""Bound optimal elliptic quotient degrees using good reductions of a plane quartic.

python screen.py quartic.json --primes 5 7 --degree 4
"""
import argparse
import json
from math import lcm
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
ATLAS = HERE.parents[2] / "applications" / "atlas"
sys.path.insert(0, str(ATLAS))
sys.path.insert(0, str(ATLAS / "original/vendor/crest/scan_completeness_xxx/torsion"))
import algebra
from frobenius_fields import certificate

s = sp.Symbol("s")


def local_profile(p, weil):
    """Use a verified genus-three Weil polynomial, with descending coefficients."""
    geometry = certificate(p, weil)
    a, b, c = weil[1:4]
    real = sp.Poly(s**3 + a*s**2 + (b-3*p)*s + c-2*p*a, s)
    radical = real.sqf_part()
    blocks = []
    for component in geometry["elliptic_components"]:
        multiplicity = component["geometric_multiplicity"]
        block = {"geometric_multiplicity": multiplicity, "weil": component["weil_polynomial"]}
        if multiplicity == 1:
            coeffs = block["weil"]
            if len(coeffs) != 3 or coeffs[0] != 1 or coeffs[2] != p:
                raise ArithmeticError("Isolated elliptic block has unexpected Frobenius polynomial.")
            trace = -coeffs[1]
            complement = radical.exquo(sp.Poly(s-trace, s))
            value = int(complement.eval(trace))
            if not value:
                raise ArithmeticError("A supposedly isolated trace occurs in its complement.")
            block.update(trace=trace, real_complement=list(map(int, complement.all_coeffs())),
                         projector_multiplier=value,
                         possible_special_optimal_degrees=[
                             int(d) for d in sp.divisors(abs(value)) if d >= 2
                         ])
        blocks.append(block)
    return {"prime": p, "weil": weil, "real_weil": list(map(int, real.all_coeffs())),
            "eigenvalue_ratio_orders": geometry["root_ratio_orders"],
            "geometric_splitting_power": geometry["even_splitting_power"],
            "elliptic_blocks": blocks}


def curve_profile(form, p):
    if not sp.isprime(p) or p <= 3:
        raise ValueError("The geometric component verifier uses primes greater than three.")
    smooth = algebra.smoothness(form, p)
    if smooth is None:
        raise ValueError(f"Smooth reduction not certified at {p}.")
    counts = [algebra.point_count(form, p, exponent) for exponent in (1, 2, 3)]
    weil = algebra.characteristic(p, counts)[::-1]
    return {**local_profile(p, weil), "smoothness": smooth, "point_counts": counts}


def prime_to(n, p):
    while n % p == 0:
        n //= p
    return n


def local_class_bound(profile, degree):
    """None denotes a repeated geometric elliptic block with no finite bound here."""
    blocks = profile["elliptic_blocks"]
    if any(block["geometric_multiplicity"] > 1 for block in blocks):
        return None
    reduced = prime_to(degree, profile["prime"])
    return sum(
        any(prime_to(d, profile["prime"]) == reduced for d in block["possible_special_optimal_degrees"])
        for block in blocks
    )


def combined(profiles, degree=4):
    if not isinstance(degree, int) or degree < 2:
        raise ValueError("Optimal degrees for a smooth genus-three curve start at two.")
    if len({pr["prime"] for pr in profiles}) != len(profiles):
        raise ValueError("Choose distinct primes for independent local restrictions.")
    bounds = [local_class_bound(pr, degree) for pr in profiles]
    finite_bounds = [v for v in bounds if v is not None]
    summary = {"degree": degree, "local_class_upper_bounds": bounds,
               "geometric_class_upper_bound": min(finite_bounds, default=None),
               "possible_optimal_degrees": None}
    if any(not pr["elliptic_blocks"] for pr in profiles):
        summary["possible_optimal_degrees"] = []
        return summary
    controls = [pr for pr, bound in zip(profiles, bounds) if bound is not None]
    if len(controls) >= 2:
        # Every prime dividing n is controlled by at least one of two distinct fibres.
        first, second = controls[:2]
        def common_multiple(pr):
            return lcm(*(abs(b["projector_multiplier"]) for b in pr["elliptic_blocks"]))
        cap = common_multiple(first) * common_multiple(second)
        summary["possible_optimal_degrees"] = [
            int(n) for n in sp.divisors(cap) if n >= 2 and
            all(local_class_bound(pr, n) != 0 for pr in controls)
        ]
    return summary


def read_terms(path):
    rows = json.loads(Path(path).read_text(encoding="utf-8"))["terms"]
    form = {}
    for row in rows:
        if (len(row) != 4 or any(type(v) is not int for v in row) or
                any(v < 0 for v in row[:3]) or sum(row[:3]) != 4):
            raise ValueError("Each term must be [x power, y power, z power, integer coefficient].")
        key = tuple(row[:3])
        form[key] = form.get(key, 0) + row[3]
    return form


def run(path, primes, degree):
    form = read_terms(path)
    profiles = [curve_profile(form, p) for p in primes]
    return {"reductions": profiles, "bounds": combined(profiles, degree)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quartic", type=Path)
    parser.add_argument("--primes", nargs="+", type=int, required=True)
    parser.add_argument("--degree", type=int, default=4)
    args = parser.parse_args()
    print(json.dumps(run(args.quartic, args.primes, args.degree), indent=2))
