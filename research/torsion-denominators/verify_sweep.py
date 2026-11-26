"""Determine the geometric minimum for every non-product neighbour in the sweep."""

import argparse
from collections import Counter
from pathlib import Path
import json
import subprocess
import sys

import sympy as sp

from jac import geometrically_distinct
from richelot import pairings, richelot, x

HERE = Path(__file__).resolve().parent
QUOTIENTS = HERE.parent / "elliptic-quotients"
KERNEL_STUDY = QUOTIENTS / "exact-completeness"
sys.path.insert(0, str(QUOTIENTS / "bielliptic-reconstruction"))
from reconstruct import involutions

SOURCES = ((1, 2, 3), (1, 2, 4), (1, 2, 5), (1, 3, 4), (2, 3, 5))
REPEATED = (1, 2, 4)


def quotient_traces(source, p):
    """Trace of the two displayed rational degree-two quotients at a good odd p."""
    roots = tuple(a*a % p for a in source)
    if 0 in roots or len(set(roots)) != 3:
        raise ValueError("The elliptic quotients have bad reduction.")
    squares = {a*a % p for a in range(1, p)}
    first, second = 0, -1  # One point at infinity on the cubic, two on the monic quartic.
    for a in range(p):
        cubic = sp.prod((a-b) for b in roots) % p
        quartic = a*cubic % p
        first -= (1 if cubic in squares else -1) if cubic else 0
        second -= (1 if quartic in squares else -1) if quartic else 0
    return int(first), int(second)


def separation_prime(source):
    for p in map(int, sp.primerange(3, 100)):
        try:
            traces = quotient_traces(source, p)
        except ValueError:
            continue
        if geometrically_distinct(*traces, p):
            return {"prime": p, "traces": list(traces), "power": 120}
    raise ArithmeticError("No certificate of geometric distinction found for the source factors.")


def normalized_neighbour(pairs):
    quadratics = [(x-a)*(x-b) for a, b in pairs]
    delta, factors = richelot(quadratics)
    if not delta:
        return None
    model = sp.Poly(delta*sp.prod(factors), x, domain=sp.ZZ).primitive()[1]
    if model.degree() not in (5, 6) or not model.is_sqf:
        raise ArithmeticError("A non-product neighbour is not a smooth genus-two model.")
    return int(delta), list(map(int, reversed(model.all_coeffs())))


def certified_kernel_models():
    subprocess.run([sys.executable, "-B", str(KERNEL_STUDY / "verify.py")], check=True)
    data = json.loads((KERNEL_STUDY / "certificates.json").read_text())
    records = {}
    for name in ("T1", "T2"):
        rec = data["curves"][name]
        if any(min(map(int, certificate["spectrum"])) != 4 for certificate in rec["certificates"]):
            raise ArithmeticError("An independent finite-kernel certificate has changed.")
        f = rec["f"]
        records[tuple(f)] = name
        records[tuple((-1)**i*c for i, c in enumerate(f))] = f"{name} with x replaced by -x"
    return records


def run():
    old = json.loads((HERE / "sweep.json").read_text())
    history = {(tuple(row["source"]), tuple(map(tuple, row["pairing"]))): row for row in old}
    repeated_degree_four = certified_kernel_models()
    seen, recovered = set(), []
    rows, sources = [], []
    for source in SOURCES:
        witness = None if source == REPEATED else separation_prime(source)
        source_rows = []
        signed_roots = [value for a in source for value in (a, -a)]
        for pairs in pairings(signed_roots):
            neighbour = normalized_neighbour(pairs)
            if neighbour is None:
                continue
            delta, coeffs = neighbour
            key = (source, tuple(pairs))
            previous = history[key]
            if coeffs != previous["coeffs"] or delta != previous["Delta"]:
                raise ArithmeticError("A previously swept model differs from the exact Richelot model.")
            seen.add(key)
            _, _, _, count = involutions(coeffs)
            if count:
                degree = 2
                certificate = "geometric elliptic involution"
            elif witness is not None:
                degree = 4
                certificate = "distinct source factors, Richelot degree-four map and no elliptic involution"
            else:
                degree = 4
                certificate = f"finite-kernel spectrum for {repeated_degree_four[tuple(coeffs)]}"
            if previous["certified_degrees"] not in (None, [degree, degree]):
                raise ArithmeticError("The exact minimum differs from a stored torsion result.")
            entry = {"source": list(source), "pairing": [list(pair) for pair in pairs],
                     "geometric_elliptic_involutions": count, "minimum": degree,
                     "certificate": certificate}
            if previous["certified_degrees"] is None:
                recovered.append(entry)
            source_rows.append(entry)
        degrees = Counter(str(row["minimum"]) for row in source_rows)
        sources.append({"source": list(source), "separation": witness, "degrees": dict(degrees)})
        rows.extend(source_rows)
    if seen != set(history):
        raise ArithmeticError("The exact and archived Richelot sweeps cover different neighbours.")
    return {"sources": sources, "rows": rows, "previously_inconclusive": recovered}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="Write the exact sweep receipt")
    args = parser.parse_args()
    path = HERE / "exact_sweep.json"
    result = run()
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    elif result != json.loads(path.read_text()):
        raise ArithmeticError("The exact sweep differs from the stored receipt.")
    counts = Counter(row["minimum"] for row in result["rows"])
    print(f"Certified all {sum(counts.values())} neighbours. Minima: {dict(sorted(counts.items()))}.")
    print(f"Resolved {len(result['previously_inconclusive'])} previously inconclusive minima.")


if __name__ == "__main__":
    main()
