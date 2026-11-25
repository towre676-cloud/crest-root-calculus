"""Verify stored maps and the full geometric spectra of the Richelot controls."""

from collections import Counter
import json
from pathlib import Path
import sys

import sympy as sp

from certificate import certify_record, x

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bielliptic-reconstruction"))
from reconstruct import involutions

C2_ROOTS = (1, -1, 2, -2, 4, -4)


def pairings(roots):
    if not roots:
        yield ()
        return
    first, *tail = roots
    for point in tail:
        for rest in pairings([p for p in tail if p != point]):
            yield ((first, point), *rest)


def richelot_rhs(pairs):
    polynomials = [sp.Poly((x-a)*(x-b), x) for a, b in pairs]
    determinant = sp.Matrix([[q.nth(i) for i in range(3)] for q in polynomials]).det()
    if not determinant:
        return None
    brackets = [v.diff()*w - v*w.diff()
                for v, w in zip(polynomials[1:]+polynomials[:1], polynomials[2:]+polynomials[:2])]
    return brackets[0]*brackets[1]*brackets[2]


def expected_spectrum(degree_two_classes, bound):
    spectrum = {}
    for m in range(1, bound//2+1, 2):
        primes = sp.factorint(m)
        if any(p % 4 != 1 for p in primes):
            continue
        count = 2**len(primes)
        if degree_two_classes:
            spectrum[2*m] = degree_two_classes*count
        if degree_two_classes < 4 and 4*m <= bound:
            spectrum[4*m] = (4-degree_two_classes)*count
    return dict(sorted(spectrum.items()))


def check_maps(receipt):
    for name, case in receipt["curves"].items():
        maps = [certify_record(case["f"], row) for row in case["maps"]]
        assert len({row[2] for row in maps}) == len(maps), name
        assert all(degree == case["minimal_degree"] for degree, *_ in maps), name
        expected = Counter({"Q": 2, "Q(i)": 2}) if name.startswith(("T1:", "T2:")) else Counter({"Q": 2})
        assert Counter(row[1] for row in maps) == expected, name
    print("Verified 12 maps, their fields of definition and their distinct differential lines.")


def check_richelot(receipt):
    targets = json.loads((HERE / "targets.json").read_text(encoding="utf-8"))
    recorded = json.loads((HERE / "spectra100.json").read_text(encoding="utf-8"))
    smooth = Counter()
    degree_counts = {}
    for pairs in pairings(C2_ROOTS):
        rhs = richelot_rhs(pairs)
        if rhs is None:
            continue
        assert rhs.degree() in (5, 6) and rhs.is_sqf
        _, _, _, number = involutions(rhs.all_coeffs()[::-1])
        degree_counts[tuple(sorted(tuple(sorted(pair)) for pair in pairs))] = 2*number
        smooth[2*number] += 1
    assert smooth == {0: 4, 2: 8, 4: 1}

    for case in targets:
        pair_key = tuple(sorted(tuple(sorted(pair)) for pair in case["pairing"]))
        assert degree_counts[pair_key] == 0
        target = sp.Poly.from_list(case["f"][::-1], x, domain=sp.QQ)
        quotient = sp.cancel(target.as_expr()/richelot_rhs(case["pairing"]).as_expr())
        assert quotient.is_Rational and quotient != 0

    exact = expected_spectrum(0, 100)
    assert len(recorded) == 2
    for key, spectrum in recorded.items():
        assert {int(d): n for d, n in spectrum.items()} == exact, key
    for name, case in receipt["curves"].items():
        if name.startswith(("T1:", "T2:")):
            assert {int(d): n for d, n in case["spectrum_through_X"].items()} == expected_spectrum(0, case["X"])
            assert case["minimal_degree"] == min(exact) == 4
            assert len(case["maps"]) == exact[4] == 4
    print(f"Richelot matchings by number of degree-two classes: {dict(sorted(smooth.items()))}")
    print(f"Degenerate matchings: {15-sum(smooth.values())}. T1 and T2 have proved spectrum {exact} through 100.")


if __name__ == "__main__":
    receipt = json.loads((HERE / "replay.json").read_text(encoding="utf-8"))
    check_maps(receipt)
    check_richelot(receipt)
