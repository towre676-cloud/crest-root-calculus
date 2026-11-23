"""Exact certificates for two four-bars. No third-party dependencies.

Run: python check.py --check receipt.json
     python check.py --write /path/to/new-receipt.json
Counts are on smooth compact curves. Only the displayed angles and j decimals
use floating point. Projective pose coordinates are encoded as integers or "inf".
"""

import argparse
import itertools
import json
import math
import random
from collections import defaultdict
from fractions import Fraction as F
from pathlib import Path


def coefficients(lengths, k=1):
    g, a, b, c = lengths
    low, high = (g - a) ** 2, (g + a) ** 2
    upper, lower = b * b + c * c + 2 * k * b * c, b * b + c * c - 2 * k * b * c
    return upper - low, upper - high, low - lower, high - lower


def invariants(lengths, k=1):
    u, v, w, z = coefficients(lengths, k)
    A, B, C = v * z, u * z + v * w, u * w
    pi, kap = A * C, 16 * k * math.prod(lengths)
    assert u * z - v * w == kap and B * B - 4 * A * C == kap**2
    disc = 16 * pi * kap**4
    if disc:
        j = F(16 * (16 * pi + kap**2) ** 3, pi * kap**4)
        # Independent binary-quartic invariants.
        I, J = B * B + 12 * A * C, 72 * A * B * C - 2 * B**3
        assert 27 * disc == 4 * I**3 - J**2
        assert j == F(6912 * I**3, 4 * I**3 - J**2)
    else:
        j = None
    return A, B, C, pi, kap, disc, j


def prime_counts(lengths):
    A, B, C, _, _, disc, _ = invariants(lengths)
    rows = []
    for p in range(3, 400, 2):
        if any(p % q == 0 for q in range(2, math.isqrt(p) + 1)) or disc % p == 0:
            continue
        # A homogeneous quartic has fibre sizes 0, 1, 2, also at infinity.
        roots = [0] * p
        for y in range(p):
            roots[y * y % p] += 1
        finite = sum(roots[(A * t**4 + B * t * t + C) % p] for t in range(p))
        rows.append([p, finite, roots[A % p], finite + roots[A % p]])
    return rows


def flips(lengths, p, expected):
    """Enumerate P1 x P1 by the original dot-product equation, then swap fibres."""
    g, a, b, c = lengths
    directions = [(i, 1) for i in range(p)] + [(1, 0)]

    def trig(n, d):
        return d * d - n * n, 2 * n * d, d * d + n * n

    trig_values = [trig(*t) for t in directions]
    poses = []
    rows, columns = defaultdict(list), defaultdict(list)
    for t, (ct, st, ht) in enumerate(trig_values):
        for u, (cu, su, hu) in enumerate(trig_values):
            closure = (
                (g * g + a * a + c * c - b * b) * ht * hu
                + 2 * g * c * ht * cu
                - 2 * g * a * ct * hu
                - 2 * a * c * (ct * cu + st * su)
            )
            if closure % p == 0:
                poses.append((t, u))
                rows[t].append(u)
                columns[u].append(t)
    assert len(poses) == expected

    def other(pair, value):
        assert value in pair and len(pair) in (1, 2)
        return pair[-1] if value == pair[0] else pair[0]

    t, u = poses[0]
    result = []
    while (t, u) not in result:
        result.append((t, u))
        u = other(rows[t], u)
        t = other(columns[u], t)
    assert (t, u) == result[0]
    return [[x if x != p else "inf" for x in pair] for pair in result]


def box_margins(lengths, k, error):
    """Worst folded and extended margins in an independent positive length box."""
    g, a, b, c = lengths
    if error < 0 or error >= min(lengths) or not 0 <= k <= 1:
        raise ValueError("The error box must preserve positive lengths and 0 ≤ k ≤ 1.")
    distance = max(0, abs(g - a) - 2 * error)
    max_fold = max(
        bb * bb + cc * cc - 2 * k * bb * cc
        for bb, cc in itertools.product((b - error, b + error), (c - error, c + error))
    )
    b0, c0 = b - error, c - error
    return (
        distance**2 - max_fold,
        b0 * b0 + c0 * c0 + 2 * k * b0 * c0 - (g + a + 2 * error) ** 2,
    )


def wing(lengths, witness_prime, other_prime):
    abc = invariants(lengths)
    A, B, C, pi, kap, disc, j = abc
    rows = prime_counts(lengths)
    counts = {row[0]: row[-1] for row in rows}
    gcd = math.gcd(*(r[-1] for r in rows))
    # Two good primes suffice even using only prime-to-p torsion reduction.
    bound = math.gcd(counts[witness_prime], counts[other_prime])
    assert gcd == bound == 4
    assert counts[other_prime] % witness_prime and counts[witness_prime] % other_prime
    orbit = flips(lengths, witness_prime, counts[witness_prime])
    assert gcd % len(orbit) != 0
    assert not math.isqrt(pi) ** 2 == pi and B > kap > 0
    # A positive transmission-angle certificate with fixed rational cosine.
    good, bad = F(0), F(min(lengths), 2)
    k = F(19, 20)
    assert min(box_margins(lengths, k, good)) > 0
    assert min(box_margins(lengths, k, bad)) < 0
    for _ in range(32):
        mid = (good + bad) / 2
        if min(box_margins(lengths, k, mid)) > 0:
            good = mid
        else:
            bad = mid
    return {
        "lengths_tenths_mm": lengths,
        "alpha0_alpha1_beta0_beta1": coefficients(lengths),
        "A_B_C": [A, B, C],
        "Pi": pi,
        "kappa": kap,
        "binary_quartic_discriminant": disc,
        "j_rational": str(j),
        "j_approx": float(j),
        "good_prime_rows_p_finite_infinity_total": rows,
        "gcd_all_orders": gcd,
        "two_prime_torsion_bound": [
            [p, counts[p]] for p in (witness_prime, other_prime)
        ],
        "rational_torsion": "Z/2",
        "flip_witness_prime": witness_prime,
        "double_flip_orbit_before_return": orbit,
        "order_mod_p": len(orbit),
        "translation_order_over_Q": "infinite",
        "uniform_error_at_cos19_over20_mm": {
            "passing_lower": str(good / 10),
            "wall_or_failure_upper": str(bad / 10),
        },
    }


def random_checks():
    rng = random.Random(8675309)
    n = 0
    while n < 400:
        # Primitive integer tuples are rational tuples after independent denominators.
        lengths = tuple(F(rng.randint(1, 70), rng.randint(1, 9)) for _ in range(4))
        if invariants(lengths)[5]:
            n += 1
    counts = {"strict_grashof": 0, "non_grashof": 0, "singular": 0}
    for _ in range(20000):
        while True:
            lengths = [rng.randint(1, 90) for _ in range(4)]
            x, y, z, w = sorted(lengths)
            if w < x + y + z:
                break
        u, v, s, t = coefficients(lengths)
        pi = u * v * s * t
        assert (pi > 0) - (pi < 0) == (x + w < y + z) - (x + w > y + z)
        counts[
            "strict_grashof" if pi > 0 else "non_grashof" if pi < 0 else "singular"
        ] += 1
    return {
        "seed": 8675309,
        "smooth_rational_invariant_tests": n,
        "quadrilateral_tests": counts,
    }


def main_receipt():
    g, a, b, c = 100, 30, 75, F("54.99994")
    u, v, w, z = coefficients((g, a, b, c))
    half_width = 2 * math.degrees(math.atan(math.sqrt(-v / u)))
    missed = sum(
        abs((degree / 2 - 180.25 + 180) % 360 - 180) < half_width
        for degree in range(720)
    )
    assert missed == 0 and w > 0 and z > 0
    residues = sorted(
        {
            ((15 * d * d + 7 * n * n) * (d * d + 9 * n * n)) % 16
            for n, d in itertools.product(range(16), repeat=2)
            if n % 2 or d % 2
        }
    )
    assert residues == [7, 12, 15] and set(residues).isdisjoint(
        {i * i % 16 for i in range(16)}
    )
    # b²+c²=25=(g+a)², so the k→0 limit is the finite special value 1728.
    exceptional = (3, 2, 3, 4)
    assert invariants(exceptional, 0)[5] == 0
    assert abs(invariants(exceptional, F(1, 1000))[-1] - 1728) < F(1, 50)
    return {
        "nominal_100_percent": wing((590, 293, 529, 700), 19, 29),
        "printable_159_percent": wing((938, 466, 841, 1113), 17, 43),
        "real_full_turn_without_rational_pose": {
            "lengths": [2, 1, 2, 2],
            "primitive_quartic_residues_mod16": residues,
        },
        "hidden_failure": {
            "lengths_mm": ["100", "30", "75", "54.99994"],
            "ground_rotation_deg": 0.25,
            "outer_margin_mm": str(b + c - g - a),
            "failure_width_deg": 2 * half_width,
            "failure_interval_deg": [180.25 - half_width, 180.25 + half_width],
            "violations_seen_on_grid_of_720": missed,
        },
        "transmission_finite_limit_example": {"lengths": exceptional, "j_limit": 1728},
        "regression_checks": random_checks(),
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--check", type=Path)
    mode.add_argument("--write", type=Path)
    args = parser.parse_args()
    data = json.dumps(main_receipt(), indent=2) + "\n"
    if args.check:
        if json.loads(data) != json.loads(args.check.read_text()):
            raise SystemExit("Receipt differs.")
        print("All four-bar checks match the receipt.")
    elif args.write:
        args.write.write_text(data)
    else:
        print(data, end="")
