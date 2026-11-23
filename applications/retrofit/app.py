"""Resolve two competing envelope upgrades with exact linear bounds."""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys


def number(value):
    if isinstance(value, bool):
        raise ValueError("Expected a finite decimal")
    return F(str(value))


def interval(values):
    if not isinstance(values, list) or len(values) != 2:
        raise ValueError("Supply [minimum, maximum]")
    a, b = map(number, values)
    if a > b:
        raise ValueError("Minimum cannot exceed maximum")
    return a, b


def dot(a, p):
    return a[0] * p[0] + a[1] * p[1]


def clip(poly, vec, upper):
    """Closed rational half-plane clip. Retain a point or a segment at equality."""
    out = []
    for p, q in zip(poly, poly[1:] + poly[:1]):
        u, v = dot(vec, p) - upper, dot(vec, q) - upper
        if u <= 0:
            out.append(p)
        if u * v < 0:
            s = u / (u - v)
            out.append(tuple(x + s * (y - x) for x, y in zip(p, q)))
    return list(dict.fromkeys(out))


def hull(points):
    points = sorted(set(points))
    if len(points) < 3:
        return points

    def chain(seq):
        out = []
        for p in seq:
            while len(out) > 1:
                a, b = out[-2:]
                if (b[0] - a[0]) * (p[1] - a[1]) > (b[1] - a[1]) * (p[0] - a[0]):
                    break
                out.pop()
            out.append(p)
        return out

    return chain(points)[:-1] + chain(reversed(points))[:-1]


def span(points, weights):
    vals = [dot(weights, p) for p in points]
    return min(vals), max(vals)


def axis(name):
    return {"fabric": (1, 0), "infiltration": (0, 1), "total": (1, 1)}[name]


def solve(d):
    hours = number(d["heating_degree_hours_Kh"])
    if hours <= 0:
        raise ValueError("heating_degree_hours_Kh must be positive")
    total = interval(d["total_W_K"])
    if total[0] < 0:
        raise ValueError("total_W_K must be nonnegative")
    limit = total[1]
    poly = [(F(0), F(0)), (limit, F(0)), (limit, limit), (F(0), limit)]
    measurements = [
        {"quantity": "total", "bounds_W_K": list(total)},
        *d.get("measured", []),
    ]
    for m in measurements:
        lo, hi = interval(m["bounds_W_K"])
        direction = axis(m["quantity"])
        poly = clip(clip(poly, direction, hi), tuple(-z for z in direction), -lo)
    if not poly:
        return {
            "status": "model_mismatch",
            "reason": "The measurements have no common nonnegative fabric and infiltration losses.",
        }
    actions = d["actions"]
    if len(actions) != 2 or len({a["id"] for a in actions}) != 2:
        raise ValueError("Give precisely two distinct upgrades")
    scores, options = [], []
    for a in actions:
        rates = [
            number(a[k])
            for k in ("fraction_fabric_saved", "fraction_infiltration_saved")
        ]
        cost = number(a["cost"])
        if cost <= 0 or not all(0 <= z <= 1 for z in rates):
            raise ValueError("Costs must be positive and savings fractions in [0,1]")
        energy = tuple(z * hours / 1000 for z in rates)
        score = tuple(z / cost for z in energy)
        scores.append(score)
        options.append(
            {
                "id": a["id"],
                "heat_saved_kWh_per_year": list(map(float, span(poly, energy))),
                "annual_kWh_per_cost_unit": list(map(float, span(poly, score))),
            }
        )
    advantage = tuple(x - y for x, y in zip(scores[0], scores[1]))
    lo, hi = span(poly, advantage)
    if lo > 0 or hi < 0:
        status, chosen = "one_upgrade_dominates", actions[int(hi < 0)]["id"]
    elif lo == hi == 0:
        status, chosen = "equal_value", None
    elif lo == 0 or hi == 0:
        status = "one_upgrade_never_worse"
        chosen = actions[int(hi == 0)]["id"]
    else:
        status, chosen = "needs_measurement", None
    # Two plausible homes can share a test reading exactly when their true
    # test values differ by at most twice its error. Optimize the spread in
    # score advantage over P-P intersected with that strip.
    differences = hull([(x - u, y - v) for x, y in poly for u, v in poly])
    tests = []
    for m in d.get("tests", []):
        ev, cost = number(m["error_W_K"]), number(m["cost"])
        if ev < 0 or cost <= 0:
            raise ValueError("Measurement error must be nonnegative and cost positive")
        v = axis(m["quantity"])
        residual = clip(clip(differences, v, 2 * ev), tuple(-x for x in v), 2 * ev)
        width = span(residual, advantage)[1]
        gain = hi - lo - width
        tests.append((gain / cost, gain, m["id"], float(width)))
    tests.sort(key=lambda k: (-k[0], -k[1], k[2]))
    if len({m[2] for m in tests}) != len(tests):
        raise ValueError("Test IDs must be unique")
    return {
        "status": status,
        "chosen_upgrade": chosen,
        "options": options,
        "possible_fabric_W_K": list(map(float, span(poly, (1, 0)))),
        "possible_infiltration_W_K": list(map(float, span(poly, (0, 1)))),
        "first_minus_second_kWh_per_year_per_cost": list(map(float, (lo, hi))),
        "next_test": tests[0][2]
        if status == "needs_measurement" and tests and tests[0][0] > 0
        else None,
        "test_comparison": [
            {
                "id": name,
                "worst_remaining_advantage_width": width,
                "guaranteed_width_reduction": float(gain),
                "reduction_per_test_cost": float(eff),
            }
            for eff, gain, name, width in tests
        ],
    }


def main():
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "input", type=Path, nargs="?", default=Path(__file__).with_name("example.json")
    )
    try:
        print(
            json.dumps(
                solve(json.loads(p.parse_args().input.read_text(), parse_float=str)),
                indent=2,
                allow_nan=False,
            )
        )
    except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError) as err:
        sys.exit("Input error: " + str(err))


if __name__ == "__main__":
    main()
