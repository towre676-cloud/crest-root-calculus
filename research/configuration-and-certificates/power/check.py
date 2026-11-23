"""Reproduce exact voltage certificates and the distinct local volt-var calculation."""

import argparse
import csv
import hashlib
import itertools
import json
import math
from fractions import Fraction as F
from pathlib import Path

HERE = Path(__file__).resolve().parent


def linear_offer_demo():
    # Inject at the two leaves of a fork. The voltage rows have been divided
    # by .02, giving 2p+q <= 3 and p+3q <= 3, p,q >= 0.
    rows = ((2, 1), (1, 3))
    primal = (F(6, 5), F(3, 5))
    dual = (F(2, 5), F(1, 5))
    assert all(a * primal[0] + b * primal[1] == 3 for a, b in rows)
    assert all(sum(dual[i] * rows[i][j] for i in range(2)) == 1 for j in range(2))
    assert sum(primal) == 3 * sum(dual)
    for p, q in itertools.product((0, primal[0]), (0, primal[1])):
        assert all(a * p + b * q <= 3 for a, b in rows)
    # An exact high-voltage AC branch violating only the *surrogate* bound.
    E, V, r = F(1), F(1059, 1000), F(1, 10)
    current = (V - E) / r
    injection = V * current
    linear = E**2 + 2 * r * injection
    limit = F(106, 100) ** 2
    assert V**2 < limit < linear and linear - V**2 == (V - E) ** 2
    return {
        "fork_optimal_leaf_injections": list(map(str, primal)),
        "fork_dual_multipliers": list(map(str, dual)),
        "fork_optimal_sum": str(sum(primal)),
        "one_line_counterexample": {
            "r": str(r),
            "active_injection": str(injection),
            "actual_v_squared": str(V * V),
            "limit_squared": str(limit),
            "linear_v_squared": str(linear),
        },
    }


def feeder_certificate():
    with (HERE / "feeder18.csv").open(newline="") as f:
        homes = list(csv.DictReader(f))
    zbase, powerbase = F(16, 25), F(250)
    transform_r, transform_z = F(".0144"), F(".06")
    xt_squared = transform_z**2 - transform_r**2
    transform_x_hi = F(".058247")
    assert F(".058246") ** 2 < xt_squared < transform_x_hi**2
    pi_hi = F(22, 7)
    names, parents, rs, xs, charging = [], [], [], [], []

    def branch(name, parent, r, x, half_charge):
        i = len(names)
        names.append(name)
        parents.append(parent)
        rs.append(r)
        xs.append(x)
        charging.append(half_charge)
        if parent >= 0:
            charging[parent] += half_charge
        return i

    def cable(name, parent, length, r_per_km, x_per_km, cnf):
        sh = pi_hi * 50 * cnf * F(1, 10**9) * length * zbase
        return branch(
            name, parent, r_per_km * length / zbase, x_per_km * length / zbase, sh
        )

    core = branch("core", -1, transform_r / 2, transform_x_hi / 2, F(0))
    main = branch("LV0", core, transform_r / 2, transform_x_hi / 2, F(0))
    for k in range(1, 11):
        main = cable(f"M{k}", main, F(".045"), F(".208"), F(".08"), 261)
        lat = main
        for j in range(1, 4):
            lat = cable(f"L{k}.{j}", lat, F(".03"), F(".642"), F(".083"), 210)
    paths = []
    for i, parent in enumerate(parents):
        paths.append({i} | (paths[parent] if parent >= 0 else set()))
    R = [[sum((rs[e] for e in u & v), F(0)) for v in paths] for u in paths]
    Xhi = [[sum((xs[e] for e in u & v), F(0)) for v in paths] for u in paths]
    source, bound = F("1.02") ** 2, F("1.06") ** 2
    existing, candidates = [F(0)] * len(names), [0] * len(names)
    for h in homes:
        i = names.index(f"L{h['main']}.{h['lateral']}")
        if h["role"] == "existing":
            existing[i] += F(h["pv_kw"])
        elif h["role"] == "candidate":
            candidates[i] += 1
        else:
            raise ValueError(h["role"])
    assert sum(candidates) == 18 and sum(existing) == 60 and len(homes) == 50
    rows = []
    for i in range(len(names)):
        h = 2 * sum(x * b for x, b in zip(Xhi[i], charging))
        c = 2 * sum(r * p for r, p in zip(R[i], existing)) / powerbase
        a = 2 * sum(r * n for r, n in zip(R[i], candidates)) / powerbase
        rows.append((a, bound * (1 - h) - source - c, h))
    # Maximize one common kW rating. Source, existing PV, and reactive upper
    # bounds are simultaneous. Distinct current/undervoltage limits are outside this LP.
    assert all(0 <= h < 1 and a > 0 and b > 0 for a, b, h in rows)
    limits = [b / a for a, b, _ in rows]
    cap = min(limits)
    i = limits.index(cap)
    assert cap < min(F(h["pv_kw"]) for h in homes if h["role"] == "candidate")
    assert all(a * cap <= b for a, b, _ in rows)
    dual = 1 / rows[i][0]
    assert rows[i][0] * dual == 1 and rows[i][1] * dual == cap
    end = names.index("L10.3")
    x_cable = Xhi[end][end] - transform_x_hi
    xt = math.sqrt(float(xt_squared))
    x_total = xt + float(x_cable)
    return {
        "source_model": "archive commit 4d311c055d98e8155d6d12d73588f747de034545; flexexport/radial.py and home_layout(32,.45,18,0)",
        "layout_sha256": hashlib.sha256(
            (HERE / "feeder18.csv").read_bytes()
        ).hexdigest(),
        "n_existing": 32,
        "n_candidates": 18,
        "existing_pv_total_kw": str(sum(existing)),
        "source_upper_v_sq": str(source),
        "voltage_limit_sq": str(bound),
        "shunt_feedback_norm_upper": str(max(h for a, b, h in rows)),
        "common_cap_kw_exact": str(cap),
        "common_cap_kw_approx": float(cap),
        "binding_bus": names[i],
        "binding_row_lhs_coefficient": str(rows[i][0]),
        "binding_row_rhs": str(rows[i][1]),
        "dual_at_binding_bus": str(dual),
        "all_feeder_rows_checked": len(rows),
        "end_bus": "L10.3",
        "transformer_x_pu_squared": str(xt_squared),
        "end_to_source_series_r_pu": str(R[end][end]),
        "end_cable_x_pu": str(x_cable),
        "transformer_fraction_of_x": xt / x_total,
        "end_to_source_r_over_x": float(R[end][end]) / x_total,
    }


def gain_checks():
    out = {}
    for name, lam in [("reported_standard", F(".477")), ("reported_steep", F("1.256"))]:
        out[name] = {
            "supplied_gain": str(lam),
            "undamped_mode": str(-lam),
            "half_relaxed_mode": str(1 - F(1, 2) * (1 + lam)),
            "relaxation_upper_exclusive": str(2 / (1 + lam)),
        }
    return out


def run():
    return {
        "rational_examples": linear_offer_demo(),
        "archived_18_home_feeder": feeder_certificate(),
        "conditional_gains": gain_checks(),
    }


if __name__ == "__main__":
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument(
        "--check", action="store_true", help="recompute and compare to receipt.json"
    )
    ap.add_argument("--write", action="store_true", help="regenerate receipt.json")
    args = ap.parse_args()
    result = run()
    path = HERE / "receipt.json"
    if args.check:
        if result != json.loads(path.read_text()):
            raise SystemExit("Mismatch in power certificate.")
        print("All rational certificates, feeder rows and gain checks reproduced.")
    if args.write:
        path.write_text(json.dumps(result, indent=2) + "\n")
    if not (args.check or args.write):
        print(json.dumps(result, indent=2))
