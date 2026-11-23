"""Keep every listed leak scenario compatible with bounded pressure-drop observations."""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
import sys


def num(x):
    if isinstance(x, bool):
        raise ValueError("Expected a finite number")
    return F(str(x))


def bounds(x):
    if not isinstance(x, list) or len(x) != 2:
        raise ValueError("An interval needs two endpoints")
    lo, hi = map(num, x)
    if lo > hi:
        raise ValueError("Reversed interval")
    return lo, hi


def prediction(row, sensor, q):
    slope = num(row["slope"][sensor])
    a, b = bounds(row.get("offset", {}).get(sensor, [0, 0]))
    p, r = sorted(slope * x for x in q)
    return p + a, r + b


def constrain(row, measurements):
    low, high = bounds(row["strength"])
    if low < 0:
        raise ValueError("Leak strength must be nonnegative")
    for sensor, obs in measurements.items():
        a, b = bounds(obs)
        slope = num(row["slope"][sensor])
        u, v = bounds(row.get("offset", {}).get(sensor, [0, 0]))
        if slope == 0:
            if a > v or b < u:
                return None
        else:
            lower, upper = sorted(((a - v) / slope, (b - u) / slope))
            low, high = max(low, lower), min(high, upper)
            if low > high:
                return None
    return low, high


def solve(d):
    sensors = d["sensors"]
    if not sensors:
        raise ValueError("Supply at least one sensor")
    observations = d.get("observations", {})
    if set(observations) - set(sensors):
        raise ValueError("Unknown sensor in observations")
    for o in observations.values():
        bounds(o)
    search = {k: num(v) for k, v in d["zones"].items()}
    if not search or any(v <= 0 for v in search.values()):
        raise ValueError("Positive disjoint search extents required")
    if d.get("extent_unit") not in ("m", "sites"):
        raise ValueError("extent_unit must be m or sites")
    if not d["scenarios"] or {row["zone"] for row in d["scenarios"]} != set(search):
        raise ValueError("Every declared zone needs at least one scenario")
    for s in sensors.values():
        if num(s["error"]) < 0 or num(s["cost"]) <= 0:
            raise ValueError("Sensor error must be nonnegative and cost positive")
    survivors = []
    for row in d["scenarios"]:
        if set(row["slope"]) != set(sensors) or set(row.get("offset", {})) - set(
            sensors
        ):
            raise ValueError(
                "Supply a response for every sensor and no unknown sensors"
            )
        for s in sensors:
            num(row["slope"][s])
            bounds(row.get("offset", {}).get(s, [0, 0]))
        q = constrain(row, observations)
        if q is not None:
            survivors.append((row, q))
    zones = sorted({r["zone"] for r, _ in survivors})
    current = sum(search[z] for z in zones)
    candidates = []
    for s, cfg in sensors.items():
        if s in observations or not survivors:
            continue
        events = {}
        err = num(cfg["error"])
        for r, q in survivors:
            lo, hi = prediction(r, s, q)
            events.setdefault(lo - err, [[], []])[0].append(r["zone"])
            events.setdefault(hi + err, [[], []])[1].append(r["zone"])
        active = {}
        worst = 0
        for x, (enter, leave) in sorted(events.items()):
            for z in enter:
                active[z] = active.get(z, 0) + 1
            worst = max(worst, sum(search[z] for z in active))
            for z in leave:
                active[z] -= 1
                if active[z] == 0:
                    del active[z]
        gain = current - worst
        candidates.append((gain / num(cfg["cost"]), gain, s, worst))
    candidates.sort(key=lambda a: (-a[0], -a[1], a[2]))
    return {
        "status": "model_mismatch"
        if not zones
        else ("single_zone" if len(zones) == 1 else "multiple_zones"),
        "zones": zones,
        "extent_unit": d["extent_unit"],
        "search_extent": float(current),
        "surviving_scenarios": [
            {"zone": row["zone"], "strength": list(map(float, q))}
            for row, q in survivors
        ],
        "next_logger": candidates[0][2]
        if candidates and candidates[0][1] > 0
        else None,
        "logger_comparison": [
            {
                "sensor": s,
                "worst_remaining_extent": float(w),
                "guaranteed_reduction": float(g),
                "reduction_per_visit_cost": float(e),
            }
            for e, g, s, w in candidates
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
    except (OSError, ValueError, KeyError, TypeError, ZeroDivisionError) as e:
        sys.exit("Input error: " + str(e))


if __name__ == "__main__":
    main()
