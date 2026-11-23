"""Replay measured pressures from operating networks against independently simulated leaks."""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path
from app import solve

HERE = Path(__file__).parent
ERROR = F("0.025")


def observation(run, s):
    k = run["sensors"].index(s)
    value = F(str(run["drop_bar"][k])) - F(str(run["drop_bar"][0]))
    return [float(value - ERROR), float(value + ERROR)]


def make_case(run, guard=F(0)):
    sensors = run["sensors"][1:]
    offsets = {
        s: [float(F(str(lo)) - guard), float(F(str(hi)) + guard)]
        for s, (lo, hi) in zip(sensors, run["background_differential_bar"][1:])
    }
    return {
        "title": "Measured controlled discharge " + run["date"],
        "extent_unit": "sites",
        "pressure_unit": "bar difference relative to " + run["sensors"][0],
        "zones": {z: 1 for z in run["hydraulic_drop_bar"]},
        "sensors": {s: {"error": float(ERROR), "cost": 1} for s in sensors},
        "observations": {sensors[0]: observation(run, sensors[0])},
        "scenarios": [
            {
                "zone": z,
                "strength": [1, 1],
                "slope": dict(
                    zip(sensors, [float(F(str(v)) - F(str(ds[0]))) for v in ds[1:]])
                ),
                "offset": offsets,
            }
            for z, ds in run["hydraulic_drop_bar"].items()
        ],
    }


def evaluate(run, guard=F(0)):
    d = make_case(run, guard)
    available = list(d["sensors"])
    first = next(iter(d["observations"]))
    initial = solve(d)
    adaptive = initial["next_logger"] or available[1]

    def observe(s):
        d["observations"] = {t: observation(run, t) for t in (first, s)}
        return solve(d)["zones"]

    fixed = observe(available[1])
    selected = observe(adaptive)
    distances = []
    for row in d["scenarios"]:
        total = 0
        for sensor, obs in d["observations"].items():
            model = row["slope"][sensor] + sum(row["offset"][sensor]) / 2
            total += abs(model - sum(obs) / 2)
        distances.append((total, row["zone"]))
    return {
        "event": run["date"],
        "truth": run["true_zone"],
        "reference": run["sensors"][0],
        "first": first,
        "second": adaptive,
        "ranked_choice": initial["next_logger"],
        "initial": initial["zones"],
        "selected": selected,
        "fixed": fixed,
        "nearest_baseline": min(distances)[1],
    }


def summarize(events):
    summary = {}
    for method in ("initial", "selected", "fixed"):
        summary[method] = {
            "true_site_retained": sum(e["truth"] in e[method] for e in events),
            "total_sites_to_check": sum(len(e[method]) for e in events),
            "empty_sets": sum(not e[method] for e in events),
        }
    summary["nearest_correct_after_selected"] = sum(
        e["truth"] == e["nearest_baseline"] for e in events
    )
    return summary


def calibration(runs):
    # Smallest single outward allowance that retains all known true sites on
    # every training logger independently, above the stated measurement error.
    margin = F(0)
    for run in runs:
        d = make_case(run)
        row = next(r for r in d["scenarios"] if r["zone"] == run["true_zone"])
        for s in d["sensors"]:
            obs_lo, obs_hi = map(lambda x: F(str(x)), observation(run, s))
            center = F(str(row["slope"][s]))
            lo, hi = map(lambda x: F(str(x)), row["offset"][s])
            margin = max(margin, center + lo - obs_hi, obs_lo - (center + hi))
    return margin


def replay(source):
    runs = source["events"]
    events = [evaluate(r) for r in runs]
    training_day = runs[0]["date"][:10]
    training = [r for r in runs if r["date"][:10] == training_day]
    held = [r for r in runs if r["date"][:10] != training_day]
    guard = calibration(training)
    return {
        "events_scored": len(events),
        "excluded_events": source["excluded_events"],
        "summary": summarize(events),
        "events": events,
        "later_day": {
            "training_day": training_day,
            "calibration_events": len(training),
            "held_out_events": len(held),
            "additional_model_halfwidth_bar": float(guard),
            "before_calibration": summarize([evaluate(r) for r in held]),
            "after_calibration": summarize([evaluate(r, guard) for r in held]),
        },
    }


if __name__ == "__main__":
    p = argparse.ArgumentParser(description=__doc__)
    p.add_argument(
        "--event",
        type=int,
        help="print an ordinary app input for the 1-based scored event",
    )
    args = p.parse_args()
    d = json.loads((HERE / "field.json").read_text())
    if args.event is not None and not 1 <= args.event <= len(d["events"]):
        p.error("event outside replay")
    print(
        json.dumps(
            replay(d) if args.event is None else make_case(d["events"][args.event - 1]),
            indent=2,
        )
    )
