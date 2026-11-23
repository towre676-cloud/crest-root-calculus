"""Replay independent illustrative crossing screens from an exact JSON input.

This prints separate mathematical results, never a bridge-level safety verdict.
All cases use metres, newtons, and newton-metres. They are independent synthetic
controls, not the specification or load model of one proposed machine.
"""

import argparse
from fractions import Fraction
import hashlib
import json
from pathlib import Path

from model import (
    Interval, cantilever_uniform_load, deployment_sequence_screen, exact,
    join_geometry_screen, ramp_slope_screen, rigid_half_kinematics,
    support_reaction_screen,
)


def interval_input(value):
    if isinstance(value, list):
        if len(value) != 2:
            raise ValueError("An interval must contain exactly two endpoints")
        return Interval(*value)
    return exact(value)


def run_cases(cases, function, interval_fields=()):
    if not isinstance(cases, list) or not cases:
        raise ValueError("Each screen family needs a nonempty case list")
    result = []
    seen = set()
    for case in cases:
        if not isinstance(case, dict):
            raise ValueError("Every case must be an object")
        inputs = dict(case)
        label = inputs.pop("label")
        if not isinstance(label, str) or not label or label in seen:
            raise ValueError("Case labels must be nonempty and unique within a family")
        seen.add(label)
        for key in interval_fields:
            if key in inputs:
                inputs[key] = interval_input(inputs[key])
        result.append({"label": label, "result": function(**inputs)})
    return result


def access_arithmetic(inputs):
    names = {"rise", "declared_max_slope", "declared_max_rise_per_run",
             "declared_intermediate_landing_length"}
    if not isinstance(inputs, dict) or set(inputs) != names:
        raise ValueError("Unexpected access-arithmetic fields")
    rise, slope, maximum, landing = (
        exact(inputs[key]) for key in (
            "rise", "declared_max_slope", "declared_max_rise_per_run",
            "declared_intermediate_landing_length"))
    if rise <= 0 or slope <= 0 or maximum <= 0 or landing < 0:
        raise ValueError("Require positive rise/slope/run-rise and nonnegative landing length")
    ratio = rise / maximum
    runs = (ratio.numerator + ratio.denominator - 1) // ratio.denominator
    ramp_run = rise / slope
    intermediate = (runs - 1) * landing
    feet = Fraction("0.3048")
    return {
        "rise_m": str(rise),
        "rise_ft": str(rise / feet),
        "sloping_portions_horizontal_run_m": str(ramp_run),
        "sloping_portions_horizontal_run_ft": str(ramp_run / feet),
        "minimum_runs_under_declared_rise_limit": runs,
        "intermediate_landings": runs - 1,
        "intermediate_landing_travel_m": str(intermediate),
        "horizontal_route_travel_before_terminal_landings_m": str(ramp_run + intermediate),
        "horizontal_route_travel_before_terminal_landings_ft": str((ramp_run + intermediate) / feet),
        "scope": "Per approach; hypothetical walking-surface rise, not road clearance. "
                 "Horizontal route travel is not a straight-line parcel requirement or compliance finding.",
    }


def evaluate(data):
    required = {"format", "case_note", "units", "cantilevers", "pad_reactions",
                "kinematics", "joins", "ramps", "access_arithmetic", "deployment_record"}
    if not isinstance(data, dict) or set(data) != required:
        raise ValueError("Unexpected study fields")
    if data["format"] != "crest-root.deployable-crossing.examples/1":
        raise ValueError("Unrecognized study format")
    if data["units"] != {"length": "m", "force": "N", "moment": "N*m"}:
        raise ValueError("This example runner declares metres, newtons and newton-metres")
    if not isinstance(data["case_note"], str):
        raise ValueError("case_note must be text")
    return {
        "format": "crest-root.deployable-crossing.screens/1",
        "scope": "Independent illustrative mathematical screens; no bridge-level readiness verdict",
        "case_note": data["case_note"],
        "units": data["units"],
        "cantilevers": run_cases(data["cantilevers"], cantilever_uniform_load,
                                  ("line_load", "reach")),
        "pad_reactions": run_cases(data["pad_reactions"], support_reaction_screen,
                                    ("total_downward_load", "total_clockwise_moment")),
        "kinematics": run_cases(data["kinematics"], rigid_half_kinematics),
        "joins": run_cases(data["joins"], join_geometry_screen,
                            ("roadside_separation", "left_reach", "right_reach",
                             "left_tip_height", "right_tip_height")),
        "ramps": run_cases(data["ramps"], ramp_slope_screen,
                            ("rise", "horizontal_run")),
        "access_arithmetic": access_arithmetic(data["access_arithmetic"]),
        "deployment_record": deployment_sequence_screen(data["deployment_record"]),
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", nargs="?", type=Path,
                        default=Path(__file__).with_name("example.json"))
    args = parser.parse_args()
    try:
        raw = args.input.read_bytes()
        data = json.loads(raw)
        result = evaluate(data)
        result["input_sha256"] = hashlib.sha256(raw).hexdigest()
    except (OSError, KeyError, TypeError, ValueError) as error:
        parser.error(str(error))
    print(json.dumps(result, indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
