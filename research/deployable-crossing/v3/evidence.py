"""Finite, exact evidence compatibility; no physical controller or opening authority.

Predictions are closed rational intervals.  A recorded interval [a,b] and
absolute error e permits [a-e,b+e].  Compatibility uses a single candidate
state across every observed channel.  Independent interval nuisance values
are allowed; correlations must instead be encoded in separate states.
"""
from copy import deepcopy
from fractions import Fraction
from hashlib import sha256
from itertools import product
import argparse
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
OUTCOMES = ("condition_established", "condition_fails", "ambiguous", "model_mismatch")
SOURCE_KINDS = {"synthetic_fixture", "independent_measurement", "mixed"}
OBSERVATION_KINDS = {"synthetic_fixture", "independent_measurement", "command_record"}


def canonical_sha256(value):
    """Hash JSON content, distinct from the SHA-256 of a serialized file."""
    return sha256(json.dumps(value, sort_keys=True, separators=(",", ":"),
                             ensure_ascii=False, allow_nan=False).encode("utf-8")).hexdigest()


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise ValueError("Exact scalars must be integers or decimal/fraction strings")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise ValueError("Invalid finite rational scalar") from exc


def interval(value, name):
    if not isinstance(value, list) or len(value) != 2:
        raise ValueError(name + " must be a two-entry closed interval")
    lo, hi = map(exact, value)
    if lo > hi:
        raise ValueError(name + " has reversed endpoints")
    return lo, hi


def digest(value, name):
    if not isinstance(value, str) or len(value) != 64 or any(c not in "0123456789abcdef" for c in value):
        raise ValueError(name + " must be a lowercase SHA-256 digest")
    return value


def names(value, name, nonempty=True):
    if (not isinstance(value, list) or (nonempty and not value)
            or any(not isinstance(v, str) or not v for v in value)
            or len(value) != len(set(value))):
        raise ValueError(name + " must contain unique nonempty names")
    return value


def _parse(case):
    if not isinstance(case, dict) or type(case.get("schema_version")) is not int or case.get("schema_version") != 1:
        raise ValueError("Evidence case schema_version must be 1")
    model, source = case["model"], case["source"]
    bound, expected = case["bindings"], case["expected_bindings"]
    keys = {"design_sha256", "model_sha256", "source_sha256"}
    if not isinstance(bound, dict) or not isinstance(expected, dict) or set(bound) != keys or set(expected) != keys:
        raise ValueError("Complete design/model/source bindings are required")
    for key in keys:
        digest(bound[key], key)
        digest(expected[key], "expected " + key)
    if bound != expected:
        raise ValueError("Supplied bindings differ from expected bindings")
    if canonical_sha256(model) != bound["model_sha256"]:
        raise ValueError("Model content hash mismatch")
    if canonical_sha256(source) != bound["source_sha256"]:
        raise ValueError("Source record content hash mismatch")
    if model["design_sha256"] != bound["design_sha256"] or source["design_sha256"] != bound["design_sha256"]:
        raise ValueError("Design hash mismatch")
    if source["model_sha256"] != bound["model_sha256"]:
        raise ValueError("Source refers to a different model")
    if source["kind"] not in SOURCE_KINDS or not isinstance(source["source_id"], str) or not source["source_id"]:
        raise ValueError("Explicit source identity and provenance kind are required")
    step = case["evaluation_step"]
    if type(step) is not int or step < 0:
        raise ValueError("evaluation_step must be a nonnegative integer")
    hazards = names(model["required_hazards"], "required_hazards")
    if names(case["expected_hazards"], "expected_hazards") != hazards:
        raise ValueError("The model omits or changes the expected hazard obligations")
    sensors = model["sensors"]
    if not isinstance(sensors, dict) or not sensors:
        raise ValueError("At least one sensor/record definition is required")
    parsed_sensors = {}
    for key, row in sensors.items():
        if not isinstance(key, str) or not key or not isinstance(row["unit"], str) or not row["unit"]:
            raise ValueError("Measurement identities and units must be nonempty strings")
        limit, planned, cost = exact(row["max_error"]), exact(row["planning_error"]), exact(row["cost"])
        if limit < 0 or planned < 0 or planned > limit or cost <= 0:
            raise ValueError("Measurement errors must be bounded/nonnegative and cost positive")
        age = row["max_age_steps"]
        if type(age) is not int or age < 0:
            raise ValueError("max_age_steps must be a nonnegative integer")
        allowed = names(row["allowed_observation_kinds"], "allowed_observation_kinds")
        if not set(allowed) <= OBSERVATION_KINDS:
            raise ValueError("Unknown observation provenance kind")
        parsed_sensors[key] = {"definition": row, "max_error": limit,
                               "planning_error": planned, "cost": cost}
    required = names(model["required_measurements"], "required_measurements")
    if not set(required) <= set(sensors):
        raise ValueError("Unknown required measurement")
    if any("command_record" in sensors[key]["allowed_observation_kinds"] for key in required):
        raise ValueError("A command record cannot be a required independent observation")
    raw_states = model["states"]
    if not isinstance(raw_states, list) or not raw_states:
        raise ValueError("An explicit nonempty candidate library is required")
    states, state_ids = [], set()
    for state in raw_states:
        state_id = state["id"]
        if not isinstance(state_id, str) or not state_id or state_id in state_ids:
            raise ValueError("Candidate state IDs must be unique nonempty strings")
        state_ids.add(state_id)
        if set(state["hazards"]) != set(hazards) or any(type(v) is not bool for v in state["hazards"].values()):
            raise ValueError("Every state must explicitly discharge or fail every declared hazard")
        if set(state["predictions"]) != set(sensors):
            raise ValueError("Every state must predict every defined measurement")
        predictions = {k: interval(v, "prediction " + k) for k, v in state["predictions"].items()}
        states.append({"id": state_id, "hazards": state["hazards"], "predictions": predictions,
                       "condition": all(state["hazards"].values())})
    observations = source["observations"]
    if not isinstance(observations, dict) or not set(observations) <= set(sensors):
        raise ValueError("Observations must use known measurement IDs")
    parsed_observations = {}
    for key, row in observations.items():
        sensor = parsed_sensors[key]
        if row["unit"] != sensor["definition"]["unit"]:
            raise ValueError("Measurement unit mismatch: " + key)
        if row["kind"] not in sensor["definition"]["allowed_observation_kinds"]:
            raise ValueError("Observation kind cannot satisfy this channel: " + key)
        if source["kind"] == "synthetic_fixture" and row["kind"] == "independent_measurement":
            raise ValueError("A synthetic source cannot claim an actual independent measurement")
        if source["kind"] == "independent_measurement" and row["kind"] == "synthetic_fixture":
            raise ValueError("An actual measurement source cannot silently include a synthetic fixture")
        error = exact(row["error"])
        if error < 0 or error > sensor["max_error"]:
            raise ValueError("Measurement error outside declared bounds: " + key)
        at = row["at_step"]
        if type(at) is not int or at < 0 or at > step or step-at > sensor["definition"]["max_age_steps"]:
            raise ValueError("Stale or future observation: " + key)
        lo, hi = interval(row["bounds"], "observation " + key)
        parsed_observations[key] = {"interval": (lo-error, hi+error), "kind": row["kind"],
                                    "raw_interval": (lo, hi), "error": error, "at_step": at}
    return model, source, bound, hazards, parsed_sensors, required, states, parsed_observations


def _overlap(a, b):
    return max(a[0], b[0]) <= min(a[1], b[1])


def solve_evidence(case):
    """Return one of four outcomes, including fail-closed malformed-input results.

    Expected hashes are caller-provided trust anchors, not authenticated facts.
    Success is universal over a supplied finite library, not over all physical
    states.  Required measurements cannot be replaced by command records.
    """
    result = {"schema_version": 1, "outcome": "model_mismatch", "surviving_state_ids": [],
              "missing_measurement_ids": [], "reasons": [], "excluded_states": [],
              "scope": "Finite explicit model and bounded intervals only; no physical opening authority",
              "physical_authorization": False}
    try:
        model, source, bound, hazards, sensors, required, states, observations = _parse(case)
    except (KeyError, TypeError, ValueError, AttributeError, OverflowError) as exc:
        result["reasons"] = ["invalid_or_unbound_evidence: " + str(exc)]
        return result
    result.update(bindings=deepcopy(bound), source_kind=source["kind"],
                  source_id=source["source_id"], required_hazards=hazards,
                  declared_candidate_count=len(states), observed_measurement_ids=sorted(observations),
                  observation_provenance={k: v["kind"] for k, v in observations.items()},
                  effective_observation_intervals={k: list(map(str, v["interval"])) for k, v in observations.items()})
    survivors = []
    for state in states:
        violations = [{"measurement_id": key, "prediction": list(map(str, state["predictions"][key])),
                       "permitted_observation": list(map(str, row["interval"]))}
                      for key, row in observations.items() if not _overlap(state["predictions"][key], row["interval"])]
        if violations:
            result["excluded_states"].append({"state_id": state["id"], "disjoint_intervals": violations})
        else:
            survivors.append(state)
    result["surviving_state_ids"] = [s["id"] for s in survivors]
    missing = sorted(set(required)-set(observations))
    result["missing_measurement_ids"] = missing
    if not survivors:
        result["reasons"].append("No declared candidate explains every observation; empty sets never pass")
    elif missing:
        result["outcome"] = "ambiguous"
        result["reasons"].append("Required evidence is absent: " + ", ".join(missing))
    elif all(s["condition"] for s in survivors):
        result["outcome"] = "condition_established"
        result["reasons"].append("Every compatible candidate discharges every declared hazard")
    elif all(not s["condition"] for s in survivors):
        result["outcome"] = "condition_fails"
        result["reasons"].append("Every compatible candidate fails at least one declared hazard")
    else:
        result["outcome"] = "ambiguous"
        result["reasons"].append("Compatible candidates disagree about the declared condition")
    result["hazard_witnesses"] = {hazard: {"discharging_state_ids": [s["id"] for s in survivors if s["hazards"][hazard]],
                                          "failing_state_ids": [s["id"] for s in survivors if not s["hazards"][hazard]]}
                                   for hazard in hazards}
    good, bad = [s for s in survivors if s["condition"]], [s for s in survivors if not s["condition"]]
    ranked = []
    for key in sorted(set(sensors)-set(observations)):
        error = sensors[key]["planning_error"]
        separated = 0
        for a, b in product(good, bad):
            pa, pb = a["predictions"][key], b["predictions"][key]
            if not _overlap((pa[0]-error, pa[1]+error), (pb[0]-error, pb[1]+error)):
                separated += 1
        if good and bad:
            ranked.append({"measurement_id": key, "opposite_decision_pairs_separated": separated,
                           "opposite_decision_pairs": len(good)*len(bad), "cost": str(sensors[key]["cost"]),
                           "planning_error": str(error), "separated_pairs_per_cost": str(Fraction(separated)/sensors[key]["cost"])})
    ranked.sort(key=lambda r: (-exact(r["separated_pairs_per_cost"]), -r["opposite_decision_pairs_separated"], r["measurement_id"]))
    result["ranked_next_measurements"] = ranked
    result["next_measurement_id"] = ranked[0]["measurement_id"] if ranked and ranked[0]["opposite_decision_pairs_separated"] else None
    return result


def _seal(case):
    """Bind synthetic fixtures; not an authentication or calibration function."""
    model_hash = canonical_sha256(case["model"])
    case["source"]["model_sha256"] = model_hash
    binding = {"design_sha256": case["model"]["design_sha256"], "model_sha256": model_hash,
               "source_sha256": canonical_sha256(case["source"])}
    case["bindings"] = binding
    case["expected_bindings"] = deepcopy(binding)
    return case


def _fixture(design_sha256):
    digest(design_sha256, "design_sha256")
    sensors = {}
    for key, unit, limit, planned, kinds in [
            ("actuator_command", "mm", "0", "0", ["command_record"]),
            ("engagement_probe", "mm", "100", "1", ["synthetic_fixture", "independent_measurement"]),
            ("support_reaction_probe", "N", "5000", "100", ["synthetic_fixture", "independent_measurement"])]:
        sensors[key] = {"unit": unit, "max_error": limit, "planning_error": planned, "cost": "1",
                        "max_age_steps": 2, "allowed_observation_kinds": kinds}
    states = []
    for lock, supports in product(("full", "partial", "missing"), ("established", "missing")):
        states.append({"id": lock + "_lock__" + supports + "_support",
                       "hazards": {"full_connector_engagement": lock == "full", "support_state_established": supports == "established"},
                       "predictions": {"actuator_command": ["1000", "1000"],
                                       "engagement_probe": {"full": ["98", "102"], "partial": ["45", "55"], "missing": ["0", "2"]}[lock],
                                       "support_reaction_probe": ["20000", "22000"] if supports == "established" else ["0", "2000"]}})
    model = {"model_id": "synthetic-lock-and-support-v1", "design_sha256": design_sha256,
             "required_hazards": ["full_connector_engagement", "support_state_established"],
             "required_measurements": ["engagement_probe", "support_reaction_probe"], "sensors": sensors, "states": states,
             "scope": "Authored synthetic signal intervals, not predictions calibrated to the reference bridge"}
    return {"schema_version": 1, "evaluation_step": 10, "expected_hazards": model["required_hazards"][:], "model": model,
            "source": {"kind": "synthetic_fixture", "source_id": "authored-lock-ambiguity-control",
                       "design_sha256": design_sha256, "model_sha256": "", "observations": {}}}


def _observation(value, unit, error="0", kind="synthetic_fixture"):
    return {"bounds": [value, value], "unit": unit, "error": error, "kind": kind, "at_step": 10}


def build_evidence_examples(design_sha256):
    base = _fixture(design_sha256)
    command = _observation("1000", "mm", kind="command_record")
    probe = _observation("100", "mm", "1")
    support = _observation("21000", "N", "100")
    rows = []

    def add(case_id, description, observations, expected, mutate=None, corrupt=None):
        case = deepcopy(base)
        case["source"]["source_id"] += ":" + case_id
        case["source"]["observations"] = deepcopy(observations)
        if mutate:
            mutate(case)
        _seal(case)
        if corrupt:
            corrupt(case)
        result = solve_evidence(case)
        if result["outcome"] != expected:
            raise AssertionError((case_id, result["outcome"], expected))
        rows.append({"case_id": case_id, "description": description, "input": case,
                     "expected_outcome": expected, "result": result})

    add("command_only", "The same commanded position fits full, partial and missing locks, with either support state", {"actuator_command": command}, "ambiguous")
    full = {"actuator_command": command, "engagement_probe": probe, "support_reaction_probe": support}
    add("independent_channels_synthetic", "Two independently defined synthetic channels isolate the full-lock/established-support candidate", full, "condition_established")
    add("partial_lock", "The engagement observation rejects the full-lock state", {**full, "engagement_probe": _observation("50", "mm", "1")}, "condition_fails")
    add("missing_lock", "The engagement observation supports a missing connector", {**full, "engagement_probe": _observation("1", "mm", "1")}, "condition_fails")
    add("missing_support_state", "A full lock reading cannot compensate for missing support establishment", {**full, "support_reaction_probe": _observation("1000", "N", "100")}, "condition_fails")
    add("missing_sensor", "A missing required support observation blocks establishment", {"actuator_command": command, "engagement_probe": probe}, "ambiguous")
    add("broad_error", "A bounded but broad error admits full, partial and missing locks", {**full, "engagement_probe": _observation("50", "mm", "55")}, "ambiguous")
    add("impossible_reading", "No library state predicts this engagement reading", {**full, "engagement_probe": _observation("200", "mm", "1")}, "model_mismatch")
    add("stale_reading", "A stale support reading is rejected", full, "model_mismatch", mutate=lambda c: c["source"]["observations"]["support_reaction_probe"].update(at_step=1))
    add("bad_unit", "Units must match the declared channel exactly", full, "model_mismatch", mutate=lambda c: c["source"]["observations"]["engagement_probe"].update(unit="m"))
    add("excess_error", "An observation cannot exceed the declared error bound", full, "model_mismatch", mutate=lambda c: c["source"]["observations"]["engagement_probe"].update(error="101"))
    add("wrong_design", "A source record bound to another assembly is rejected", full, "model_mismatch", mutate=lambda c: c["source"].update(design_sha256="0"*64))
    add("changed_model_bytes", "Changed model bytes without a matching trusted binding are rejected", full, "model_mismatch", corrupt=lambda c: c["model"]["states"][0]["predictions"]["engagement_probe"].__setitem__(0, "99"))
    add("omitted_hazard", "A candidate cannot silently omit a declared hazard", full, "model_mismatch", mutate=lambda c: c["model"]["states"][0]["hazards"].pop("support_state_established"))
    return {"schema_version": 1, "design_sha256": design_sha256, "source_kind": "synthetic_fixture", "case_count": len(rows),
            "physical_measurement_count": 0, "cases": rows,
            "claims": {"physical_authorization": False, "actual_bridge_sensor_validation": False,
                       "finite_library_exhaustive_physical_state_space": False},
            "meaning": "Exact compatibility examples demonstrate ambiguity and fail-closed behavior; all signal intervals are authored fixtures"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", nargs="?", type=Path)
    parser.add_argument("--design", type=Path, default=ROOT/"design.json")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    result = solve_evidence(json.loads(args.case.read_text(encoding="utf-8"))) if args.case else build_evidence_examples(canonical_sha256(json.loads(args.design.read_text(encoding="utf-8"))))
    output = json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False) + "\n"
    if args.write:
        (ROOT/"evidence-report.json").write_text(output, encoding="utf-8")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
