"""Focused exact evidence, transition and serialized-receipt verification."""
from copy import deepcopy
from hashlib import sha256
import argparse
import json
from pathlib import Path

from evidence import OUTCOMES, _seal, build_evidence_examples, canonical_sha256, solve_evidence
from operation import PHYSICAL_GATES, check_graph, run_operation


ROOT = Path(__file__).resolve().parent


def verify_all(check_reports=True):
    design_hash = sha256(b"explicit synthetic verification design; not a bridge").hexdigest()
    report = build_evidence_examples(design_hash)
    rows = {row["case_id"]: row for row in report["cases"]}
    checks = []

    def require(value, name):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    require({row["result"]["outcome"] for row in rows.values()} == set(OUTCOMES), "all_four_outcomes")
    require(len(rows["command_only"]["result"]["surviving_state_ids"]) == 6, "identical_command_keeps_all_six_states")
    require(rows["command_only"]["result"]["next_measurement_id"] == "engagement_probe", "independent_channel_separates_more_opposite_decisions")
    require(rows["independent_channels_synthetic"]["result"]["surviving_state_ids"] == ["full_lock__established_support"], "two_channels_isolate_full_condition")
    require(rows["impossible_reading"]["result"]["surviving_state_ids"] == [], "empty_set_is_model_mismatch")
    require(all(row["result"]["physical_authorization"] is False for row in rows.values()), "no_evidence_fixture_authorizes_physical_operation")
    base = rows["independent_channels_synthetic"]["input"]
    endpoint = deepcopy(base)
    endpoint["source"]["observations"]["engagement_probe"].update(bounds=["44", "44"], error="1")
    endpoint_result = solve_evidence(_seal(endpoint))
    require(endpoint_result["surviving_state_ids"] == ["partial_lock__established_support"], "closed_intervals_keep_exact_touching_endpoint")
    wide = rows["broad_error"]["result"]
    require(set(rows["independent_channels_synthetic"]["result"]["surviving_state_ids"]) <= set(wide["surviving_state_ids"]), "wider_error_does_not_invent_uniqueness")
    require(len(wide["surviving_state_ids"]) == 3, "broad_error_retains_all_three_lock_conditions")

    # Independent finite grid oracle for the synthetic integer signal windows.
    # Enumerating the permitted observed values is intentionally different from
    # the implementation's endpoint comparisons.
    for value in (-2, 0, 3, 43, 44, 50, 56, 97, 100, 103, 104, 200):
        sample = deepcopy(base)
        sample["source"]["observations"]["engagement_probe"].update(bounds=[str(value), str(value)], error="1")
        got = solve_evidence(_seal(sample))["surviving_state_ids"]
        observed = set(range(value-1, value+2))
        oracle = []
        for state in sample["model"]["states"]:
            if not state["hazards"]["support_state_established"]:
                continue
            lo, hi = map(int, state["predictions"]["engagement_probe"])
            if observed & set(range(lo, hi+1)):
                oracle.append(state["id"])
        require(got == oracle, "independent_integer_grid_oracle_" + str(value))

    mutators = {
        "negative_error": lambda c: c["source"]["observations"]["engagement_probe"].update(error="-1"),
        "float_error": lambda c: c["source"]["observations"]["engagement_probe"].update(error=0.5),
        "boolean_error": lambda c: c["source"]["observations"]["engagement_probe"].update(error=True),
        "nonfinite_error": lambda c: c["source"]["observations"]["engagement_probe"].update(error="NaN"),
        "reversed_bounds": lambda c: c["source"]["observations"]["engagement_probe"].update(bounds=["101", "99"]),
        "unknown_measurement": lambda c: c["source"]["observations"].update(unknown=c["source"]["observations"]["engagement_probe"]),
        "future_measurement": lambda c: c["source"]["observations"]["engagement_probe"].update(at_step=11),
        "empty_library": lambda c: c["model"].update(states=[]),
        "empty_obligations": lambda c: c["model"].update(required_hazards=[]),
        "missing_expected_hazard": lambda c: c.update(expected_hazards=["full_connector_engagement"]),
        "unknown_required_measurement": lambda c: c["model"].update(required_measurements=["unregistered"]),
        "command_substituted_for_measurement": lambda c: c["model"].update(required_measurements=["actuator_command"]),
        "command_provenance_on_probe": lambda c: c["source"]["observations"]["engagement_probe"].update(kind="command_record"),
        "synthetic_claimed_as_actual": lambda c: c["source"]["observations"]["engagement_probe"].update(kind="independent_measurement"),
        "omitted_prediction": lambda c: c["model"]["states"][0]["predictions"].pop("engagement_probe"),
        "duplicate_state_id": lambda c: c["model"]["states"][1].update(id=c["model"]["states"][0]["id"]),
        "boolean_schema": lambda c: c.update(schema_version=True),
    }
    for name, mutate in mutators.items():
        case = deepcopy(base)
        mutate(case)
        result = solve_evidence(_seal(case))
        require(result["outcome"] == "model_mismatch", "malformed_" + name)
    for bad in (None, [], {}, {"schema_version": 1}, {"schema_version": "1"}):
        require(solve_evidence(bad)["outcome"] == "model_mismatch", "malformed_top_level_" + str(len(checks)))
    tampered = deepcopy(base)
    tampered["source"]["observations"]["engagement_probe"]["bounds"] = ["0", "0"]
    require(solve_evidence(tampered)["outcome"] == "model_mismatch", "source_hash_tampering")
    stale_model = deepcopy(base)
    stale_model["source"]["model_sha256"] = "0"*64
    stale_model["bindings"]["source_sha256"] = stale_model["expected_bindings"]["source_sha256"] = canonical_sha256(stale_model["source"])
    require(solve_evidence(stale_model)["outcome"] == "model_mismatch", "source_bound_to_stale_model")

    # Artificial receipt objects test graph logic only.  Published operation
    # output below is rebuilt separately from actual assembly/motion receipts.
    assembly = {"schema_version": 1, "design_sha256": design_hash, "physical_approval": False,
                "source_kind": "synthetic_unit_test_receipt", "checks": {
                    "mass_conserved": True, "all_expected_models_solved": True,
                    "invalid_modes_have_mechanism": True, "controlled_release_verified": True}}
    assembly["assembly"] = {"design_sha256": design_hash, "fixture": "synthetic unit test geometry"}
    assembly["capture_reference"] = {"fixture": "synthetic unit test capture"}
    assembly["geometry_and_mass_source_sha256"] = canonical_sha256(assembly["assembly"])
    motion = {"schema_version": 1, "design_sha256": design_hash,
              "geometry_and_mass_source_sha256": canonical_sha256(assembly["assembly"]),
              "loaded_capture_source_sha256": canonical_sha256(assembly["capture_reference"]),
              "structural_report_sha256": canonical_sha256(assembly),
              "physical_approval": False, "source_kind": "synthetic_unit_test_receipt", "checks": {
                  "continuous_support_inside": True, "continuous_half_separation": True,
                  "nominal_capture_within_envelope": True, "loaded_capture_within_envelope": True, "stage_endpoints_match": True,
                  "release_support_inside": True, "service_support_inside": True, "service_reference_clearance_satisfied": True}}
    good = run_operation(report, assembly, motion)
    require(good["safety"]["declared_graph_invariants_hold"], "normal_graph_release_invariants")
    require(good["liveness"]["model_service_state_reached"], "model_service_reachable_with_all_conditional_guards")
    require(good["liveness"]["accepted_hold_reached"] and good["liveness"]["unexpected_deadlock_count"] == 0, "hold_is_accepted_without_being_service")
    require(good["public_opening"] == {"allowed": False, "unresolved_gates": list(PHYSICAL_GATES), "meaning": "Software-model decisions never confer public opening or physical operating authorization"}, "public_opening_remains_blocked_when_model_service_reached")
    unsafe = good["verification_controls"]["unsafe_release_graph"]["safety"]
    require(not unsafe["declared_graph_invariants_hold"] and unsafe["violations"][0]["witness_transition_ids"] == ["UNSAFE_release_before_locks_or_supports"], "unsafe_graph_concrete_release_witness")
    require(good["verification_controls"]["missing_hold_escape"]["liveness"]["unexpected_deadlock_count"] == 1, "deadlock_distinguished_from_accepted_hold")
    require("transferred" not in good["verification_controls"]["failed_release_certificate"]["reachable_modes"], "failed_release_certificate_blocks_transfer")
    for case_id in ("command_only", "partial_lock", "missing_sensor", "broad_error", "stale_reading", "wrong_design", "impossible_reading"):
        blocked = run_operation(report, assembly, motion, evidence_case_id=case_id)
        require("locked" not in blocked["reachable_modes"] and "service" not in blocked["reachable_modes"], "evidence_blocks_lock_" + case_id)
    forged = deepcopy(report)
    forged["cases"][0]["result"]["outcome"] = "condition_established"
    require(run_operation(forged, assembly, motion)["guards"]["bindings_valid"]["status"] == "failed", "forged_evidence_result_rejected")
    for filename, original in (("assembly", assembly), ("motion", motion)):
        other = deepcopy(original)
        other["design_sha256"] = "0"*64
        run = run_operation(report, other if filename == "assembly" else assembly, other if filename == "motion" else motion)
        require(run["reachable_modes"] == ["hold", "staged"], filename + "_wrong_design_blocks_progress")
    absent = deepcopy(motion)
    del absent["checks"]["continuous_support_inside"]
    require(run_operation(report, assembly, absent)["guards"]["continuous_support_inside"]["status"] == "unresolved", "missing_boolean_does_not_default_true")
    malformed = deepcopy(motion)
    malformed["checks"]["continuous_support_inside"] = "true"
    require(run_operation(report, assembly, malformed)["guards"]["continuous_support_inside"]["status"] == "failed", "string_true_does_not_satisfy_guard")
    missing = run_operation(report)
    require(missing["reachable_modes"] == ["hold", "staged"], "missing_continuous_receipts_fail_closed")
    mismatched_geometry = deepcopy(assembly)
    mismatched_geometry["assembly"]["fixture"] = "changed geometry without matching source binding"
    require(run_operation(report, mismatched_geometry, motion)["guards"]["bindings_valid"]["status"] == "failed", "whole_geometry_mass_content_binding_checked")
    mismatched_capture = deepcopy(motion)
    mismatched_capture["loaded_capture_source_sha256"] = "0"*64
    require(run_operation(report, assembly, mismatched_capture)["guards"]["bindings_valid"]["status"] == "failed", "loaded_capture_structural_source_binding_checked")
    capture_failure = deepcopy(motion)
    capture_failure["checks"]["loaded_capture_within_envelope"] = False
    require("captured" not in run_operation(report, assembly, capture_failure)["reachable_modes"], "loaded_capture_failure_blocks_nominal_capture")
    unknown_schema = deepcopy(motion)
    unknown_schema["schema_version"] = 999
    require(run_operation(report, assembly, unknown_schema)["guards"]["bindings_valid"]["status"] == "failed", "unknown_receipt_schema_fails_closed")
    wrong_obligation = deepcopy(report)
    selected_row = next(row for row in wrong_obligation["cases"] if row["case_id"] == "independent_channels_synthetic")
    selected_case = selected_row["input"]
    selected_case["expected_hazards"] = ["unrelated_condition"]
    selected_case["model"]["required_hazards"] = ["unrelated_condition"]
    for state in selected_case["model"]["states"]:
        state["hazards"] = {"unrelated_condition": True}
    selected_row["result"] = solve_evidence(_seal(selected_case))
    require(run_operation(wrong_obligation, assembly, motion)["guards"]["locks_condition_established"]["status"] == "failed", "unrelated_established_condition_cannot_discharge_lock_obligations")
    for flag, forbidden in (("release_support_inside", "transferred"), ("service_support_inside", "service"), ("service_reference_clearance_satisfied", "service")):
        for value in (False, None):
            blocked_motion = deepcopy(motion)
            if value is None:
                del blocked_motion["checks"][flag]
            else:
                blocked_motion["checks"][flag] = value
            require(forbidden not in run_operation(report, assembly, blocked_motion)["reachable_modes"], "support_clearance_guard_" + flag + "_" + str(value))
    stale_structure = deepcopy(motion)
    stale_structure["structural_report_sha256"] = "0"*64
    require(run_operation(report, assembly, stale_structure)["guards"]["bindings_valid"]["status"] == "failed", "whole_structural_report_binding_checked")

    if check_reports:
        load = lambda name: json.loads((ROOT/name).read_text(encoding="utf-8"))
        real_design = canonical_sha256(load("design.json"))
        stored_evidence = load("evidence-report.json")
        require(stored_evidence == build_evidence_examples(real_design), "published_evidence_exact_rebuild")
        stored_operation = load("operation-report.json")
        rebuilt = run_operation(stored_evidence, load("assembly-report.json"), load("motion-report.json"), design_sha256=real_design)
        source_hashes = stored_operation.pop("source_file_sha256")
        require(source_hashes == {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in ("design.json", "evidence-report.json", "assembly-report.json", "motion-report.json")}, "published_operation_source_file_hashes")
        require(stored_operation == rebuilt, "published_operation_exact_rebuild")
    return {"status": "passed", "check_count": len(checks), "checks": checks,
            "synthetic_evidence_cases": report["case_count"], "physical_measurements": 0,
            "published_receipts_checked": check_reports}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-only", action="store_true", help="Skip published-file checks while generators are being built")
    args = parser.parse_args()
    print(json.dumps(verify_all(not args.unit_only), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
