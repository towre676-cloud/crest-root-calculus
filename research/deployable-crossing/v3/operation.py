"""Finite hybrid-model guard/reachability checker, not a physical controller.

Continuous geometry/support conclusions are imported from bound receipts;
this program checks their use at discrete model transitions.  It neither
proves the supplied graph physically complete nor authorizes public opening.
"""
from collections import deque
from copy import deepcopy
from hashlib import sha256
import argparse
import json
from pathlib import Path

from evidence import canonical_sha256, digest, solve_evidence


ROOT = Path(__file__).resolve().parent
MODES = ("staged", "stabilized", "deploying", "captured", "locked", "transferred", "service", "hold", "recovery")
PHYSICAL_GATES = (
    "site_selection_and_permission", "structural_capacity_and_foundations", "accessible_route",
    "calibrated_field_measurements", "hardware_and_actuation_model", "deployment_work_zone",
    "inspection_and_emergency_procedure",
)


def _guard(value, evidence):
    status = "satisfied" if value is True else "failed" if value is False else "unresolved"
    return {"status": status, "evidence": evidence}


def _receipt_check(receipt, key, bound):
    if not bound:
        return _guard(None, "Receipt is missing, malformed or bound to another design")
    checks = receipt.get("checks")
    value = checks.get(key) if isinstance(checks, dict) else None
    if value is not None and type(value) is not bool:
        return _guard(False, "A receipt check must be a JSON boolean: " + key)
    return _guard(value, "Bound receipt checks." + key)


def _binding(receipt, design_sha256):
    if receipt is None:
        return None
    if not isinstance(receipt, dict):
        return False
    return (type(receipt.get("schema_version")) is int and receipt["schema_version"] == 1
            and receipt.get("design_sha256") == design_sha256
            and receipt.get("physical_approval") is False)


def _evidence_binding(report, design_sha256, case_id):
    try:
        if (not isinstance(report, dict) or type(report.get("schema_version")) is not int or report.get("schema_version") != 1
                or report.get("design_sha256") != design_sha256
                or type(report.get("case_count")) is not int or report["case_count"] != len(report["cases"])):
            return False, None, "Evidence report schema, design or count mismatch"
        ids, selected = set(), None
        for row in report["cases"]:
            if row["case_id"] in ids:
                return False, None, "Duplicate evidence case ID"
            ids.add(row["case_id"])
            if row["input"]["expected_bindings"]["design_sha256"] != design_sha256:
                return False, None, "Evidence case expects another design"
            actual = solve_evidence(row["input"])
            if actual != row["result"] or actual["outcome"] != row["expected_outcome"]:
                return False, None, "Evidence result fails exact replay: " + row["case_id"]
            if row["case_id"] == case_id:
                selected = actual
        if selected is None:
            return False, None, "Selected lock-evidence case is absent"
        return True, selected, "Every evidence case was replayed against its bound input"
    except (KeyError, TypeError, ValueError, AttributeError):
        return False, None, "Malformed evidence report"


def _graph(unsafe=False, omit_hold=False):
    rows = [
        {"id": "stabilize", "from": "staged", "to": "stabilized",
         "requires": ["bindings_valid", "mass_conserved", "assembly_models_solved", "continuous_support_inside"],
         "effects": {"supports_established": True}},
        {"id": "deploy", "from": "stabilized", "to": "deploying",
         "requires": ["bindings_valid", "continuous_support_inside", "stage_endpoints_match", "state_supports_established"], "effects": {}},
        {"id": "capture", "from": "deploying", "to": "captured",
         "requires": ["continuous_support_inside", "continuous_half_separation", "nominal_capture_within_envelope", "loaded_capture_within_envelope"], "effects": {}},
        {"id": "lock", "from": "captured", "to": "locked",
         "requires": ["locks_condition_established", "continuous_support_inside", "state_supports_established"],
         "effects": {"locks_secured": True}},
        {"id": "release_cantilever_restraint", "from": "locked", "to": "transferred",
         "requires": ["controlled_release_verified", "release_support_inside", "locks_condition_established", "continuous_support_inside", "state_supports_established", "state_locks_secured"],
         "effects": {"cantilever_restraint_released": True}},
        {"id": "enter_model_service", "from": "transferred", "to": "service",
         "requires": ["assembly_models_solved", "service_support_inside", "service_reference_clearance_satisfied", "stage_endpoints_match", "locks_condition_established", "state_locks_secured", "state_supports_established", "state_restraint_released"], "effects": {}},
        {"id": "begin_recovery", "from": "hold", "to": "recovery", "requires": ["recovery_path_verified"], "effects": {}},
        {"id": "finish_recovery", "from": "recovery", "to": "staged", "requires": ["recovery_path_verified"],
         "effects": {"cantilever_restraint_released": False, "locks_secured": False, "supports_established": False}},
    ]
    if not omit_hold:
        for mode in MODES:
            if mode not in ("hold", "service"):
                rows.append({"id": "hold_from_" + mode, "from": mode, "to": "hold", "requires": [], "effects": {}})
    if unsafe:
        rows.append({"id": "UNSAFE_release_before_locks_or_supports", "from": "staged", "to": "transferred",
                     "requires": [], "effects": {"cantilever_restraint_released": True}})
    return rows


def _state_tuple(state):
    return (state["mode"], state["supports_established"], state["locks_secured"], state["cantilever_restraint_released"])


def _violations(state):
    problems = []
    if state["cantilever_restraint_released"] and not (state["supports_established"] and state["locks_secured"]):
        problems.append("Cantilever restraint released before both model locks and support establishment")
    if state["mode"] == "service" and not all(state[k] for k in ("supports_established", "locks_secured", "cantilever_restraint_released")):
        problems.append("Model service requires support, lock and completed-transfer state")
    return problems


def check_graph(transitions, guards):
    """Exhaust a finite state graph and return paths, violations and deadlocks."""
    allowed_flags = {"supports_established", "locks_secured", "cantilever_restraint_released"}
    state_guards = {"state_supports_established": "supports_established", "state_locks_secured": "locks_secured", "state_restraint_released": "cantilever_restraint_released"}
    ids = set()
    for edge in transitions:
        if (edge["id"] in ids or edge["from"] not in MODES or edge["to"] not in MODES
                or not set(edge["effects"]) <= allowed_flags
                or any(type(v) is not bool for v in edge["effects"].values())
                or any(k not in guards and k not in state_guards for k in edge["requires"])):
            raise ValueError("Malformed transition or unknown guard")
        ids.add(edge["id"])
    initial = {"mode": "staged", "supports_established": False, "locks_secured": False, "cantilever_restraint_released": False}
    queue = deque([(initial, [])])
    seen = {_state_tuple(initial)}
    visited, blocked, violations, deadlocks, terminal, paths = [], [], [], [], [], {}
    checked = 0
    while queue:
        state, path = queue.popleft()
        visited.append(state)
        paths.setdefault(state["mode"], path)
        problems = _violations(state)
        if problems:
            violations.append({"state": state, "problems": problems, "witness_transition_ids": path})
            continue
        enabled = 0
        for edge in transitions:
            if edge["from"] != state["mode"]:
                continue
            checked += 1
            missing = []
            for key in edge["requires"]:
                status = ("satisfied" if state[state_guards[key]] else "failed") if key in state_guards else guards[key]["status"]
                if status != "satisfied":
                    missing.append({"guard": key, "status": status})
            if missing:
                blocked.append({"from_state": state, "transition_id": edge["id"], "unsatisfied_guards": missing})
                continue
            enabled += 1
            nxt = {**state, **edge["effects"], "mode": edge["to"]}
            key = _state_tuple(nxt)
            if key not in seen:
                seen.add(key)
                queue.append((nxt, path+[edge["id"]]))
        if not enabled:
            row = {"state": state, "witness_transition_ids": path}
            (terminal if state["mode"] in ("hold", "service") else deadlocks).append(row)
    service = "service" in paths
    return {"reachable_modes": sorted(paths), "reachable_states": visited, "witness_paths": paths,
            "blocked_transitions": blocked, "checked_transition_count": checked,
            "safety": {"declared_graph_invariants_hold": not violations, "violations": violations,
                       "scope": "All states reachable in this finite graph for this fixed guard valuation"},
            "liveness": {"model_service_state_reached": service, "accepted_hold_reached": "hold" in paths,
                         "unexpected_deadlock_count": len(deadlocks), "unexpected_deadlocks": deadlocks,
                         "accepted_terminal_states": terminal,
                         "universal_service_liveness_proved": False,
                         "meaning": "Existential reachability and deadlock reporting; accepted hold is not service, and recovery is unresolved"}}


def run_operation(evidence_report, assembly_receipt=None, motion_receipt=None, *, design_sha256=None,
                  evidence_case_id="independent_channels_synthetic", unsafe_graph=False, include_controls=True):
    """Bind real assembly/motion receipts and replay evidence before guard use.

    Content digests use canonical JSON.  This consumes recorded continuous
    conclusions; their mathematical replay belongs to the receipt generators.
    Evidence examples are synthetic guard fixtures, not a physical time log.
    """
    design_sha256 = design_sha256 or (evidence_report.get("design_sha256") if isinstance(evidence_report, dict) else None)
    digest(design_sha256, "design_sha256")
    assembly_bound, motion_bound = _binding(assembly_receipt, design_sha256), _binding(motion_receipt, design_sha256)
    geometry_hash = None
    capture_hash = None
    if assembly_bound:
        model = assembly_receipt.get("assembly")
        if not isinstance(model, dict) or model.get("design_sha256") != design_sha256:
            assembly_bound = False
        else:
            geometry_hash = canonical_sha256(model)
            capture_hash = canonical_sha256(assembly_receipt.get("capture_reference"))
            if assembly_receipt.get("geometry_and_mass_source_sha256") != geometry_hash:
                assembly_bound = False
    if motion_bound:
        if assembly_bound is None:
            motion_bound = None
        elif (not assembly_bound or motion_receipt.get("geometry_and_mass_source_sha256") != geometry_hash
              or motion_receipt.get("loaded_capture_source_sha256") != capture_hash
              or motion_receipt.get("structural_report_sha256") != canonical_sha256(assembly_receipt)):
            motion_bound = False
    evidence_bound, selected, explanation = _evidence_binding(evidence_report, design_sha256, evidence_case_id)
    binding_values = (assembly_bound, motion_bound, evidence_bound)
    valid = False if False in binding_values else None if None in binding_values else True
    guards = {"bindings_valid": _guard(valid, {"assembly": assembly_bound, "motion": motion_bound, "evidence": evidence_bound, "evidence_replay": explanation})}
    for guard, key in (("mass_conserved", "mass_conserved"), ("assembly_models_solved", "all_expected_models_solved"),
                       ("invalid_modes_have_mechanism", "invalid_modes_have_mechanism"), ("controlled_release_verified", "controlled_release_verified")):
        guards[guard] = _receipt_check(assembly_receipt, key, assembly_bound)
    for key in ("continuous_support_inside", "continuous_half_separation", "nominal_capture_within_envelope", "loaded_capture_within_envelope", "stage_endpoints_match", "release_support_inside", "service_support_inside", "service_reference_clearance_satisfied"):
        guards[key] = _receipt_check(motion_receipt, key, motion_bound)
    lock_obligations = {"full_connector_engagement", "support_state_established"}
    locks = (selected["outcome"] == "condition_established"
             and lock_obligations <= set(selected.get("required_hazards", []))) if evidence_bound and selected else None
    guards["locks_condition_established"] = _guard(locks, {"case_id": evidence_case_id, "outcome": selected["outcome"] if selected else None,
                                                             "source_kind": selected.get("source_kind") if selected else None,
                                                             "meaning": "Synthetic evidence fixture for a conditional model transition"})
    guards["recovery_path_verified"] = _guard(None, "No reverse deployment/failed-hardware recovery path is supplied")
    transitions = _graph(unsafe_graph)
    analysis = check_graph(transitions, guards)
    result = {"schema_version": 1, "design_sha256": design_sha256,
              "receipt_bindings": {"assembly_receipt_content_sha256": canonical_sha256(assembly_receipt) if assembly_receipt is not None else None,
                                   "motion_receipt_content_sha256": canonical_sha256(motion_receipt) if motion_receipt is not None else None,
                                   "evidence_report_content_sha256": canonical_sha256(evidence_report),
                                   "geometry_and_mass_source_sha256": geometry_hash,
                                   "loaded_capture_source_sha256": capture_hash},
              "guards": guards, "transitions": transitions, **analysis,
              "public_opening": {"allowed": False, "unresolved_gates": list(PHYSICAL_GATES),
                                 "meaning": "Software-model decisions never confer public opening or physical operating authorization"},
              "physical_approval": False, "evidence_is_reused_as_guard_fixture": True,
              "graph_is_physically_complete": False,
              "scope": "Finite model transitions composed with bound continuous-model receipts and synthetic evidence; no actuator commands are issued"}
    if include_controls:
        unsafe = check_graph(_graph(unsafe=True), guards)
        missing = deepcopy(guards)
        missing["continuous_support_inside"] = _guard(None, "Deliberately omitted support certificate")
        missing_run = check_graph(_graph(), missing)
        deadlock = check_graph(_graph(omit_hold=True), missing)
        release = deepcopy(guards)
        release["controlled_release_verified"] = _guard(False, "Deliberately failed release certificate")
        release_run = check_graph(_graph(), release)
        result["verification_controls"] = {
            "unsafe_release_graph": {"source_kind": "deliberately_unsafe_graph_fixture", "safety": unsafe["safety"]},
            "missing_support_certificate": {"reachable_modes": missing_run["reachable_modes"], "liveness": missing_run["liveness"]},
            "missing_hold_escape": {"liveness": deadlock["liveness"]},
            "failed_release_certificate": {"reachable_modes": release_run["reachable_modes"], "liveness": release_run["liveness"]},
        }
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--evidence", type=Path, default=ROOT/"evidence-report.json")
    parser.add_argument("--assembly", type=Path, default=ROOT/"assembly-report.json")
    parser.add_argument("--motion", type=Path, default=ROOT/"motion-report.json")
    parser.add_argument("--design", type=Path, default=ROOT/"design.json")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    read = lambda path: json.loads(path.read_text(encoding="utf-8"))
    result = run_operation(read(args.evidence), read(args.assembly), read(args.motion), design_sha256=canonical_sha256(read(args.design)))
    result["source_file_sha256"] = {"design.json": sha256(args.design.read_bytes()).hexdigest(),
                                     "evidence-report.json": sha256(args.evidence.read_bytes()).hexdigest(),
                                     "assembly-report.json": sha256(args.assembly.read_bytes()).hexdigest(),
                                     "motion-report.json": sha256(args.motion.read_bytes()).hexdigest()}
    output = json.dumps(result, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False) + "\n"
    if args.write:
        (ROOT/"operation-report.json").write_text(output, encoding="utf-8")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
