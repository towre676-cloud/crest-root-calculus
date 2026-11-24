"""Independent corner enumeration and optional pinned Safe Roads comparison."""
from copy import deepcopy
from datetime import date
from decimal import Decimal, localcontext
from fractions import Fraction
from hashlib import sha256
import argparse
import importlib.util
from itertools import product
import json
from pathlib import Path
import sys

from evidence import canonical_sha256, exact
from outcomes import OUTCOMES, ROOT, SOURCE_MODEL, UNITS, build_outcome_report, solve_outcome, synthetic_controls


def _corner_oracle(case):
    names = list(UNITS)
    values = []
    for corner in product(*(tuple(map(exact, case["inputs"][k]["bounds"])) for k in names)):
        v = dict(zip(names, corner))
        share = v["availability"]*v["uptake_when_available"]
        baseline = v["baseline_exposure"]*v["baseline_risk"]
        after = v["intervention_exposure"]*(share*v["facility_risk"]+(1-share)*v["fallback_risk"])+v["deployment_harm"]+v["displacement_harm"]
        values.append((baseline, after, baseline-after, share))
    return {name: [str(min(row[i] for row in values)), str(max(row[i] for row in values))]
            for i, name in enumerate(("baseline_expected_harm", "intervention_expected_harm", "net_harm_avoided", "facility_use_share"))}


def _source_comparison(path, controls):
    path = Path(path)
    if sha256(path.read_bytes()).hexdigest() != SOURCE_MODEL["sha256"]:
        raise AssertionError("Pinned Safe Roads source bytes do not match")
    spec = importlib.util.spec_from_file_location("safe_roads_pinned_reference", path)
    module = importlib.util.module_from_spec(spec)
    sys.modules[spec.name] = module
    spec.loader.exec_module(module)
    count = 0
    for control in controls:
        case = control["input"]
        period = module.Period(*map(date.fromisoformat, case["period"]))
        kwargs = {}
        for name, row in case["inputs"].items():
            with localcontext() as context:
                context.prec = 100
                bounds = [Decimal(exact(x).numerator)/Decimal(exact(x).denominator) for x in row["bounds"]]
            ev = row["evidence"]
            evidence = module.Evidence(ev["status"], ev["source"], date.fromisoformat(ev["as_of"]), ev["interpretation"])
            kwargs[name] = module.LocalInput(row["input_id"], module.Interval(*bounds), row["unit"], row["site_id"], period, row["outcome"], evidence)
        result = module.local_crossing_harm(**kwargs, extra_harm_disjointness_evidence=case["extra_harm_disjointness_evidence"])
        for key in ("baseline_expected_harm", "intervention_expected_harm", "net_harm_avoided", "facility_use_share"):
            actual = [str(Fraction(result[key].low)), str(Fraction(result[key].high))]
            if actual != control["result"][key]:
                raise AssertionError("Recovered Safe Roads result differs: " + key)
        count += 1
    return count


def verify_all(check_report=True, source_model=None):
    checks = []

    def require(value, name):
        if not value:
            raise AssertionError(name)
        checks.append(name)

    design = sha256(b"explicit synthetic outcomes verification design").hexdigest()
    controls = synthetic_controls(design)
    for row in controls:
        oracle = _corner_oracle(row["input"])
        require(all(row["result"][key] == value for key, value in oracle.items()), "independent_512_corner_oracle_" + row["id"])
    require(controls[0]["result"]["net_harm_avoided"] == ["58/125", "58/125"], "exact_synthetic_death_control")
    require(controls[1]["result"]["net_harm_avoided"] == ["116/25", "116/25"], "exact_separate_serious_injury_control")
    require(controls[2]["result"]["net_range_sign"] == "contains_zero", "uncertainty_can_leave_benefit_sign_unresolved")
    site = {"schema_version": 1, "design_sha256": design, "physical_approval": False,
            "sites": [{"id": "source-site-a", "label": "Source site A"}, {"id": "source-site-b", "label": "Source site B"}]}
    report = build_outcome_report(design, site)
    real = [outcome for row in report["site_results"] for outcome in row["outcomes"]]
    require(len(real) == 4 and all(x["result"]["status"] == "unknown" for x in real), "each_source_site_has_two_separate_unknown_outcomes")
    require(all(x["result"]["net_harm_avoided"] is None and len(x["result"]["missing_inputs"]) == 11 for x in real), "unknown_period_nine_inputs_and_extra_harm_boundary_are_not_zero")
    require(report["site_report_sha256"] == canonical_sha256(site), "whole_site_report_binding")
    require(all(set(x["input"]["inputs"]) == set(UNITS) for x in real), "all_nine_unknown_inputs_explicit")
    require(all(row["result"]["causal_effect_established"] is False for row in controls), "scenario_math_is_not_a_causal_effect")
    base = controls[0]["input"]
    for name in UNITS:
        partial = deepcopy(base)
        partial["inputs"][name]["bounds"] = None
        partial["inputs"][name]["evidence"]["status"] = "unknown"
        out = solve_outcome(partial)
        require(out["status"] == "unknown" and out["net_harm_avoided"] is None, "missing_" + name + "_never_defaults_zero")
    mutations = {
        "wrong_unit": lambda c: c["inputs"]["baseline_exposure"].update(unit="fraction"),
        "wrong_site": lambda c: c["inputs"]["baseline_risk"].update(site_id="different-site"),
        "mixed_outcomes": lambda c: c["inputs"]["baseline_risk"].update(outcome="serious_nonfatal_injuries"),
        "different_design": lambda c: c["inputs"]["facility_risk"].update(design_sha256="0"*64),
        "different_period": lambda c: c["inputs"]["facility_risk"].update(period=["2028-01-01", "2029-01-01"]),
        "unknown_source_model": lambda c: c.update(source_model_sha256="0"*64),
        "duplicate_id": lambda c: c["inputs"]["baseline_risk"].update(input_id=c["inputs"]["facility_risk"]["input_id"]),
        "negative_exposure": lambda c: c["inputs"]["baseline_exposure"].update(bounds=["-1", "0"]),
        "probability_over_one": lambda c: c["inputs"]["baseline_risk"].update(bounds=["0", "2"]),
        "uptake_over_one": lambda c: c["inputs"]["uptake_when_available"].update(bounds=["1", "2"]),
        "reversed_range": lambda c: c["inputs"]["facility_risk"].update(bounds=["1", "0"]),
        "float_input": lambda c: c["inputs"]["availability"].update(bounds=[0.5, 0.5]),
        "omitted_extra_harm": lambda c: c["inputs"].pop("deployment_harm"),
        "unsupported_outcome": lambda c: c.update(outcome="all_harm"),
        "empty_nonoverlap_evidence": lambda c: c.update(extra_harm_disjointness_evidence=""),
        "missing_provenance": lambda c: c["inputs"]["availability"]["evidence"].update(source=None),
        "boolean_schema": lambda c: c.update(schema_version=True),
    }
    for name, mutate in mutations.items():
        case = deepcopy(base)
        mutate(case)
        require(solve_outcome(case)["status"] == "invalid_input", "reject_" + name)
    for value in (None, [], {}, {"schema_version": 1}):
        require(solve_outcome(value)["status"] == "invalid_input", "reject_malformed_top_level_" + str(len(checks)))
    no_boundary = deepcopy(base)
    no_boundary["extra_harm_disjointness_evidence"] = None
    require(solve_outcome(no_boundary)["status"] == "unknown", "unknown_nonoverlap_blocks_computation")
    adverse = deepcopy(base)
    adverse["inputs"]["deployment_harm"]["bounds"] = ["2", "2"]
    require(solve_outcome(adverse)["net_range_sign"] == "negative", "explicit_extra_harm_can_reverse_synthetic_sign")
    bad_site = deepcopy(site)
    bad_site["design_sha256"] = "0"*64
    try:
        build_outcome_report(design, bad_site)
    except ValueError:
        require(True, "mismatched_site_report_rejected")
    else:
        raise AssertionError("Mismatched site report was accepted")
    source_comparisons = _source_comparison(source_model, controls) if source_model else 0
    if source_model:
        require(source_comparisons == 3, "three_controls_match_pinned_recovered_source")
    if check_report:
        load = lambda name: json.loads((ROOT/name).read_text(encoding="utf-8"))
        stored = load("outcome-report.json")
        hashes = stored.pop("source_file_sha256")
        require(hashes == {name: sha256((ROOT/name).read_bytes()).hexdigest() for name in ("design.json", "site-report.json")}, "published_source_file_hashes")
        require(stored == build_outcome_report(canonical_sha256(load("design.json")), load("site-report.json")), "published_outcome_report_exact_rebuild")
        require(stored["summary"]["computed_site_benefits"] == 0 and stored["summary"]["positive_site_benefit_claims"] == 0, "no_real_site_benefit_claim")
    return {"status": "passed", "check_count": len(checks), "checks": checks,
            "independent_corner_evaluations": len(controls)*512, "recovered_source_comparisons": source_comparisons,
            "source_model_sha256": SOURCE_MODEL["sha256"], "published_report_checked": check_report}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--unit-only", action="store_true")
    parser.add_argument("--source-model", type=Path, help="Optional original recovered models.py; bytes must match the pinned SHA-256")
    args = parser.parse_args()
    print(json.dumps(verify_all(not args.unit_only, args.source_model), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
