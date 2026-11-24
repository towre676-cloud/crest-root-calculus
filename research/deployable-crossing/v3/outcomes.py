"""Exact local harm accounting adapted from the recovered Safe Roads model.

This is scenario accounting, not a causal estimator.  Missing site inputs
remain unknown; no geometric, structural or simulated operational pass is
converted into a public-health benefit.  Deaths and serious nonfatal injuries
are separate outcomes and are never added here.
"""
from copy import deepcopy
from datetime import date
from hashlib import sha256
import argparse
import json
from pathlib import Path

from evidence import canonical_sha256, digest, exact, interval


ROOT = Path(__file__).resolve().parent
SOURCE_MODEL = {
    "archive": "Safe_Roads_Recovery_2026-11-23 (revision 2)",
    "path": "gap_closure/model_repair/models.py",
    "function": "local_crossing_harm",
    "sha256": "d2c6a478d8276f376367ae80dfff1e439a9a214844e22d93ce3413979283412f",
    "adaptation": "Self-contained Fraction arithmetic for the same interval-accounting formula; explicit unknown inputs and design/site bindings added",
    "attribution": "Inspected recovered Safe Roads research source; no empirical effect is imported",
}
OUTCOMES = ("deaths", "serious_nonfatal_injuries")
UNITS = {
    "baseline_exposure": "person_crossings", "intervention_exposure": "person_crossings",
    "availability": "fraction", "uptake_when_available": "fraction",
    "baseline_risk": "harm_events/person_crossing", "facility_risk": "harm_events/person_crossing",
    "fallback_risk": "harm_events/person_crossing", "deployment_harm": "harm_events",
    "displacement_harm": "harm_events",
}


def _label(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ValueError("Explicit " + name + " required")
    return value


def _period(value):
    if value is None:
        return None
    if not isinstance(value, list) or len(value) != 2 or any(not isinstance(v, str) for v in value):
        raise ValueError("Period must be a half-open pair of ISO dates")
    start, end = map(date.fromisoformat, value)
    if start >= end or [start.isoformat(), end.isoformat()] != value:
        raise ValueError("A canonical positive-length date period is required")
    return value


def _parse(case):
    if not isinstance(case, dict) or type(case.get("schema_version")) is not int or case["schema_version"] != 1:
        raise ValueError("Outcome case schema_version must be 1")
    design = digest(case["design_sha256"], "design_sha256")
    site = _label(case["site_id"], "site_id")
    outcome = case["outcome"]
    if outcome not in OUTCOMES:
        raise ValueError("Deaths and serious nonfatal injuries must be modeled separately")
    if case["source_model_sha256"] != SOURCE_MODEL["sha256"]:
        raise ValueError("Accounting source model hash mismatch")
    period = _period(case["period"])
    missing = [] if period is not None else ["period"]
    boundary = case["extra_harm_disjointness_evidence"]
    if boundary is None:
        missing.append("extra_harm_disjointness_evidence")
    else:
        _label(boundary, "deployment/displacement/main-site nonoverlap evidence")
    inputs = case["inputs"]
    if not isinstance(inputs, dict) or set(inputs) != set(UNITS):
        raise ValueError("All nine input fields are required; unknown values must be explicit nulls")
    parsed, ids, statuses = {}, set(), []
    for name, unit in UNITS.items():
        row = inputs[name]
        input_id = _label(row["input_id"], "input_id")
        if input_id in ids:
            raise ValueError("Duplicate local input identifier")
        ids.add(input_id)
        if row["unit"] != unit:
            raise ValueError("Exposure, event and proportion units cannot be interchanged: " + name)
        if (row["site_id"], row["outcome"], row["period"], row["design_sha256"]) != (site, outcome, period, design):
            raise ValueError("Local inputs must share design, site, period and outcome: " + name)
        evidence = row["evidence"]
        if evidence["status"] not in ("unknown", "observed", "estimated", "scenario"):
            raise ValueError("Unknown evidence status")
        _label(evidence["interpretation"], "input interpretation")
        if row["bounds"] is None:
            if evidence["status"] != "unknown":
                raise ValueError("An absent input must be labeled unknown: " + name)
            missing.append(name)
            continue
        if evidence["status"] == "unknown":
            raise ValueError("Unknown evidence cannot support numeric bounds: " + name)
        _label(evidence["source"], "input evidence source")
        if not isinstance(evidence["as_of"], str) or date.fromisoformat(evidence["as_of"]).isoformat() != evidence["as_of"]:
            raise ValueError("Explicit canonical evidence date required")
        lo, hi = interval(row["bounds"], name)
        if lo < 0:
            raise ValueError("Local input bounds cannot be negative: " + name)
        if unit in ("fraction", "harm_events/person_crossing") and hi > 1:
            raise ValueError("Fractions and single-person event probabilities cannot exceed one: " + name)
        parsed[name] = lo, hi
        statuses.append(evidence["status"])
    return parsed, missing, statuses


def solve_outcome(case):
    result = {"status": "invalid_input", "missing_inputs": [], "baseline_expected_harm": None,
              "intervention_expected_harm": None, "net_harm_avoided": None, "facility_use_share": None,
              "physical_approval": False, "causal_effect_established": False,
              "interpretation": "Input-range accounting enclosure, not a confidence interval, causal estimate or measured benefit"}
    try:
        p, missing, statuses = _parse(case)
    except (KeyError, TypeError, ValueError, AttributeError) as exc:
        result["reason"] = str(exc)
        return result
    result.update(design_sha256=case["design_sha256"], site_id=case["site_id"], outcome=case["outcome"],
                  period=case["period"], missing_inputs=missing,
                  source_model_sha256=SOURCE_MODEL["sha256"], input_content_sha256=canonical_sha256(case))
    if missing:
        result.update(status="unknown", reason="No site benefit is computed while any input, period or extra-harm boundary is unknown")
        return result
    a, u = p["availability"], p["uptake_when_available"]
    share = a[0]*u[0], a[1]*u[1]
    n0, n1 = p["baseline_exposure"], p["intervention_exposure"]
    r0, rb, rf = p["baseline_risk"], p["facility_risk"], p["fallback_risk"]
    mix_lo = min(t*rb[0]+(1-t)*rf[0] for t in share)
    mix_hi = max(t*rb[1]+(1-t)*rf[1] for t in share)
    baseline = n0[0]*r0[0], n0[1]*r0[1]
    d, s = p["deployment_harm"], p["displacement_harm"]
    intervention = n1[0]*mix_lo+d[0]+s[0], n1[1]*mix_hi+d[1]+s[1]
    avoided = baseline[0]-intervention[1], baseline[1]-intervention[0]
    status_order = {"observed": 0, "estimated": 1, "scenario": 2}
    result.update(status="conditional_range", baseline_expected_harm=list(map(str, baseline)),
                  intervention_expected_harm=list(map(str, intervention)), net_harm_avoided=list(map(str, avoided)),
                  facility_use_share=list(map(str, share)), input_evidence_status=max(statuses, key=status_order.get),
                  net_range_sign="positive" if avoided[0] > 0 else "negative" if avoided[1] < 0 else "contains_zero",
                  extra_harm_boundary=case["extra_harm_disjointness_evidence"])
    return result


def _case(design_sha256, site_id, outcome, values=None):
    period = None if values is None else ["2027-01-01", "2028-01-01"]
    case = {"schema_version": 1, "design_sha256": design_sha256, "site_id": site_id, "outcome": outcome,
            "period": period, "source_model_sha256": SOURCE_MODEL["sha256"],
            "extra_harm_disjointness_evidence": None if values is None else "Synthetic control explicitly assigns deployment, displaced-site and main-site events to disjoint categories", "inputs": {}}
    for name, unit in UNITS.items():
        case["inputs"][name] = {"input_id": site_id+":"+outcome+":"+name, "bounds": None if values is None else values[name],
                                "unit": unit, "site_id": site_id, "outcome": outcome, "period": deepcopy(period),
                                "design_sha256": design_sha256,
                                "evidence": {"status": "unknown" if values is None else "scenario", "source": None if values is None else "Authored exact arithmetic verification fixture",
                                             "as_of": None if values is None else "2026-11-24",
                                             "interpretation": "Not established for this site" if values is None else "Synthetic arithmetic only; not a Tillamook estimate"}}
    return case


def synthetic_controls(design_sha256):
    point = lambda x: [x, x]
    common = {"baseline_exposure": point("100000"), "intervention_exposure": point("110000"),
              "availability": point("4/5"), "uptake_when_available": point("3/4")}
    controls = []
    for outcome, scale in (("deaths", 1), ("serious_nonfatal_injuries", 10)):
        values = {**common, "baseline_risk": point(str(exact(scale)/100000)), "facility_risk": point(str(exact(scale)/1000000)),
                  "fallback_risk": point(str(exact(scale)/100000)), "deployment_harm": point(str(exact(scale)/100)),
                  "displacement_harm": point(str(exact(scale)/50))}
        case = _case(design_sha256, "synthetic-arithmetic-control", outcome, values)
        controls.append({"id": "point_" + outcome, "input": case, "result": solve_outcome(case)})
    interval_values = {"baseline_exposure": ["100", "200"], "intervention_exposure": ["100", "300"],
                       "availability": ["1/2", "1"], "uptake_when_available": ["1/2", "1"],
                       "baseline_risk": ["0", "1/1000"], "facility_risk": ["1/1000000", "1/100000"],
                       "fallback_risk": ["1/10000", "1/1000"], "deployment_harm": ["0", "1/50"],
                       "displacement_harm": ["0", "3/100"]}
    case = _case(design_sha256, "synthetic-arithmetic-control", "deaths", interval_values)
    controls.append({"id": "interval_crosses_zero", "input": case, "result": solve_outcome(case)})
    return controls


def build_outcome_report(design_sha256, site_report):
    digest(design_sha256, "design_sha256")
    if type(site_report.get("schema_version")) is not int or site_report.get("schema_version") != 1 or site_report.get("design_sha256") != design_sha256 or site_report.get("physical_approval") is not False:
        raise ValueError("Site report must match the declared design and research-only schema")
    sites = site_report["sites"]
    if not isinstance(sites, list) or not sites:
        raise ValueError("Explicit source sites are required")
    seen, rows = set(), []
    for site in sites:
        site_id = _label(site["id"], "site_id")
        if site_id in seen:
            raise ValueError("Duplicate source site ID")
        seen.add(site_id)
        outcomes = []
        for outcome in OUTCOMES:
            case = _case(design_sha256, site_id, outcome)
            outcomes.append({"outcome": outcome, "input": case, "result": solve_outcome(case)})
        rows.append({"site_id": site_id, "label": _label(site["label"], "site label"), "outcomes": outcomes})
    return {"schema_version": 1, "design_sha256": design_sha256, "site_report_sha256": canonical_sha256(site_report),
            "source_model": deepcopy(SOURCE_MODEL), "site_results": rows, "synthetic_controls": synthetic_controls(design_sha256),
            "summary": {"site_count": len(rows), "site_outcome_evaluations": len(rows)*2, "unknown_site_outcomes": len(rows)*2,
                        "computed_site_benefits": 0, "positive_site_benefit_claims": 0, "physical_approval": False},
            "input_units": UNITS, "physical_approval": False,
            "formula": "N0*r0 - N1*(a*u*rb + (1-a*u)*rf) - D - S, separately for each outcome",
            "availability_definition": "Fraction of intervention crossing opportunities for which the facility is available, not an assumed clock-time uptime",
            "uptake_definition": "Fraction choosing the facility conditional on its availability",
            "scope": "No geometry or model-service result supplies exposure, uptake, risk, deployment harm, displacement harm or causal evidence"}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=ROOT/"design.json")
    parser.add_argument("--site", type=Path, default=ROOT/"site-report.json")
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    read = lambda path: json.loads(path.read_text(encoding="utf-8"))
    report = build_outcome_report(canonical_sha256(read(args.design)), read(args.site))
    report["source_file_sha256"] = {"design.json": sha256(args.design.read_bytes()).hexdigest(), "site-report.json": sha256(args.site.read_bytes()).hexdigest()}
    output = json.dumps(report, indent=2, sort_keys=True, ensure_ascii=True, allow_nan=False)+"\n"
    if args.write:
        (ROOT/"outcome-report.json").write_text(output, encoding="utf-8")
    else:
        print(output, end="")


if __name__ == "__main__":
    main()
