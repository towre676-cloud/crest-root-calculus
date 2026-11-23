"""Check all planned measurements, frozen inputs, exact outputs and solver identities."""
import argparse
import hashlib
import json
from pathlib import Path
import random
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def output_digest(value):
    return hashlib.sha256(json.dumps(value,sort_keys=True).encode()).hexdigest()


def check_sources(hashes):
    for source,value in hashes.items():
        if digest(HERE.parent/source)!=value:
            raise ArithmeticError("A timed research dependency changed: "+source)


def verify_study(name):
    path = HERE/(name+"-results.json")
    protocol_path = HERE/(name+"-results-protocol.json")
    panel_path = HERE/(name+"-panel.json")
    result = json.loads(path.read_text())
    protocol = json.loads(protocol_path.read_text())
    panel = json.loads(panel_path.read_text())
    if not result["completed"] or result["protocol_sha256"]!=digest(protocol_path):
        raise ArithmeticError("The study is incomplete or its protocol changed.")
    if protocol["panel_sha256"]!=digest(panel_path):
        raise ArithmeticError("The frozen input panel changed.")
    worker = HERE/{"map":"map_worker.py","quartic":"worker.py",
                   "frobenius":"frobenius_benchmark.py"}[name]
    if protocol["worker_sha256"]!=digest(worker):
        raise ArithmeticError("The timed solver adapter changed.")
    check_sources(protocol.get("research_source_sha256",{}))
    if protocol.get("runner_sha256") and protocol["runner_sha256"]!=digest(HERE/"map_benchmark.py"):
        raise ArithmeticError("The map study runner changed.")
    if name=="map":
        reference_path = HERE/protocol["reference_file"]
        if digest(reference_path)!=protocol["reference_sha256"]:
            raise ArithmeticError("The frozen exact reference changed.")
        from map_benchmark import prepare_reference
        expected_outputs = prepare_reference(reference_path,panel_path)
        if expected_outputs!=result["outputs"]:
            raise ArithmeticError("Study answers differ from the verified exact reference.")
    for resumption in result.get("resumptions",[]):
        check_sources(resumption["research_source_sha256"])
        if not resumption.get("dependencies_unchanged"):
            raise ArithmeticError("A resumed study lacks a final dependency check.")
    if name!="quartic" and not result.get("dependencies_unchanged"):
        raise ArithmeticError("The study lacks a final dependency check.")
    cases = {row["id"]:row for row in panel["cases"]}
    expected = {(case,method,replicate) for case in cases for method in protocol["methods"]
                for replicate in range(1,protocol["replicates"]+1)}
    observed = set()
    rng = random.Random({"quartic":76123,"map":76124,"frobenius":76125}[name])
    schedule = []
    for replicate in range(1,protocol["replicates"]+1):
        order = [(case["id"],method) for case in panel["cases"] for method in protocol["methods"]]
        rng.shuffle(order)
        schedule.extend((case,method,replicate) for case,method in order)
    if [(r["case"],r["method"],r["replicate"]) for r in result["runs"]]!=schedule:
        raise ArithmeticError("The measured rows differ from the seeded plan.")
    for row in result["runs"]:
        key = (row["case"],row["method"],row["replicate"])
        if key not in expected or key in observed:
            raise ArithmeticError("Unexpected or duplicated run.")
        observed.add(key)
        if row["status"]=="error":
            raise ArithmeticError("A solver error remains in the receipt.")
        if row["status"]=="timeout":
            if row["limit_seconds"]!=protocol["per_process_wall_limit_seconds"]:
                raise ArithmeticError("A run used a different resource limit.")
            continue
        if row["status"] not in {"complete","wrong_answer"} or min(row["solve_cpu_seconds"],row["solve_wall_seconds"],
                                             row["process_wall_seconds"])<=0:
            raise ArithmeticError("Malformed completion timing.")
        output = result["supports"][row["case"]] if name=="quartic" else result["outputs"][row["case"]]
        if row["status"]=="wrong_answer":
            if name!="map" or row["expected_output"]!=output or row["actual_output"]==output:
                raise ArithmeticError("Incorrectly classified wrong answer.")
            output = row["actual_output"]
        if name=="frobenius":
            if row["result_sha256"]!=output_digest(output):
                raise ArithmeticError("A returned point-count output hash disagrees.")
        else:
            basis = output["rational_maps"] if name=="map" else output["basis"]
            count = output["geometric_maps"] if name=="map" else output["geometric_marks"]
            if row["support_sha256"]!=output_digest(basis) or row["geometric_marks"]!=count:
                raise ArithmeticError("A returned support hash or geometric count disagrees.")
    if observed!=expected:
        raise ArithmeticError("The receipt omits planned runs.")
    if name=="map":
        x = sp.Symbol("x")
        for case in cases.values():
            output = result["outputs"][case["id"]]
            F = sum(sp.Rational(a)*x**j for j,a in enumerate(case["source"]))
            G = sum(sp.Rational(a)*x**j for j,a in enumerate(case["target"]))
            if len(output["rational_maps"])>output["geometric_maps"]:
                raise ArithmeticError("Rational points exceed the geometric support.")
            for mapping in output["rational_maps"]:
                U = sum(sp.Rational(a)*x**j for j,a in enumerate(mapping["u"]))
                V = sum(sp.Rational(a)*x**j for j,a in enumerate(mapping["v_over_y"]))
                if sp.Poly(U,x).degree()!=case["degree"] or sp.expand(G.subs(x,U)-F*V**2)!=0:
                    raise ArithmeticError("A published map fails its input curve identity.")
    if name=="frobenius":
        if protocol["cohorts"]!=panel["cohorts"] or protocol["field_element_budget"]!=panel["field_element_budget"]:
            raise ArithmeticError("The frozen counting cohorts or table budget changed.")
        members = [case for group in panel["cohorts"].values() for case in group]
        if len(set(members))!=len(members) or set(members)!=set(cases):
            raise ArithmeticError("The counting cohorts do not partition the frozen inputs.")
        for case in cases.values():
            output = result["outputs"][case["id"]]
            numerator,counts = output["zeta_numerator"],output["counts"]
            g,p = case["genus"],case["prime"]
            if len(numerator)!=2*g+1 or numerator[0]!=1 or len(counts)!=g:
                raise ArithmeticError("Wrong local numerator degree or normalization.")
            if any(numerator[2*g-j]!=p**(g-j)*numerator[j] for j in range(g+1)):
                raise ArithmeticError("The local numerator violates reciprocity.")
            traces = [p**(k+1)+1-count for k,count in enumerate(counts)]
            if any(k*numerator[k]+traces[k-1]+sum(numerator[j]*traces[k-j-1]
                   for j in range(1,k))!=0 for k in range(1,g+1)):
                raise ArithmeticError("The extension counts disagree with Newton identities.")
    return {"study":name,"runs":len(observed),"complete":sum(r["status"]=="complete" for r in result["runs"]),
            "timeouts":sum(r["status"]=="timeout" for r in result["runs"]),
            "wrong_answers":sum(r["status"]=="wrong_answer" for r in result["runs"]),"cases":len(cases)}


def verify_aborted_map():
    """Audit preserved pre-correction data without treating it as a completed study."""
    raw = json.loads((HERE/"map-aborted-results.json").read_text())
    protocol_path = HERE/"map-aborted-protocol.json"
    protocol = json.loads(protocol_path.read_text())
    if raw["completed"] or raw["protocol_sha256"]!=digest(protocol_path):
        raise ArithmeticError("The preserved aborted protocol changed.")
    if protocol["worker_sha256"]!=digest(HERE/"map-worker-before-correction.py"):
        raise ArithmeticError("The preserved stock adapter changed.")
    if protocol["panel_sha256"]!=digest(HERE/"map-panel.json"):
        raise ArithmeticError("The restarted study changed its input panel.")
    check_sources(protocol["research_source_sha256"])
    cases = json.loads((HERE/"map-panel.json").read_text())["cases"]
    schedule = [(case["id"],method) for case in cases for method in protocol["methods"]]
    random.Random(76124).shuffle(schedule)
    if [(r["case"],r["method"]) for r in raw["runs"]]!=schedule[:len(raw["runs"])]:
        raise ArithmeticError("The preserved aborted prefix differs from its plan.")
    expected = json.loads((HERE/"map-results.json").read_text())["outputs"]
    retrospective = []
    for row in raw["runs"]:
        if row["status"]!="complete":
            continue
        actual = raw["outputs"][row["case"]]
        if output_digest(actual["rational_maps"])!=row["support_sha256"] or actual["geometric_maps"]!=row["geometric_marks"]:
            raise ArithmeticError("An archived output differs from its original measurement.")
        if actual!=expected[row["case"]]:
            retrospective.append({"case":row["case"],"method":row["method"],"replicate":row["replicate"],
                                  "actual_output":actual,"expected_output":expected[row["case"]]})
    return {"recorded_runs_preserved":len(raw["runs"]),
            "wrong_answers_identified_in_preserved_prefix":retrospective,
            "interruption":"The following g2-d3-twist SymPy run triggered the exact-output assertion before its timing row was serialized. No timing has been invented for that unrecorded row."}


def verify_regression_identities():
    """Replay mathematical witnesses without requiring a future library to retain its old bug."""
    from verify_sympy_regression import coefficient_recurrence
    evidence = json.loads((HERE/"sympy-regression.json").read_text())
    x = sp.Symbol("x")
    poly = lambda values: sp.Poly(sum(sp.Rational(a)*x**i for i,a in enumerate(values)),x,domain=sp.QQ)
    for key in ("minimal_counterexample","benchmark_counterexample"):
        row = evidence[key]
        F,U,G = map(poly,(row["source"],row["exact_inner"],row["exact_outer"]))
        if F.gcd(F.diff()).degree() or sp.expand(G.as_expr().subs(x,U.as_expr())-F.as_expr())!=0:
            raise ArithmeticError("An archived decomposition identity failed.")
        corrected = coefficient_recurrence(F.monic(),U.degree())
        if corrected!=poly(row["corrected_monic_zero_constant_inner"]):
            raise ArithmeticError("The corrected recurrence failed its exact witness.")
    row = evidence["benchmark_counterexample"]
    case = next(c for c in json.loads((HERE/"map-panel.json").read_text())["cases"] if c["id"]==row["id"])
    if poly(case["source"])!=poly(row["source"]) or poly(case["target"])!=poly(row["target"]):
        raise ArithmeticError("The regression witness no longer matches its frozen map input.")
    F,U,G = map(poly,(row["source"],row["exact_inner"],row["target"]))
    if sp.expand(G.as_expr().subs(x,U.as_expr())-sp.Rational(row["vertical_scale_squared"])*F.as_expr())!=0:
        raise ArithmeticError("The geometric-map witness failed.")
    expected = json.loads((HERE/"map-results.json").read_text())["outputs"][row["id"]]
    if (row["geometric_maps"],row["rational_maps"])!=(expected["geometric_maps"],len(expected["rational_maps"])):
        raise ArithmeticError("The regression map counts differ from the verified reference.")
    return {"exact_identity_witnesses":2,"archived_stock_version":evidence["sympy_version"],
            "current_version_required_to_retain_bug":False}


def verify_preliminary_counts():
    raw = json.loads((HERE/"frobenius-preliminary-results.json").read_text())
    path = HERE/"frobenius-preliminary-protocol.json"
    protocol = json.loads(path.read_text())
    old_panel = json.loads((HERE/"frobenius-preliminary-panel.json").read_text())
    panel = json.loads((HERE/"frobenius-panel.json").read_text())
    result = json.loads((HERE/"frobenius-results.json").read_text())
    if not raw["completed"] or raw["protocol_sha256"]!=digest(path):
        raise ArithmeticError("The preliminary counting study was not preserved.")
    if protocol["panel_sha256"]!=digest(HERE/"frobenius-preliminary-panel.json") or protocol["worker_sha256"]!=digest(HERE/"frobenius-preliminary-worker.py"):
        raise ArithmeticError("The preliminary counting inputs or adapter changed.")
    check_sources(protocol["research_source_sha256"])
    if panel["cases"][:len(old_panel["cases"])]!=old_panel["cases"]:
        raise ArithmeticError("An original counting input changed during extension.")
    if any(result["outputs"].get(case)!=value for case,value in raw["outputs"].items()):
        raise ArithmeticError("The final counting study changed an earlier exact polynomial.")
    expected = {(case["id"],method,rep) for case in old_panel["cases"] for method in protocol["methods"] for rep in (1,2)}
    observed = set()
    for row in raw["runs"]:
        key = (row["case"],row["method"],row["replicate"])
        if key in observed or key not in expected or row["status"]!="complete":
            raise ArithmeticError("Unexpected archived counting measurement.")
        observed.add(key)
        if row["result_sha256"]!=output_digest(raw["outputs"][row["case"]]):
            raise ArithmeticError("An archived counting output changed.")
    if observed!=expected:
        raise ArithmeticError("The preliminary counting receipt lost measurements.")
    return {"original_inputs_preserved":len(old_panel["cases"]),"original_runs_preserved":len(raw["runs"]),
            "exact_polynomials_unchanged":True}


def recompute_counts():
    sys.path.insert(0,str(HERE.parent/"bounded-correspondences"))
    from frobenius_filter import curve,local_data
    counts = json.loads((HERE/"frobenius-results.json").read_text())["outputs"]
    protocol = json.loads((HERE/"frobenius-results-protocol.json").read_text())
    for case in json.loads((HERE/"frobenius-panel.json").read_text())["cases"]:
        output = local_data(curve(case["source"]),case["prime"],protocol["field_element_budget"])
        if output!=counts[case["id"]]:
            raise ArithmeticError("A local polynomial failed exact finite-field recounting.")


def recompute():
    sys.path.insert(0,str(HERE.parent/"chart-one-cofactors"))
    from solve import describe
    from proof import verify,crt_targets
    from worker import canonical_support
    quartic = json.loads((HERE/"quartic-results.json").read_text())["supports"]
    for case in json.loads((HERE/"quartic-panel.json").read_text())["cases"]:
        report = describe(case["quartic"],case["mark"],details=True)
        verify(case["quartic"],case["mark"],report,raw=True)
        basis = canonical_support([p.as_expr() for label,p in crt_targets(report["fibres"])])
        if basis!=quartic[case["id"]]["basis"] or report["geometric_marks"]!=quartic[case["id"]]["geometric_marks"]:
            raise ArithmeticError("A quartic output failed independent cofactor replay.")
    sys.path.insert(0,str(HERE.parent/"bounded-correspondences"))
    from composition import recover_constant,verify_constant
    from map_worker import solve_case as map_solve
    maps = json.loads((HERE/"map-results.json").read_text())["outputs"]
    for case in json.loads((HERE/"map-panel.json").read_text())["cases"]:
        report = recover_constant(case["source"],case["target"])
        verify_constant(report)
        ordered = sorted(report["rational_maps"],key=lambda value:json.dumps(value,sort_keys=True))
        if ordered!=maps[case["id"]]["rational_maps"] or report["geometric_maps"]!=maps[case["id"]]["geometric_maps"]:
            raise ArithmeticError("A map output failed exact reconstruction replay.")
        corrected = map_solve(case,"corrected-dense-decomposition")
        if corrected!=maps[case["id"]]:
            raise ArithmeticError("The independently corrected dense method disagrees with the exact reference.")
    recompute_counts()


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute",action="store_true")
    parser.add_argument("--recompute-counts",action="store_true",help="Recount the counting panel without repeating quartic and map reconstruction.")
    args = parser.parse_args()
    results = [verify_study(name) for name in ("quartic","map","frobenius")]
    if args.recompute:
        recompute()
    elif args.recompute_counts:
        recompute_counts()
    print(json.dumps({"studies":results,"aborted_map_audit":verify_aborted_map(),"regression_identities":verify_regression_identities(),"preliminary_counts":verify_preliminary_counts(),
                      "recomputed_scope":"all" if args.recompute else "counting" if args.recompute_counts else "certificates"},indent=2))
