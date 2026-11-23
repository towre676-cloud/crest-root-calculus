"""Audit the entire frozen comparison and independently reconstruct every answer."""
import argparse
from collections import Counter
import copy
import hashlib
import json
import math
from pathlib import Path
import random

from model import expand_case
from verify import independent_coefficients, verify_controls
from spectral_curve import verify_all

HERE = Path(__file__).resolve().parent


def digest(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def output_digest(values):
    return hashlib.sha256(json.dumps([str(int(x)) for x in values], separators=(",", ":")).encode()).hexdigest()


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def audit_data(raw, protocol, reference, panel):
    require(raw["completed"] is True and raw["dependencies_unchanged"] is True, "Study was not completed with a dependency check")
    cases = {case["id"]:case for case in panel["cases"]}
    require(len(cases)==len(panel["cases"]),"Duplicate input identifiers")
    require(set(reference)==set(cases),"Reference does not cover exactly the panel")
    expected_schedule=[]
    rng=random.Random(protocol["schedule_seed"])
    for replicate in range(1,protocol["replicates"]+1):
        order=[(case["id"],method) for case in panel["cases"] for method in protocol["methods"]]
        rng.shuffle(order)
        expected_schedule.extend({"case":case,"method":method,"replicate":replicate} for case,method in order)
    require(protocol["schedule"]==expected_schedule,"Protocol schedule differs from the seeded design")
    require([{key:row[key] for key in ("case","method","replicate")} for row in raw["runs"]]==expected_schedule,
            "Measured rows differ from the full seeded schedule")
    periods={}
    for key,case in cases.items():
        a,b=expand_case(case)
        periods[key]=next(m for m in range(1,len(a)+1) if len(a)%m==0
                         and all(a[i]==a[i%m] and b[i]==b[i%m] for i in range(len(a))))
    observations=0
    for row in raw["runs"]:
        require(row["status"] in {"complete","timeout","error","wrong_answer","timer_floor"},
                "Unknown row status")
        require(row["cold"]["iterations"]==1,"Cold observation was not a single solve")
        if row["cold"]["status"]!="complete":
            require(row["status"]==row["cold"]["status"] and not row["batches"],"Malformed cold failure")
            require("kernel_seconds_per_solve" not in row,"Failed cold solve has a successful normalized timing")
        for measure in [row["cold"],*row["batches"]]:
            observations+=1
            require(measure["iterations"]>=1,"Nonpositive batch size")
            for key in ("process_wall_seconds","process_cpu_seconds"):
                require(math.isfinite(measure[key]) and measure[key]>0,"Malformed process timing")
            if measure["status"]=="timeout":
                require(measure["limit_seconds"]==protocol["per_process_wall_limit_seconds"],"Wrong timeout limit")
                continue
            if measure["status"]=="error":
                continue
            if measure["status"]=="wrong_answer":
                require("actual_coefficients" in measure or "expected_period" in measure,"Missing incorrect output")
                continue
            require(measure["status"]=="complete","Unknown observation status")
            require(measure["output_sha256"]==output_digest(reference[row["case"]]),"Returned output hash differs from exact reference")
            a,_=expand_case(cases[row["case"]])
            expected_period=len(a) if row["method"]=="pari-dense" else periods[row["case"]]
            require(measure["period"]==expected_period,"Incorrect discovered-period metadata")
            require(math.isfinite(measure["kernel_cpu_seconds"]) and measure["kernel_cpu_seconds"]>=0,"Malformed kernel timing")
            clock=measure["cpu_seconds"] if "cpu_seconds" in measure else measure["cpu_milliseconds"]/1000
            require(measure["kernel_cpu_seconds"]==clock,"Recorded clock changed")
        if row["cold"]["status"]=="complete":
            require(row["batches"],"Successful cold solve has no calibration observation")
            for previous in row["batches"][:-1]:
                require(previous["status"]=="complete","Measurement continued after a failed batch")
                require(previous["kernel_cpu_seconds"]<protocol["minimum_batch_cpu_seconds"],
                        "Unnecessary batch replaced a sufficient measurement")
            final=row["batches"][-1]
            if final["status"]!="complete":
                require(row["status"]==final["status"],"Top-level status differs from the failed final batch")
                require("kernel_seconds_per_solve" not in row,"Failed final batch has a successful normalized timing")
            elif final["kernel_cpu_seconds"]>=protocol["minimum_batch_cpu_seconds"]:
                require(row["status"]=="complete","Complete observations were reclassified as a failure")
            else:
                require(row["status"]=="timer_floor"
                        and final["iterations"]==protocol["maximum_batch_iterations"],
                        "Insufficient final batch is not a genuine timer-floor outcome")
                require("kernel_seconds_per_solve" not in row,"Timer-floor outcome has a successful normalized timing")
        if row["status"]=="complete":
            require(row["batches"] and all(x["status"]=="complete" for x in row["batches"]),"Successful row contains a failed batch")
            final=row["batches"][-1]
            require(final["kernel_cpu_seconds"]>=protocol["minimum_batch_cpu_seconds"],"Batch is below the clock-resolution floor")
            require(row["kernel_seconds_per_solve"]==final["kernel_cpu_seconds"]/final["iterations"],"Normalized CPU timing changed")
    return {"inputs":len(cases),"planned_rows":len(expected_schedule),"cold_and_batch_observations":observations,
            "statuses":dict(Counter(row["status"] for row in raw["runs"])),
            "actual_periods":periods}


def audit(path=HERE/"results.json", controls=True):
    protocol_path=path.with_name(path.stem+"-protocol.json")
    raw=json.loads(path.read_text())
    protocol=json.loads(protocol_path.read_text())
    panel=json.loads((HERE/"panel.json").read_text())
    refpath=path.parent/raw["reference_file"]
    reference=json.loads(refpath.read_text())
    require(raw["protocol_sha256"]==digest(protocol_path),"Protocol hash changed")
    require(protocol["panel_sha256"]==digest(HERE/"panel.json"),"Panel hash changed")
    require(protocol["reference_sha256"]==digest(refpath),"Reference hash changed")
    for name,sha in protocol["source_sha256"].items():
        require(digest(HERE/name)==sha,"Timed dependency changed: "+name)
    for case in panel["cases"]:
        require([int(x) for x in reference[case["id"]]]==independent_coefficients(case),"Reference failed full-ring independent determinant replay")
    summary=audit_data(raw,protocol,reference,panel)
    summary["all_coefficients_independently_recomputed"]=True
    if controls:
        summary["model_controls"]=verify_controls()
        spectral=verify_all()
        summary["spectral_curve_controls"]={"universal_identity_checks":spectral["universal_identity_checks"],
                                             "examples":len(spectral["examples"]),
                                             "checks_passed":spectral["checks_passed"]}
        tampered=[]
        changed=copy.deepcopy(raw);changed["runs"].pop();tampered.append(changed)
        changed=copy.deepcopy(raw);changed["runs"][0]["cold"]["output_sha256"]="0"*64;tampered.append(changed)
        successful=next(i for i,row in enumerate(raw["runs"]) if row["status"]=="complete")
        changed=copy.deepcopy(raw);changed["runs"][successful]["kernel_seconds_per_solve"]*=2;tampered.append(changed)
        changed=copy.deepcopy(raw);changed["runs"][successful]["cold"]["period"]+=1;tampered.append(changed)
        changed=copy.deepcopy(raw);changed["runs"][successful]["status"]="error";tampered.append(changed)
        changed=copy.deepcopy(raw);changed["runs"][successful]["cold"]["status"]="error";tampered.append(changed)
        for changed in tampered:
            try: audit_data(changed,protocol,reference,panel)
            except ArithmeticError: pass
            else: raise ArithmeticError("Tampered benchmark receipt was accepted")
        summary["altered_study_receipts_rejected"]=len(tampered)
    return summary


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--results",type=Path,default=HERE/"results.json")
    args=parser.parse_args()
    print(json.dumps(audit(args.results),indent=2))
