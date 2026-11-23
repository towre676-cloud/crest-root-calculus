"""Frozen exact-output comparison with PARI/GP; serial processes and CPU affinity.

Every planned observation includes a cold single solve, then an automatically
batched compute measurement. GP's integer millisecond clock is used only for
batches of at least 100 ms. No timing outcome is discarded or replaced.
"""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
import os
from pathlib import Path
import platform
import random
import resource
import signal
import subprocess
import sys
import time

from model import expand_case
from verify import independent_coefficients

HERE = Path(__file__).resolve().parent
METHODS = ["crest-cpp", "pari-dense", "pari-periodic-linear", "pari-periodic-doubling"]
DEPENDENCIES = ["model.py", "solver.cpp", "pari.gp", "panel.py", "benchmark.py", "verify.py"]


def digest(path):
    return hashlib.sha256(Path(path).read_bytes()).hexdigest()


def output_digest(coefficients):
    values = [str(int(x)) for x in coefficients]
    return hashlib.sha256(json.dumps(values, separators=(",", ":")).encode()).hexdigest()


def write_json(path, value):
    temp = path.with_suffix(path.suffix + ".tmp")
    temp.write_text(json.dumps(value, indent=2) + "\n")
    temp.replace(path)


def cpu_usage():
    u = resource.getrusage(resource.RUSAGE_CHILDREN)
    return u.ru_utime + u.ru_stime


def invocation(case, method, iterations, gp, native, cpu, limit):
    a, b = expand_case(case)
    if method == "crest-cpp":
        command = [str(native)]
        source = f"{len(a)} {iterations}\n" + " ".join(map(str,a)) + "\n" + " ".join(map(str,b)) + "\n"
    else:
        command = [str(gp), "-f", "-q", "-s", "256000000", "--default", "parisizemax=1073741824"]
        source = (f"read({json.dumps(str(HERE/'pari.gp'))});\n"
                  + f"a={a};b={b};method={json.dumps(method)};iterations={iterations};\n"
                  + "gettime();\n"
                  + "for(k=1,iterations,result=cr_solve(a,b,method));\n"
                  + "elapsed=gettime();\n"
                  + 'print("{\\\"period\\\":",result[2],",\\\"iterations\\\":",iterations,",\\\"cpu_milliseconds\\\":",elapsed,",\\\"coefficients\\\":",Vecrev(Vec(result[1])),"}");\n'
                  + "quit;\n")

    def limits():
        os.setsid()
        resource.setrlimit(resource.RLIMIT_AS, (4 * 1024**3, 4 * 1024**3))
        os.sched_setaffinity(0, {cpu})

    before_cpu, before = cpu_usage(), time.perf_counter()
    process = subprocess.Popen(command, stdin=subprocess.PIPE, stdout=subprocess.PIPE,
                               stderr=subprocess.PIPE, text=True, preexec_fn=limits)
    try:
        stdout, stderr = process.communicate(source, timeout=limit)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        stdout, stderr = process.communicate()
        return {"status": "timeout", "iterations": iterations, "limit_seconds": limit,
                "process_wall_seconds": time.perf_counter()-before,
                "process_cpu_seconds": cpu_usage()-before_cpu,
                "stderr": stderr[-3000:]}
    wall, cpu_total = time.perf_counter()-before, cpu_usage()-before_cpu
    try:
        if process.returncode:
            raise RuntimeError(f"Exit status {process.returncode}")
        lines = [line for line in stdout.splitlines() if line.startswith("{")]
        if len(lines) != 1 or "***" in stderr:
            raise RuntimeError("No unique successful solver record")
        raw = json.loads(lines[0])
        values = raw.pop("coefficients")
        if not all(type(v) is int or (type(v) is str and v.lstrip("-").isdigit()) for v in values):
            raise ArithmeticError("Noninteger output")
        values = [int(v) for v in values]
        if len(values) != len(a)+1 or values[-1] != 1:
            raise ArithmeticError("Wrong characteristic polynomial degree or leading coefficient")
        seconds = raw.get("cpu_seconds", raw.get("cpu_milliseconds", 0)/1000)
        return {"status": "complete", **raw, "kernel_cpu_seconds": seconds,
                "process_wall_seconds": wall, "process_cpu_seconds": cpu_total,
                "output_sha256": output_digest(values), "coefficients": values}
    except Exception as error:
        return {"status": "error", "iterations": iterations, "reason": str(error),
                "process_wall_seconds": wall, "process_cpu_seconds": cpu_total,
                "stdout": stdout[-3000:], "stderr": stderr[-3000:]}


def checked_invocation(case, method, iterations, protocol, expected):
    row = invocation(case, method, iterations, Path(protocol["gp_path"]),
                     Path(protocol["native_path"]), protocol["worker_cpu"],
                     protocol["per_process_wall_limit_seconds"])
    if row["status"] == "complete":
        coefficients = row.pop("coefficients")
        if coefficients != expected:
            row["status"] = "wrong_answer"
            row["actual_coefficients"] = [str(x) for x in coefficients]
        a, b = expand_case(case)
        expected_period = len(a) if method == "pari-dense" else next(
            m for m in range(1,len(a)+1) if len(a)%m == 0
            and all(a[i]==a[i%m] and b[i]==b[i%m] for i in range(len(a))))
        if row["period"] != expected_period:
            row["status"] = "wrong_answer"
            row["expected_period"] = expected_period
    return row


def make_protocol(args, panel, reference):
    rng = random.Random(76127)
    schedule = []
    for replicate in range(1,args.replicates+1):
        order = [(case["id"],method) for case in panel["cases"] for method in METHODS]
        rng.shuffle(order)
        schedule.extend({"case": case,"method": method,"replicate": replicate} for case,method in order)
    return {"schema": 1, "created_utc": datetime.now(timezone.utc).isoformat(),
            "panel_sha256": digest(HERE/"panel.json"), "reference_sha256": digest(reference),
            "source_sha256": {name:digest(HERE/name) for name in DEPENDENCIES},
            "gp_path": str(args.gp.resolve()), "gp_sha256": digest(args.gp),
            "pari_version": subprocess.check_output([str(args.gp), "--version-short"],text=True).strip(),
            "pari_build": json.loads((HERE/"pari-build.json").read_text()) if (HERE/"pari-build.json").exists() else None,
            "native_path": str(args.native.resolve()), "native_sha256": digest(args.native),
            "native_build_command": "g++ -O3 -std=c++17 solver.cpp -lgmpxx -lgmp -o solver",
            "compiler": subprocess.check_output(["g++","--version"],text=True).splitlines()[0],
            "python": sys.version, "platform": platform.platform(),
            "cpu_model": next((line.split(":",1)[1].strip() for line in Path('/proc/cpuinfo').read_text().splitlines() if line.startswith('model name')), "unknown"),
            "affinity_available": sorted(os.sched_getaffinity(0)), "worker_cpu": args.cpu,
            "methods": METHODS, "replicates": args.replicates, "schedule_seed": 76127,
            "per_process_wall_limit_seconds": args.limit, "address_space_bytes": 4*1024**3,
            "minimum_batch_cpu_seconds": 0.1, "maximum_batch_iterations": 1048576,
            "timing": "Each row has a cold single solve in a fresh native process, followed by fresh-process batches doubled until at least 100 ms solver CPU. Every solve validates and discovers period anew. Input text parsing and output serialization are excluded from kernel CPU and included in cold process wall/CPU. No output caching. PARI dense matrix construction is included; common input is the expanded diagonal/bond arrays. All calibration batches are retained. GP clock resolution is 1 ms; native CLOCK_PROCESS_CPUTIME_ID uses nanoseconds.",
            "scope": "Complete exact integer characteristic polynomial; synthetic cyclic Jacobi models, not measured materials or molecular energy data",
            "schedule": schedule}


def check_sources(protocol, include_binaries=True):
    if digest(HERE/"panel.json") != protocol["panel_sha256"]:
        raise ArithmeticError("Frozen panel changed")
    for name, value in protocol["source_sha256"].items():
        if digest(HERE/name) != value:
            raise ArithmeticError("Timed source changed: " + name)
    if include_binaries:
        for name in ("gp","native"):
            if digest(protocol[name+"_path"]) != protocol[name+"_sha256"]:
                raise ArithmeticError("Timed executable changed: " + name)


def run(args):
    target = args.output.resolve()
    protocol_path = target.with_name(target.stem+"-protocol.json")
    reference_path = target.with_name(target.stem+"-reference.json")
    if any(path.exists() for path in (target,protocol_path,reference_path)):
        raise RuntimeError("Refusing to replace an existing study; use a new output name.")
    panel = json.loads((HERE/"panel.json").read_text())
    cases = {case["id"]:case for case in panel["cases"]}
    # Independent determinant recurrence, not the timed transfer/Chebyshev algorithm.
    reference = {key:independent_coefficients(case) for key,case in cases.items()}
    write_json(reference_path, {key:[str(x) for x in values] for key,values in reference.items()})
    protocol = make_protocol(args,panel,reference_path)
    write_json(protocol_path,protocol)  # Frozen before the first timed invocation.
    results = {"schema":1, "completed":False, "protocol_sha256":digest(protocol_path),
               "runs":[], "reference_file": reference_path.name}
    write_json(target,results)
    for i, planned in enumerate(protocol["schedule"]):
        case, method = cases[planned["case"]], planned["method"]
        row = dict(planned)
        row["cold"] = checked_invocation(case,method,1,protocol,reference[case["id"]])
        row["batches"] = []
        row["status"] = row["cold"]["status"]
        if row["status"] == "complete":
            seconds = row["cold"]["kernel_cpu_seconds"]
            if seconds >= protocol["minimum_batch_cpu_seconds"]:
                iterations = 1
            else:
                iterations = min(1048576,max(2,math.ceil(0.12/max(seconds,0.0001))))
            while True:
                batch = checked_invocation(case,method,iterations,protocol,reference[case["id"]])
                row["batches"].append(batch)
                if batch["status"] != "complete":
                    row["status"] = batch["status"]
                    break
                if batch["kernel_cpu_seconds"] >= protocol["minimum_batch_cpu_seconds"]:
                    row["kernel_seconds_per_solve"] = batch["kernel_cpu_seconds"]/iterations
                    break
                if iterations >= protocol["maximum_batch_iterations"]:
                    row["status"] = "timer_floor"
                    break
                iterations = min(protocol["maximum_batch_iterations"],2*iterations)
        results["runs"].append(row)
        write_json(target,results)
        print(f"{i+1}/{len(protocol['schedule'])} {case['id']} {method} r{planned['replicate']} {row['status']}",flush=True)
    check_sources(protocol)
    if digest(reference_path) != protocol["reference_sha256"]:
        raise ArithmeticError("Frozen reference changed")
    results["completed"] = True
    results["dependencies_unchanged"] = True
    results["finished_utc"] = datetime.now(timezone.utc).isoformat()
    write_json(target,results)


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp",type=Path,required=True)
    parser.add_argument("--native",type=Path,required=True)
    parser.add_argument("--cpu",type=int,default=min(os.sched_getaffinity(0)))
    parser.add_argument("--limit",type=float,default=15)
    parser.add_argument("--replicates",type=int,default=3)
    parser.add_argument("--output",type=Path,default=HERE/"results.json")
    run(parser.parse_args())
