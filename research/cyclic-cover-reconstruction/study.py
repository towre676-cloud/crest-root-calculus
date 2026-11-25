"""Freeze, measure and replay the exact cyclic-cover comparison."""

import argparse
from collections import Counter
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import secrets
import shutil
import signal
import statistics
import subprocess
import sys
import time

import sympy as sp

from cyclic import HERE, polynomial, solve_stratum
from finite import rational_points
from fixtures import build_panel
from raw import raw_solve


METHODS = ("derivative", "native-formal-root", "native-original-identity")
REPETITIONS = 2
LIMIT = 4
ROOT = HERE.parent.parent
SOURCE_FILES = [
    HERE / name for name in ("cyclic.py", "fixtures.py", "raw.py", "study.py")
] + [
    HERE.parent / "coefficient-quotients" / name
    for name in ("discover.py", "native_control.py", "benchmark.py", "baselines.py")
] + [HERE.parent / "bounded-correspondences" / "finite.py", HERE.parent / "solver-comparison" / "map_worker.py",
     HERE.parent / "solver-comparison" / "worker.py"]
DEFAULT = HERE / "study"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def encode(value):
    return json.dumps(value, sort_keys=True, separators=(",", ":")).encode("utf-8")


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def source_hashes():
    return {path.relative_to(ROOT).as_posix(): digest(path.read_bytes()) for path in SOURCE_FILES}


def output(row):
    return {name: row[name] for name in ("geometric_points", "maps")}


def schedule(panel, seed):
    rows = [(name, method, repetition) for name in panel for method in METHODS
            for repetition in range(1, REPETITIONS+1)]
    random.Random(seed).shuffle(rows)
    return rows


def freeze(directory):
    directory.mkdir(parents=True, exist_ok=True)
    target = directory / "protocol.json"
    if target.exists():
        raise FileExistsError(target)
    hashes = source_hashes()
    write(target, {"source_sha256": hashes, "panel_seed": secrets.randbits(64),
                   "shuffle_seed": secrets.randbits(64), "limit_seconds": LIMIT, "repetitions": REPETITIONS})
    inputs, provenance = build_panel(read(target)["panel_seed"])
    write(directory / "inputs.json", inputs)
    write(directory / "provenance.json", provenance)


def worker(payload):
    # Imports used only by the timed Linux worker. Offline audit works on Windows.
    from unittest.mock import patch
    from benchmark import cpu
    from native_control import native_basis
    os.environ.pop("CREST_SINGULAR", None)
    start = cpu()
    case = payload["input"]
    F, k = polynomial(case["coefficients"]), case["exponent"]
    method = payload["method"]

    def backend(equations, variables):
        dimension, basis = native_basis(payload["executable"], equations, variables)
        with patch("finite.shutil.which", return_value=None):
            certificate = rational_points(basis, variables)
        if certificate["algebra_dimension"] != dimension:
            raise ArithmeticError("Native and local quotient dimensions disagree.")
        return certificate

    if method == "native-original-identity":
        answer = raw_solve(F, k, *case["stratum"], basis_solver=backend)
    else:
        answer = output(solve_stratum(F, k, *case["stratum"], derivative=method == "derivative",
                                      basis_solver=backend))
    return {"solver_cpu": cpu()-start, "output": answer}


def measure(payload):
    from benchmark import cpu
    start = time.perf_counter()
    cpu_start = cpu()
    process = subprocess.Popen([sys.executable, "-B", str(__file__), "worker"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, start_new_session=True)
    try:
        stdout, stderr = process.communicate(encode(payload).decode("utf-8"), timeout=LIMIT)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        return {"status": "timeout", "wall": time.perf_counter()-start, "process_cpu": cpu()-cpu_start}
    observations = {"wall": time.perf_counter()-start, "process_cpu": cpu()-cpu_start}
    if process.returncode:
        return {"status": "process_error", **observations, "error": stderr[-3000:]}
    return {"status": "completed", **observations, **json.loads(stdout)}


def summarize(rows, inputs):
    summary = {}
    for method in METHODS:
        observed = [r for r in rows if r["method"] == method]
        comparisons = {}
        for metric in ("solver_cpu", "process_cpu", "wall"):
            ratios = []
            for name in inputs:
                samples = [[r[metric] for r in rows if r["case"] == name and r["method"] == route
                            and r["status"] == "correct"] for route in ("derivative", method)]
                if all(len(sample) == REPETITIONS for sample in samples):
                    ratios.append(statistics.median(samples[1])/statistics.median(samples[0]))
            comparisons[metric] = {"paired_inputs": len(ratios), "derivative_wins": sum(v > 1 for v in ratios),
                                   "comparator_wins": sum(v < 1 for v in ratios),
                                   "median_comparator_over_derivative": round(statistics.median(ratios), 6) if ratios else None}
        summary[method] = {"rows": dict(Counter(r["status"] for r in observed)),
                           "inputs_completed_twice": sum(sum(r["case"] == name and r["status"] == "correct"
                                                             for r in observed) == REPETITIONS for name in inputs),
                           "paired": comparisons}
    return summary


def run(directory):
    protocol = read(directory / "protocol.json")
    if protocol["source_sha256"] != source_hashes():
        raise ArithmeticError("Source differs from the pre-panel freeze.")
    inputs = read(directory / "inputs.json")
    provenance = read(directory / "provenance.json")
    if (inputs, provenance) != build_panel(protocol["panel_seed"]):
        raise ArithmeticError("Inputs differ from the frozen panel.")
    reference = {name: output(solve_stratum(polynomial(case["coefficients"]), case["exponent"], *case["stratum"]))
                 for name, case in inputs.items()}
    write(directory / "reference.json", reference)
    data_hashes = {name: digest((directory / name).read_bytes()) for name in
                   ("protocol.json", "inputs.json", "provenance.json", "reference.json")}
    executable = os.environ.get("CREST_QUOTIENT_SINGULAR") or shutil.which("Singular")
    if not executable:
        raise RuntimeError("Install Singular or set CREST_QUOTIENT_SINGULAR to its executable.")
    version = subprocess.run([executable, "--version"], check=True, capture_output=True, text=True).stdout.splitlines()[0]
    selected = sorted(os.sched_getaffinity(0))[-1]
    os.sched_setaffinity(0, {selected})
    rows = []
    plan = schedule(inputs, protocol["shuffle_seed"])
    for i, (case, method, repetition) in enumerate(plan, 1):
        result = measure({"input": inputs[case], "method": method, "executable": executable})
        if result["status"] == "completed":
            result["status"] = "correct" if result["output"] == reference[case] else "wrong_answer"
            result["output_sha256"] = digest(encode(result["output"]))
            if result["status"] == "correct":
                del result["output"]
        rows.append({"case": case, "method": method, "repetition": repetition, **result})
        if i % 20 == 0:
            print(f"Measured {i}/{len(plan)} scheduled attempts.", flush=True)
    if (protocol["source_sha256"] != source_hashes()
            or any(digest((directory / name).read_bytes()) != value for name, value in data_hashes.items())):
        raise ArithmeticError("A frozen file changed during timing.")
    receipt = {"data_sha256": data_hashes, "singular_version": version,
               "singular_sha256": digest(Path(executable).read_bytes()), "selected_cpu": selected,
               "host": platform.platform(), "python": platform.python_version(), "sympy": sp.__version__,
               "rows": rows, "summary": summarize(rows, inputs)}
    write(directory / "measurements.json", receipt)


def audit(directory=DEFAULT):
    protocol, inputs, provenance, reference, receipt = (read(directory / f"{name}.json") for name in
                                                        ("protocol", "inputs", "provenance", "reference", "measurements"))
    if (protocol["source_sha256"] != source_hashes() or (inputs, provenance) != build_panel(protocol["panel_seed"])
            or protocol["limit_seconds"] != LIMIT or protocol["repetitions"] != REPETITIONS):
        raise ArithmeticError("The source or protocol does not match the frozen panel.")
    if (set(receipt["data_sha256"]) != {"protocol.json", "inputs.json", "provenance.json", "reference.json"}
            or any(digest((directory / name).read_bytes()) != value for name, value in receipt["data_sha256"].items())):
        raise ArithmeticError("One of the timed input files changed.")
    if [(row["case"], row["method"], row["repetition"]) for row in receipt["rows"]] != schedule(inputs, protocol["shuffle_seed"]):
        raise ArithmeticError("The archived observations do not match the frozen schedule.")
    for name, case in inputs.items():
        exact = output(solve_stratum(polynomial(case["coefficients"]), case["exponent"], *case["stratum"]))
        if exact != reference[name]:
            raise ArithmeticError(f"The independent replay differs on {name}.")
        planted = provenance[name]["planting"]
        if planted is not None and planted not in exact["maps"]:
            raise ArithmeticError("A planted map is missing.")
    for row in receipt["rows"]:
        status = row["status"]
        if status not in {"correct", "wrong_answer", "timeout", "process_error"} or row["wall"] <= 0 or row["process_cpu"] < 0:
            raise ArithmeticError("An observation has an invalid status or duration.")
        if status == "correct" and row["output_sha256"] != digest(encode(reference[row["case"]])):
            raise ArithmeticError("An exact output failed its recorded hash.")
        if status == "wrong_answer" and (row["output"] == reference[row["case"]]
                                        or row["output_sha256"] != digest(encode(row["output"]))):
            raise ArithmeticError("An incorrect output was misclassified.")
    if receipt["summary"] != summarize(receipt["rows"], inputs):
        raise ArithmeticError("The timing summary changed.")
    return {"inputs": len(inputs), "observations": len(receipt["rows"]),
            "geometric_solutions": sum(row["geometric_points"] for row in reference.values()),
            "summary": receipt["summary"]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("action", choices=("freeze", "run", "audit", "worker"))
    parser.add_argument("--directory", type=Path, default=DEFAULT)
    args = parser.parse_args()
    if args.action == "worker":
        print(json.dumps(worker(json.load(sys.stdin))))
    elif args.action == "freeze":
        freeze(args.directory)
    elif args.action == "run":
        run(args.directory)
    else:
        print(json.dumps(audit(args.directory), indent=2))


if __name__ == "__main__":
    main()
