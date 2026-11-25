"""Measure an isolated, frozen coefficient panel, retaining every scheduled attempt."""

import argparse
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import resource
import signal
import subprocess
import sys
import time

import sympy as sp

from baselines import dense, gp, singular
from discover import (coefficients, constant_maps, genus, nonconstant_degrees,
                      normalization, polynomial, variable_maps)

HERE = Path(__file__).resolve().parent
METHODS = ("crest", "crest-no-sieve", "corrected-dense", "singular", "gp", "stock-sympy")
VARIABLE_METHODS = ("derivative", "without-derivative")
REPETITIONS = 3
TIMEOUTS = {"constant": 12, "nonconstant": 4}


def cpu():
    return sum(usage.ru_utime + usage.ru_stime for usage in
               (resource.getrusage(resource.RUSAGE_SELF), resource.getrusage(resource.RUSAGE_CHILDREN)))


def encode(data):
    return json.dumps(data, sort_keys=True, separators=(",", ":")).encode("utf-8")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def read(path):
    return json.loads(path.read_text(encoding="utf-8"))


def write(path, data):
    path.write_text(json.dumps(data, indent=2) + "\n", encoding="utf-8")


def map_list(rows):
    return [{key: row[key] for key in ("degree", "target_degree", "inner", "outer")}
            for row in rows if "outer" in row]


def assert_identities(result, values):
    f, a, p = map(polynomial, (result["source"], result["square"], values))
    if a.LC() != 1 or f * a**2 != p or f.gcd(f.diff()).degree() != 0:
        raise ArithmeticError("Incorrect squarefree normalization.")
    seen = set()
    for mapping in result["maps"]:
        u, g = map(polynomial, (mapping["inner"], mapping["outer"]))
        d, m = mapping["degree"], mapping["target_degree"]
        if (d in seen or u.degree() != d or g.degree() != m or u.LC() != 1 or u.nth(0) != 0
                or d < 2 or not 1 <= genus(g) < genus(f) or g.compose(u) != f
                or g.gcd(g.diff()).degree() != 0):
            raise ArithmeticError("Incorrect constant-multiplier quotient.")
        seen.add(d)


def worker(payload):
    values, method = payload["coefficients"], payload["method"]
    if method in VARIABLE_METHODS:
        started = cpu()
        source = normalization(polynomial(values))[1]
        start = cpu()
        result = variable_maps(source, *payload["stratum"], derivative=method == "derivative")
        elapsed = cpu() - start
        output = {k: result[k] for k in ("geometric_points", "maps")}
        return {"output": output, "phases": {"normalize_cpu": start - started,
                                               "search_and_certificate_cpu": elapsed}}

    started = cpu()
    if method == "gp":
        square, source, maps, native = gp(values)
        phases = {"native_and_adapter_cpu": cpu() - started,
                  "gp_process_normalize_cpu_ms_resolution": native[0],
                  "gp_process_degree_search_cpu_ms_resolution": native[1]}
    else:
        a, f = normalization(polynomial(values))
        square, source = coefficients(a), coefficients(f)
        after_normalization = cpu()
        if method.startswith("crest"):
            maps = map_list(constant_maps(f, sieve=method == "crest"))
        elif method == "singular":
            maps = singular(f)
        else:
            maps = dense(f, corrected=method == "corrected-dense")
        phases = {"normalize_cpu": after_normalization - started, "degree_search_cpu": cpu() - after_normalization}
    result = {"square": square, "source": source, "maps": maps}
    before_check = cpu()
    assert_identities(result, values)
    phases["identity_check_cpu"] = cpu() - before_check
    return {"output": result, "phases": phases}


def constant_reference(values):
    a, f = normalization(polynomial(values))
    return {"square": coefficients(a), "source": coefficients(f), "maps": dense(f)}


def expected_and_schedule(inputs, provenance, shuffle_seed):
    rng = random.Random(shuffle_seed)
    schedule, reference = [], {}
    for identifier, case in inputs.items():
        info = provenance["cases"][identifier]
        if info["study"] == "constant":
            reference[identifier] = constant_reference(case["coefficients"])
            methods, stratum = METHODS, None
        else:
            stratum = info["degree_stratum"]
            source = normalization(polynomial(case["coefficients"]))[1]
            if tuple(stratum) not in nonconstant_degrees(source):
                raise ArithmeticError("Invalid predeclared cover degree stratum.")
            methods = VARIABLE_METHODS
        for method in methods:
            for repetition in range(1, REPETITIONS + 1):
                schedule.append((identifier, info["study"], method, repetition, stratum))
    rng.shuffle(schedule)
    return schedule, reference


def invoke(payload, timeout):
    start = time.perf_counter()
    started_cpu = cpu()
    process = subprocess.Popen([sys.executable, "-B", str(__file__), "--worker"],
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE,
                               text=True, start_new_session=True)
    try:
        stdout, stderr = process.communicate(encode(payload).decode("utf-8"), timeout=timeout)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        return {"status": "timeout", "wall": time.perf_counter() - start, "process_cpu": cpu() - started_cpu}
    row = {"wall": time.perf_counter() - start, "process_cpu": cpu() - started_cpu}
    if process.returncode:
        return {**row, "status": "process_error", "error": stderr[-2400:]}
    return {**row, "status": "completed", **json.loads(stdout)}


def study(directory, freeze, protocol):
    inputs, provenance = (read(directory / name) for name in ("inputs.json", "provenance.json"))
    schedule, reference = expected_and_schedule(inputs, provenance, protocol["shuffle_seed"])
    reference_path = directory / "constant-reference.json"
    write(reference_path, reference)
    protocol["file_sha256"][reference_path.name] = digest(reference_path.read_bytes())
    rows = []
    for index, (case, family, method, repetition, stratum) in enumerate(schedule, 1):
        outcome = invoke({**inputs[case], "method": method, "stratum": stratum}, TIMEOUTS[family])
        if outcome["status"] == "completed" and family == "constant" and outcome["output"] != reference[case]:
            outcome["status"] = "wrong_answer"
        if "output" in outcome:
            outcome["output_sha256"] = digest(encode(outcome["output"]))
        if outcome["status"] == "completed" and family == "constant":
            del outcome["output"]
        rows.append({"case": case, "study": family, "method": method, "repetition": repetition,
                     "stratum": stratum, **outcome})
        if index % 50 == 0:
            print(f"Measured {index}/{len(schedule)} scheduled attempts.", flush=True)
    for relative, hash_value in freeze["sha256"].items():
        if digest((HERE.parent.parent / relative).read_bytes()) != hash_value:
            raise ArithmeticError(f"Frozen source changed during timing: {relative}")
    write(directory / "measurements.json", {"protocol": protocol, "rows": rows})


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--worker", action="store_true")
    parser.add_argument("--directory", type=Path, default=HERE / "study")
    cli = parser.parse_args()
    if cli.worker:
        print(json.dumps(worker(json.load(sys.stdin))))
        return
    freeze = read(cli.directory / "freeze.json")
    for relative, hash_value in freeze["sha256"].items():
        if digest((HERE.parent.parent / relative).read_bytes()) != hash_value:
            raise ArithmeticError(f"Frozen source changed before timing: {relative}")
    cpus = sorted(os.sched_getaffinity(0))
    selected = cpus[-1]
    os.sched_setaffinity(0, {selected})
    protocol = {"shuffle_seed": freeze["schedule_seed"], "selected_cpu": selected, "limits_seconds": TIMEOUTS,
                "repetitions": REPETITIONS, "host": platform.platform(), "python": platform.python_version(),
                "sympy": sp.__version__,
                "file_sha256": {n: digest((cli.directory / n).read_bytes())
                                for n in ("freeze.json", "inputs.json", "provenance.json")},
                "executables": {}}
    install = os.environ.get("CREST_DECOMP_LIB")
    if install:
        protocol["decomp_library_sha256"] = digest(Path(install).read_bytes())
    import shutil
    for label, env, program in [("gp", "CREST_QUOTIENT_GP", "gp"),
                                 ("singular", "CREST_QUOTIENT_SINGULAR", "Singular")]:
        executable = os.environ.get(env) or shutil.which(program)
        if not executable:
            raise RuntimeError(f"The frozen comparison requires {program}.")
        completed = subprocess.run([executable, "--version"], capture_output=True, text=True, check=True)
        protocol["executables"][label] = {"sha256": digest(Path(executable).read_bytes()),
            "version": (completed.stdout or completed.stderr).splitlines()[0].strip()}
    study(cli.directory, freeze, protocol)


if __name__ == "__main__":
    main()
