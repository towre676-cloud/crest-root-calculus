"""Post-study native Singular control on the same twenty nonconstant strata."""

import argparse
import json
import os
from pathlib import Path
import random
import shutil
import signal
import statistics
import subprocess
import sys
import time

import sympy as sp

from benchmark import HERE, cpu, digest, encode, read, write
from discover import X, coefficients, left_divide, normalization, polynomial, variable_incidence, variable_maps
from finite import rational_points


METHODS = ("derivative", "singular-formal-root", "singular-original-identity")
LIMIT = 4
REPEATS = 2
SCHEDULE_SEED = 775144
PATH = HERE / "study" / "native-control.json"


def original_incidence(f, d, m, e):
    horizontal, vertical, outer = (sp.symbols(f"{prefix}0:{length}")
                                    for prefix, length in (("a", d-1), ("b", e), ("c", m)))
    domain = sp.QQ.poly_ring(*(horizontal + vertical + outer))
    u = sp.Poly.from_list([1, *horizontal, 0], X, domain=domain)
    v = sp.Poly.from_list([1, *vertical], X, domain=domain)
    g = sp.Poly.from_list([f.LC(), *outer], X, domain=domain)
    return horizontal + vertical + outer, (g.compose(u) - f.set_domain(domain)*v**2).coeffs(), u, v


def native_basis(executable, equations, variables):
    rendered = []
    for equation in equations:
        terms = sp.Poly(equation, *variables, domain=sp.QQ).terms()
        rendered.append("+".join(f"({c})" + "".join(f"*{v}^{j}" for v, j in zip(variables, exponents)
                                                     if j)
                                 for exponents, c in terms))
    names = ",".join(map(str, variables))
    script = (f"ring r=0,({names}),dp;\nshort=0; option(redSB); ideal I={','.join(rendered)};\n"
              'ideal G=slimgb(I); int dimension=vdim(G); print("DIM "+string(dimension));\n'
              'int j; for(j=1;j<=size(G);j++){print("B "+string(G[j]));}\nquit;\n')
    run = subprocess.run([executable, "-q", "--no-rc"], input=script, capture_output=True, text=True, check=True)
    lines = run.stdout.splitlines()
    dimensions = [int(line[4:]) for line in lines if line.startswith("DIM ")]
    if "?" in run.stdout or "much slower" in run.stdout or len(dimensions) != 1 or dimensions[0] < 0:
        raise ArithmeticError(run.stdout[-2200:] + run.stderr[-1300:])
    context = dict(zip(map(str, variables), variables, strict=True))
    basis = [sp.sympify(line[2:].replace("^", "**"), locals=context) for line in lines if line.startswith("B ")]
    if not basis:
        raise ArithmeticError("Singular returned no finite basis.")
    return dimensions[0], basis


def output(row):
    return {"geometric_points": row["geometric_points"], "maps": sorted(row["maps"], key=encode)}


def worker(payload):
    if "CREST_SINGULAR" in os.environ:
        del os.environ["CREST_SINGULAR"]
    start = cpu()
    f = normalization(polynomial(payload["coefficients"]))[1]
    d, m, e = payload["stratum"]
    method = payload["method"]
    if method == "derivative":
        result = output(variable_maps(f, d, m, e))
    else:
        if method == "singular-formal-root":
            variables, equations, u, v = variable_incidence(f, d, m, e, derivative=False)
        else:
            variables, equations, u, v = original_incidence(f, d, m, e)
        dimension, basis = native_basis(payload["executable"], equations, variables)
        cert = rational_points(basis, variables)
        if cert["algebra_dimension"] != dimension:
            raise ArithmeticError("The native and local quotient dimensions disagree.")
        maps = []
        for point in cert["rational_points"]:
            assignment = dict(zip(variables, map(sp.Rational, point), strict=True))
            a = sp.Poly(u.as_expr().subs(assignment), X, domain=sp.QQ)
            b = sp.Poly(v.as_expr().subs(assignment), X, domain=sp.QQ)
            g = left_divide(f*b**2, a)
            if (g is None or g.degree() != m or g.gcd(g.diff()).degree() != 0
                    or g.compose(a) != f*b**2 or not a.diff().rem(b).is_zero):
                raise ArithmeticError("A native point failed the original map identity.")
            maps.append({"inner": coefficients(a), "outer": coefficients(g), "vertical": coefficients(b)})
        result = output({"geometric_points": cert["geometric_points"], "maps": maps})
    return {"solver_cpu": cpu() - start, "output": result}


def schedule(inputs, provenance):
    cases = [(name, provenance["cases"][name]["degree_stratum"]) for name in inputs
             if provenance["cases"][name]["study"] == "nonconstant"]
    rows = [(name, dimensions, method, repetition) for name, dimensions in cases
            for method in METHODS for repetition in range(1, REPEATS+1)]
    random.Random(SCHEDULE_SEED).shuffle(rows)
    return rows


def run_worker(payload):
    started = time.perf_counter()
    cpu_started = cpu()
    process = subprocess.Popen([sys.executable, "-B", str(__file__), "--worker"], start_new_session=True,
                               stdin=subprocess.PIPE, stdout=subprocess.PIPE, stderr=subprocess.PIPE, text=True)
    try:
        stdout, stderr = process.communicate(encode(payload).decode("utf-8"), timeout=LIMIT)
    except subprocess.TimeoutExpired:
        os.killpg(process.pid, signal.SIGKILL)
        process.communicate()
        return {"status": "timeout", "wall": time.perf_counter() - started, "process_cpu": cpu() - cpu_started}
    observations = {"wall": time.perf_counter()-started, "process_cpu": cpu()-cpu_started}
    if process.returncode:
        return {"status": "process_error", **observations, "error": stderr[-1800:]}
    return {"status": "completed", **observations, **json.loads(stdout)}


def references():
    records = read(HERE / "study" / "measurements.json")["rows"]
    answers = {}
    for row in records:
        if row["method"] == "derivative" and row["status"] == "completed":
            answer = output(row["output"])
            if answers.setdefault(row["case"], answer) != answer:
                raise ArithmeticError("Prior derivative repetitions disagree.")
    return answers


def summarize(rows):
    answers = references()
    result = {}
    for method in METHODS:
        candidates = [row for row in rows if row["method"] == method]
        statuses = {status: sum(row["status"] == status for row in candidates)
                    for status in {row["status"] for row in candidates}}
        successful = 0
        for row in candidates:
            if row["status"] == "completed":
                if row["output"] != answers[row["case"]]:
                    raise ArithmeticError(f"Fresh exact answer disagrees on {row['case']} ({method}).")
                successful += 1
        paired = {}
        for metric in ("solver_cpu", "process_cpu", "wall"):
            ratios = []
            for case in answers:
                both = [[r[metric] for r in rows if r["case"] == case and r["status"] == "completed"
                         and r["method"] == route] for route in ("derivative", method)]
                if all(len(a) == REPEATS for a in both):
                    ratios.append(statistics.median(both[1])/statistics.median(both[0]))
            paired[metric] = {"paired": len(ratios), "derivative_wins": sum(r > 1 for r in ratios),
                              "native_wins": sum(r < 1 for r in ratios),
                              "median_ratio": round(statistics.median(ratios), 6) if ratios else None}
        result[method] = {"statuses": statuses, "verified_rows": successful, "versus_derivative": paired}
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    group = parser.add_mutually_exclusive_group(required=True)
    group.add_argument("--worker", action="store_true")
    group.add_argument("--run", action="store_true")
    group.add_argument("--verify", action="store_true")
    parser.add_argument("--receipt", type=Path, default=PATH)
    args = parser.parse_args()
    if args.worker:
        print(json.dumps(worker(json.load(sys.stdin))))
        return
    source_files = [HERE / "study" / name for name in ("inputs.json", "provenance.json", "measurements.json", "freeze.json")]
    source_files.append(Path(__file__))
    source_hashes = {path.name: digest(path.read_bytes()) for path in source_files}
    inputs, provenance = (read(path) for path in source_files[:2])
    plan = schedule(inputs, provenance)
    if args.verify:
        prior = read(args.receipt)
        rows = prior["rows"]
        seen = [(r["case"], r["stratum"], r["method"], r["repetition"]) for r in rows]
        if (prior["source_hashes"] != source_hashes or seen != plan
                or prior["summary"] != summarize(rows)
                or prior["limit_seconds"] != LIMIT or prior["repetitions"] != REPEATS
                or prior["schedule_seed"] != SCHEDULE_SEED):
            raise ArithmeticError("The native control does not replay.")
        print(f"Verified {len(rows)} post-study native-control observations.")
        return
    if args.receipt.exists():
        raise FileExistsError(args.receipt)
    executable = os.environ.get("CREST_QUOTIENT_SINGULAR") or shutil.which("Singular")
    if not executable:
        raise RuntimeError("Install Singular or set CREST_QUOTIENT_SINGULAR.")
    available = sorted(os.sched_getaffinity(0))
    os.sched_setaffinity(0, {available[-1]})
    rows = []
    for i, (case, stratum, method, repetition) in enumerate(plan, 1):
        row = run_worker({**inputs[case], "stratum": stratum, "method": method, "executable": executable})
        rows.append({"case": case, "stratum": stratum, "method": method, "repetition": repetition, **row})
        if i % 20 == 0:
            print(f"Measured {i}/{len(plan)} native-control observations.", flush=True)
    if any(digest(path.read_bytes()) != source_hashes[path.name] for path in source_files):
        raise ArithmeticError("Source changed during native timing.")
    write(args.receipt, {"source_hashes": source_hashes,
                         "native_executable_sha256": digest(Path(executable).read_bytes()),
                         "selected_cpu": available[-1], "limit_seconds": LIMIT, "repetitions": REPEATS,
                         "schedule_seed": SCHEDULE_SEED, "rows": rows, "summary": summarize(rows)})


if __name__ == "__main__":
    main()
