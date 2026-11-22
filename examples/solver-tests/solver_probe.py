"""Exercise a geometric degree-three quotient solver on the Fermat benchmark family."""
import argparse
import json
from pathlib import Path
import subprocess
import sys
import tempfile

from compressed_lift import lift_step

PRIME = 13
ROOT = Path(__file__).resolve().parent
BASE_TERMS = ((4, 0, 0, 1), (0, 4, 0, 1), (0, 0, 4, 1))
DIRECTION = ((2, 0, 2, -5), (0, 2, 2, -4))


def positive(value):
    number = int(value)
    if number < 1:
        raise argparse.ArgumentTypeError("choose a positive integer")
    return number


def make_cases(delays):
    cases = {}
    for delay in [0, *sorted(set(delays))]:
        name = f"delay-{delay}" if delay else "fermat"
        terms = [list(t) for t in BASE_TERMS]
        if delay:
            terms.extend([i, j, k, c * PRIME**delay] for i, j, k, c in DIRECTION)
        cases[name] = {"delay": delay, "terms": terms}
    return cases


def starting_classes():
    pairs = json.loads((ROOT / "seeds.json").read_text(encoding="utf-8"))
    marked_points = {(state[0], state[1]) for pair in pairs for state in pair}
    direct_points = {(x, y) for x in range(PRIME) for y in range(PRIME)
                     if (x**4 + y**4 + 1) % PRIME == 0}
    if (len(pairs) != 16 or any(len(pair) != 2 for pair in pairs)
            or len(marked_points) != 32 or marked_points != direct_points):
        raise ValueError("the 16 stored pairs must cover all 32 distinct affine Fermat marks")
    return pairs


def trace(terms, pairs, precision):
    # Each input quartic has Z=1 on this frame. Keep both marks of every quotient.
    affine = {(i, j): c for i, j, _, c in terms}
    classes = pairs
    counts = [len(classes)]
    for exponent in range(1, precision):
        advanced = []
        for pair in classes:
            children = []
            for state in pair:
                rank, lifts = lift_step(affine, "I", state, PRIME**exponent, PRIME)
                if rank != 2 or len(lifts) > 1:
                    raise ValueError("benchmark left the regular two-parameter chart")
                children.extend(lifts)
            if len(children) == 1:
                raise ValueError("the paired marks disagree on survival")
            if children:
                advanced.append(children)
        classes = advanced
        counts.append(len(classes))
        if not classes:
            break
    return counts


def expected_counts(delay, precision):
    if delay and delay < precision:
        return [16] * delay + [0]
    return [16] * precision


def exact_count(case):
    # Oracle supplied by the Fermat classification and delayed-obstruction theorem.
    return 0 if case["delay"] else 16


def run_benchmark(cases):
    pairs = starting_classes()
    precision = 1 + max(case["delay"] for case in cases.values())
    moduli = [PRIME**n for n in range(1, precision + 1)]
    print("Candidate quotient classes after each local test")
    widths = [max(6, len(str(m))) for m in moduli]
    headings = " ".join(str(m).rjust(w) for m, w in zip(moduli, widths))
    print(f"{'curve':<12} {headings}   exact oracle")
    passed = True
    for name, case in cases.items():
        counts = trace(case["terms"], pairs, precision)
        correct = counts == expected_counts(case["delay"], precision)
        passed &= correct
        cells = [str(n) for n in counts] + ["-"] * (precision - len(counts))
        line = " ".join(x.rjust(w) for x, w in zip(cells, widths))
        print(f"{name:<12} {line}   {exact_count(case):>5}   {'PASS' if correct else 'FAIL'}")
    print("A dash follows an empty candidate set. The final column comes from the exact results.")
    return 0 if passed else 1


def export_cases(cases):
    return {"scope": "geometric", "map_degree": 3, "prime": PRIME,
            "cases": [{"id": name, "terms": case["terms"]} for name, case in cases.items()]}


def premature(cases, precision):
    # An intentionally wrong stopping rule, for exercising a new test installation.
    pairs = starting_classes()
    answers = {}
    for name, case in cases.items():
        last = trace(case["terms"], pairs, precision)[-1]
        answers[name] = {"scope": "geometric", "complete": True, "lower": last, "upper": last}
    return answers


def issue_for(answer, expected):
    if not isinstance(answer, dict):
        return "expected an object with scope, complete, lower and upper"
    if answer.get("scope") != "geometric":
        return "expected geometric scope"
    lo, hi = answer.get("lower"), answer.get("upper")
    if type(lo) is not int or type(hi) is not int or not 0 <= lo <= hi:
        return "bounds must be ordered nonnegative integers"
    if answer.get("complete") is not True:
        return f"incomplete (bounds {lo}..{hi})"
    if lo != hi:
        return f"claims completion with unequal bounds {lo}..{hi}"
    if lo != expected:
        return f"reported {lo}; expected {expected}"
    return None


def unique_keys(pairs):
    data = {}
    for key, value in pairs:
        if key in data:
            raise ValueError(f"duplicate JSON key {key!r}")
        data[key] = value
    return data


def check_answers(cases, submitted):
    if not isinstance(submitted, dict):
        raise ValueError("the outer JSON object must be keyed by case ID")
    failed = False
    for name, case in cases.items():
        issue = (issue_for(submitted[name], exact_count(case))
                 if name in submitted else "missing result")
        if issue:
            failed = True
            print(f"FAIL  {name}: {issue}")
        else:
            print(f"PASS  {name}: complete geometric count {exact_count(case)}")
    for name in sorted(submitted.keys() - cases.keys()):
        failed = True
        print(f"FAIL  {name}: unexpected case ID")
    return int(failed)


def invoke_solver(cases, command, timeout):
    # Each case gets its own equation. The solver cannot rely on an input case label.
    answers = {}
    with tempfile.TemporaryDirectory(prefix="quartic-benchmark-") as directory:
        equation = Path(directory) / "quartic.json"
        for name, case in cases.items():
            equation.write_text(json.dumps({"terms": case["terms"]}), encoding="utf-8")
            invocation = [str(equation) if arg == "{input}" else arg for arg in command]
            try:
                job = subprocess.run(invocation, capture_output=True, timeout=timeout)
                if job.returncode:
                    detail = job.stderr.decode("utf-8", errors="replace").strip()
                    raise ValueError(f"exit {job.returncode}: {detail[-300:]}")
                answers[name] = json.loads(job.stdout.decode("utf-8-sig"),
                                           object_pairs_hook=unique_keys)
            except (OSError, ValueError, subprocess.TimeoutExpired) as exc:
                print(f"{name}: solver error: {exc}", file=sys.stderr)
    return answers


def output_json(data, path):
    payload = json.dumps(data, indent=2, allow_nan=False) + "\n"
    if path is None:
        sys.stdout.write(payload)
    else:
        path.write_text(payload, encoding="utf-8")


def main(argv=None):
    parser = argparse.ArgumentParser(description=__doc__)
    subs = parser.add_subparsers(dest="command", required=True)
    for command, help_text in (
        ("run", "recompute the actual marked local lifts"),
        ("export", "emit the quartics to pass to your solver"),
        ("early", "simulate a solver that guesses from its last local count"),
        ("check", "check a saved set of solver reports"),
        ("solve", "run your solver on each quartic and check its reports"),
    ):
        p = subs.add_parser(command, help=help_text)
        p.add_argument("--delays", type=positive, nargs="+", default=[1, 2, 3],
                       metavar="K", help="positive perturbation exponents (default: 1 2 3)")
        if command in {"early", "export"}:
            p.add_argument("--output", type=Path, help="write JSON here instead of stdout")
        if command == "early":
            p.add_argument("--precision", type=positive, default=3, metavar="N",
                           help="stop at modulus 13**N (default: 3)")
        if command == "solve":
            p.add_argument("--timeout", type=positive, metavar="SECONDS",
                           help="optional time limit for each solver call")
            p.add_argument("solver", nargs=argparse.REMAINDER,
                           help="-- program and arguments; {input} receives the quartic JSON path")
        if command == "check":
            p.add_argument("answers", type=Path)
    argv = sys.argv[1:] if argv is None else argv
    args = parser.parse_args(argv or ["run"])
    cases = make_cases(args.delays)
    try:
        if args.command == "run":
            return run_benchmark(cases)
        if args.command == "check":
            with args.answers.open(encoding="utf-8-sig") as f:
                return check_answers(cases, json.load(f, object_pairs_hook=unique_keys))
        if args.command == "solve":
            command = args.solver[1:] if args.solver[:1] == ["--"] else args.solver
            if not command or "{input}" not in command:
                parser.error("solve needs a command and a separate {input} argument after --")
            return check_answers(cases, invoke_solver(cases, command, args.timeout))
        data = export_cases(cases) if args.command == "export" else premature(cases, args.precision)
        output_json(data, args.output)
        return 0
    except (ValueError, OSError, ArithmeticError) as exc:
        parser.exit(2, f"error: {exc}\n")


if __name__ == "__main__":
    sys.exit(main())
