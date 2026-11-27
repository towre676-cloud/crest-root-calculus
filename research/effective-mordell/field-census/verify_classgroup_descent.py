"""Replay 539 dedicated field certificates with cypari2 and the included arithmetic module."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from classgroup_descent import classify, pari


def replay(row):
    result = classify(row["polynomial"])
    expected = row["result"]
    if json.loads(json.dumps(result)) != expected:
        raise ArithmeticError(f"Certificate disagreement for {row['polynomial']}")
    return result["status"]


def quadratic_check():
    verified = 0
    for d in range(-200, 201):
        if d in (0, 1) or not pari.issquarefree(abs(d)):
            continue
        result = classify([-d, 0, 1])
        for sign in (-1, 1):
            # The untwisted rational quotient has rank zero.
            twist = pari.ellinit([0, sign * d**3])
            lower = int(twist.ellrank()[0])
            if lower > result["rank_upper_bounds"][str(sign)]:
                raise ArithmeticError(f"Known rank exceeds the bound at {d}, {sign}")
        verified += 1
    if verified != 243:
        raise ArithmeticError("Wrong quadratic-field panel cardinality")
    print("Independent rational-twist rank checks passed on", verified, "quadratic fields.")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--workers", type=int, default=1)
    cli.add_argument("--quadratic-check", action="store_true")
    args = cli.parse_args()
    if args.workers < 1:
        cli.error("The worker count must be positive.")
    rows = json.loads(Path(__file__).with_name("classgroup-field-receipts.json").read_text())
    if args.workers == 1:
        totals = Counter(map(replay, rows))
    else:
        with ProcessPoolExecutor(args.workers) as executor:
            totals = Counter(executor.map(replay, rows))
    if totals != {"complete": 123, "partial": 416}:
        raise ArithmeticError(f"Unexpected classification totals {totals}")
    print("539 field certificates verified.", dict(totals))
    if args.quadratic_check:
        quadratic_check()
