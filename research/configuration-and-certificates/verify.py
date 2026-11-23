"""Recompute all companion research receipts, without changing them."""

from decimal import Decimal
import json
from pathlib import Path
import subprocess
import sys

HERE = Path(__file__).resolve().parent


def comparison_value(text):
    # Exact quantities are integers or rational strings. Floats are display approximations.
    return json.loads(
        text, parse_float=lambda value: Decimal(format(float(value), ".12g"))
    )


def main():
    for case in ("mechanisms", "sixbar", "spectra", "power"):
        directory = HERE / case
        command = [sys.executable, str(directory / "check.py")]
        data = HERE.parent.parent / "assay_selector" / "data"
        if case == "spectra" and data.is_dir():
            command.extend(["--source", str(data)])
        if not __debug__:
            raise SystemExit("Run Python without -O so the exact checks execute.")
        run = subprocess.run(
            command, cwd=directory, check=True, stdout=subprocess.PIPE, text=True
        )
        expected = comparison_value(
            (directory / "receipt.json").read_text(encoding="utf-8")
        )
        if comparison_value(run.stdout) != expected:
            raise SystemExit(f"{case}: computed result differs from receipt.json")
        print(f"PASS  {case}", flush=True)


if __name__ == "__main__":
    main()
