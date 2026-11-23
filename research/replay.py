"""Replay the 23 November research release from either a checkout or the ZIP."""
import json
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ELLIPTIC = HERE / "elliptic-quotients"


def run(script, receipt=None):
    result = subprocess.run([sys.executable, "-B", str(script)], cwd=script.parent,
                            capture_output=True, text=True, encoding="utf-8",
                            env={**os.environ, "PYTHONIOENCODING": "utf-8"})
    if result.returncode:
        raise RuntimeError(f"{script.parent.name}/{script.name} failed:\n{result.stderr}\n{result.stdout}")
    if receipt is not None:
        expected = json.loads((script.parent/receipt).read_text(encoding="utf-8"))
        if json.loads(result.stdout) != expected:
            raise ArithmeticError(f"The receipt for {script.parent.name} differs")
    print(f"Verified {script.parent.name}/{script.name}", flush=True)


def main():
    run(HERE.parent / "applications/atlas/original/replay.py")
    # The blind spectrum replay also recomputes both individual bielliptic proofs.
    run(ELLIPTIC / "blind-degree-spectra/verify.py", "receipt.json")
    run(ELLIPTIC / "involution-atlas/verify.py")
    run(ELLIPTIC / "primitive-degree-four/verify.py", "receipt.json")
    run(ELLIPTIC / "a3-spectrum/verify.py", "receipt.json")
    run(ELLIPTIC / "a3-spectrum/verify_cubic_lifting.py")
    run(ELLIPTIC / "primitive-degree-spectrum/verify.py")
    run(HERE / "projective-graphs/verify.py", "receipt.json")
    run(HERE / "cubic-tower-factors/verify.py")
    run(ELLIPTIC / "degree-four-incidence/verify_projective.py", "receipt.json")
    run(ELLIPTIC / "degree-four-jets/verify.py", "receipt.json")
    run(ELLIPTIC / "polarized-reduction/verify.py")
    run(HERE / "cubic-tower-higher-supports/prime_rank.py")
    run(HERE / "cubic-tower-higher-supports/mod_three.py")
    print("All fourteen research verification groups passed.")


if __name__ == "__main__":
    main()
