"""Replay the 23 November research release from either a checkout or the ZIP."""
import argparse
import json
import os
from pathlib import Path
import subprocess
import sys


HERE = Path(__file__).resolve().parent
ELLIPTIC = HERE / "elliptic-quotients"


def run(script, receipt=None, args=()):
    result = subprocess.run([sys.executable, "-B", str(script), *args], cwd=script.parent,
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
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--workers", type=int, default=min(4,os.cpu_count() or 1),
                        help="Number of processes for independent local polynomial checks")
    parser.add_argument("--recount", action="store_true",
                        help="Recompute every stored higher-support matrix from the curve equations")
    options=parser.parse_args()
    if options.workers<1:
        parser.error("Choose at least one worker")
    local_args=["--workers",str(options.workers)]
    if options.recount:
        local_args.append("--recount")
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
    run(ELLIPTIC / "flex-fiber-benchmark/verify.py")
    run(ELLIPTIC / "cartier-degree-four/verify.py")
    run(HERE / "cubic-tower-cartier/verify.py")
    run(HERE / "cubic-tower-geometric-separation/verify.py",args=local_args)
    run(HERE / "cubic-tower-size-five/verify.py",args=local_args)
    run(HERE / "cubic-tower-size-five/six_index.py",args=local_args)
    run(HERE / "cubic-tower-stage-seven/verify.py",args=local_args)
    run(HERE / "cubic-tower-size-five/high_orbits.py",args=["--workers",str(options.workers)])
    run(HERE / "chart-one-cofactors/verify.py", args=["--workers",str(options.workers)])
    run(HERE / "point-free-riemann-roch/verify.py")
    run(HERE / "local-global-correspondences/verify.py")
    run(ELLIPTIC / "index-four-descent/verify.py", "receipt.json")
    run(HERE / "genus37-rank/verify.py", "receipt.json")
    run(HERE / "monograph-publication/verify.py", "receipt.json")
    run(HERE / "bounded-correspondences/verify.py")
    run(HERE / "local-pattern-obstruction/verify.py", "receipt.json")
    run(HERE / "solver-comparison/verify.py")
    run(HERE / "periodic-spectra/audit_native.py")
    print("All 32 research verification groups passed.")


if __name__ == "__main__":
    main()
