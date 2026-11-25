"""Freeze the discovery and generation code before drawing held-out seeds."""

import argparse
from datetime import datetime, timezone
from pathlib import Path
import secrets
import subprocess

from benchmark import HERE, digest, write


ROOT = HERE.parent.parent
SOURCE_FILES = (
    "discover.py", "baselines.py", "decompose.gp", "panel.py", "benchmark.py", "freeze.py",
    "../bounded-correspondences/finite.py", "../bounded-correspondences/composition.py",
    "../bounded-correspondences/maps.py", "../solver-comparison/map_worker.py",
    "../solver-comparison/worker.py", "../periodic-spectra/spectral_curve.py",
)


def source_paths():
    return [(HERE / name).resolve().relative_to(ROOT).as_posix() for name in SOURCE_FILES]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("directory", type=Path)
    args = parser.parse_args()
    args.directory.mkdir(parents=True, exist_ok=True)
    target = args.directory / "freeze.json"
    if target.exists():
        raise FileExistsError(target)
    head = subprocess.check_output(["git", "rev-parse", "HEAD"], cwd=ROOT, text=True).strip()
    frozen = {name: digest((ROOT / name).read_bytes()) for name in source_paths()}
    write(target, {"discovery_commit": head, "sha256": frozen,
                   "created_utc": datetime.now(timezone.utc).isoformat(),
                   "panel_seed": secrets.randbits(64), "schedule_seed": secrets.randbits(64)})


if __name__ == "__main__":
    main()
