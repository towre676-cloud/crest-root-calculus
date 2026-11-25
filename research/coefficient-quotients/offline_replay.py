"""Replay archived quotient studies without the POSIX timing dependency."""

import argparse
from pathlib import Path
import runpy
import sys
from types import ModuleType


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("study", choices=("frozen", "native"))
    study = parser.parse_args().study

    # The frozen timing sources import resource, which is absent on Windows.
    # The read-only modes never call the timing API.
    try:
        import resource  # noqa: F401
    except ModuleNotFoundError:
        sys.modules["resource"] = ModuleType("resource")

    name, args = ("verify.py", []) if study == "frozen" else ("native_control.py", ["--verify"])
    script = Path(__file__).with_name(name)
    sys.argv = [str(script), *args]
    runpy.run_path(str(script), run_name="__main__")


if __name__ == "__main__":
    main()
