"""Deterministically package the complete research release and periodic study.

Run from any working directory. Only the two named outputs are built;
historical exports and the introductory PDF are read or left untouched.
No mathematics is rerun here: the archive checks verify bytes and manifests.
"""
import argparse
import ast
from collections import Counter
import hashlib
import json
from pathlib import Path, PurePosixPath
import re
import tempfile
from zipfile import ZIP_DEFLATED, ZipFile, ZipInfo


ROOT = Path(__file__).resolve().parents[1]
STAMP = (2026, 11, 23, 0, 0, 0)
DATE = "2026-11-23"
CURRENT_DATE = "2026-11-24"
EXPECTED_GROUPS = 36
ARCHIVED_FULL = ROOT / "downloads" / f"crest-root-open-research-{DATE}.zip"
PACKETS = {
    "full": ("crest-root-open-research-current.zip", "crest-root-research/"),
    "periodic": (f"crest-root-periodic-spectra-{DATE}.zip", "crest-root-periodic-spectra/"),
}
SKIP_DIRS = {"__pycache__", "node_modules", "venv", "build", "dist"}
SKIP_SUFFIXES = {
    ".pyc", ".pyo", ".o", ".obj", ".a", ".so", ".dll", ".dylib",
    ".exe", ".class", ".wasm", ".bin", ".out", ".log", ".tmp", ".bak",
    ".zip", ".tar", ".gz", ".tgz", ".bz2", ".xz", ".7z",
}
INTRODUCTION = "downloads/crest-root-introduction.pdf"


def digest(data):
    return hashlib.sha256(data).hexdigest()


def source_paths(kind):
    directories = ("research", "applications/atlas") if kind == "full" else ("research/periodic-spectra",)
    result = []
    for directory in directories:
        base = ROOT / directory
        if not base.is_dir():
            raise FileNotFoundError(base)
        for path in base.rglob("*"):
            relative = path.relative_to(ROOT)
            if any(part in SKIP_DIRS or part.startswith(".") for part in relative.parts[:-1]):
                continue
            if path.name.startswith(".") and path.name != ".nojekyll":
                continue
            if path.is_symlink():
                raise ValueError(f"Archive source must not be a symlink: {relative}")
            if not path.is_file() or path.suffix.lower() in SKIP_SUFFIXES:
                continue
            with path.open("rb") as handle:
                magic = handle.read(4)
            if magic.startswith((b"\x7fELF", b"MZ")) or magic in {
                b"\xfe\xed\xfa\xce", b"\xce\xfa\xed\xfe",
                b"\xfe\xed\xfa\xcf", b"\xcf\xfa\xed\xfe",
            }:
                continue
            result.append(path)
    if kind == "full":
        result.append(ROOT / INTRODUCTION)
    return sorted(result, key=lambda path: path.relative_to(ROOT).as_posix())


def replay_group_count():
    tree = ast.parse((ROOT / "research/replay.py").read_text(encoding="utf-8"))
    main = next(node for node in tree.body if isinstance(node, ast.FunctionDef) and node.name == "main")
    return sum(isinstance(node, ast.Call) and isinstance(node.func, ast.Name)
               and node.func.id == "run" for node in ast.walk(main))


def readiness():
    problems = []
    count = replay_group_count()
    if count != EXPECTED_GROUPS:
        problems.append(f"Expected {EXPECTED_GROUPS} replay groups; found {count}.")
    study = ROOT / "research/periodic-spectra"
    for name in ("README.md", "results.md", "native-results.md", "proof.md", "spectral-curve.md"):
        if not (study / name).is_file():
            problems.append(f"Missing periodic study document: {name}.")
    for name in ("results.json", "native-results.json"):
        path = study / name
        if not path.is_file():
            problems.append(f"Missing completed study: {name}.")
            continue
        data = json.loads(path.read_text(encoding="utf-8"))
        if data.get("completed") is not True or data.get("dependencies_unchanged") is not True:
            problems.append(f"Study is incomplete or lacks its final dependency check: {name}.")
    return count, problems


def packet_readme(kind, groups):
    if kind == "periodic":
        return """# Crest–Root periodic spectra, 23 November 2026

This standalone packet contains the exact periodic Jacobi-ring challenge,
proofs, source code, complete measurements, frozen protocols and references,
the compiled-PARI follow-up, and independent verification.

Read [the study guide](research/periodic-spectra/README.md),
[the original results](research/periodic-spectra/results.md), and
[the compiled-PARI comparison](research/periodic-spectra/native-results.md).
The [spectral-curve companion](research/periodic-spectra/spectral-curve.md)
proves and checks the classical normalization identity for repeated cells.

After extracting, open a terminal in this folder. Windows commands:

```sh
py -m pip install "sympy>=1.12,<2"
py research/periodic-spectra/audit_native.py
py research/periodic-spectra/spectral_curve.py
```

On macOS or Linux, substitute `python3` for `py`. The stored study audit runs
offline and independently recomputes every output coefficient. No compiled
solver binary is needed to check the stored mathematics and measurements.
New timing runs require the documented C/C++ compilers, GMP, PARI/GP and a
POSIX environment; use the study guide for exact build and run commands.

The archive's MANIFEST.json records the size and SHA-256 of every other file.
The complete research packet additionally contains all 32 research
verification groups and the introductory PDF.
"""
    return f"""# Crest–Root open research, {CURRENT_DATE}

Read [the research guide](research/README.md), the
[57-item register](research/backlog-2026-11-23.md), or the
[nine-page introductory essay](downloads/crest-root-introduction.pdf).
This packet preserves the earlier research content and adds the
[periodic spectral challenge](research/periodic-spectra/README.md), its exact
proofs, original and compiled-PARI comparisons, and normalization examples.
The [branch reconciliation](research/branch-reconciliation-2026-11-23.md)
records the integration of all research branches and publication source.

It includes the canonical 148-section monograph, finite hyperelliptic map
reconstruction, local-to-global certificates, the equal-count obstruction,
earlier solver comparisons, the marked chart-I solver, point-free
Riemann–Roch, genus-two correspondences, genus-37 rank certificates, quartic
degree spectra, the historical atlas and higher-support tower certificates.
The [coefficient-only quotient guide](research/coefficient-quotients/README.md)
includes the full fixed-degree search, degree-two covers with arbitrary
infinity fibres, and the archived comparisons with Singular and PARI/GP.
The [cyclic-cover study](research/cyclic-cover-reconstruction/README.md)
extends source-only reconstruction to mixed Kummer branches and records a
new native Singular comparison.

After extracting, open a terminal in this folder. Windows commands:

```sh
py -m pip install numpy "sympy>=1.12,<2"
py research/replay.py --workers 4
```

On macOS or Linux, substitute `python3` for `py`. The replay launches the
**{groups} verification groups** listed in research/replay.py. Default checks
run offline. With PARI/GP installed, the rank replay also regenerates the
number-field and local arithmetic behind the upper bounds. Without GP, it
checks the recorded descent matrices and lower-bound divisor certificates.

The periodic study can be audited separately:

```sh
py research/periodic-spectra/audit_native.py
py research/periodic-spectra/spectral_curve.py
```

The finite-map certificates require Python and SymPy; Singular accelerates
new general eliminations. Repeating measured studies requires their stated
tools and a POSIX environment. All recorded outcomes, exact output checks,
resource limits, dependency hashes and timing provenance remain included.
See the [earlier solver instructions](research/solver-comparison/README.md)
and [periodic study guide](research/periodic-spectra/README.md).

The archive's MANIFEST.json records the size and SHA-256 of every other file.
The repository-native builder is included. To rebuild both packets from
this extracted source or a repository checkout:

```sh
py research/build_packet.py --list
py research/build_packet.py
py research/build_packet.py --check
```

Fixed entry ordering, timestamps, permissions and compression settings make
the output deterministic for unchanged inputs and the same Python/zlib
toolchain. The builder verifies every member against its source bytes and
an extracted copy. It excludes caches, executables, logs and archive
recursion. Existing historical download exports and the source PDF are
left untouched.
"""


def entries(kind, groups):
    content = {path.relative_to(ROOT).as_posix(): path.read_bytes() for path in source_paths(kind)}
    content["README.md"] = packet_readme(kind, groups).encode("utf-8")
    manifest = {
        "schema": 1, "packet": kind, "release_date": CURRENT_DATE if kind == "full" else DATE,
        "manifest_scope": "Every archive file except MANIFEST.json itself; paths are relative to the packet root.",
        "files": {name: {"bytes": len(data), "sha256": digest(data)}
                  for name, data in sorted(content.items())},
    }
    content["MANIFEST.json"] = (json.dumps(manifest, indent=2, ensure_ascii=False) + "\n").encode("utf-8")
    return content


def check_previous_full_archive(target, content):
    """Reject dropped historical members or a changed introductory PDF."""
    if not target.exists():
        return {"previous_archive": "absent"}
    prefix = PACKETS["full"][1]
    with ZipFile(target) as archive:
        old = {name[len(prefix):]: archive.read(name) for name in archive.namelist()
               if not name.endswith("/") and name.startswith(prefix)}
        if len(old) != sum(not info.is_dir() for info in archive.infolist()):
            raise ArithmeticError("Previous packet has an unexpected root or duplicate member.")
    missing = sorted(set(old) - set(content))
    if missing:
        raise ArithmeticError("Previous research files would be omitted: " + ", ".join(missing))
    if INTRODUCTION in old and old[INTRODUCTION] != content[INTRODUCTION]:
        raise ArithmeticError("The preserved introductory PDF changed.")
    return {"previous_files_retained": len(old),
            "previous_files_byte_identical": sum(content[name] == data for name, data in old.items()),
            "updated_existing_paths": sorted(name for name, data in old.items() if content[name] != data)}


def verify_archive(target, prefix, content):
    with ZipFile(target) as archive:
        expected_names = [prefix + name for name in sorted(content)]
        if archive.namelist() != expected_names or archive.testzip() is not None:
            raise ArithmeticError("Archive entries, ordering or CRC checks disagree.")
        for info in archive.infolist():
            name = info.filename
            if (PurePosixPath(name).is_absolute() or ".." in PurePosixPath(name).parts
                    or info.date_time != STAMP or info.external_attr != 0o100644 << 16):
                raise ArithmeticError("Unsafe archive member or nondeterministic metadata.")
            relative = name[len(prefix):]
            if archive.read(info) != content[relative]:
                raise ArithmeticError("Archived bytes differ: " + relative)
        with tempfile.TemporaryDirectory(prefix="crest-packet-check-") as directory:
            extraction = Path(directory)
            archive.extractall(extraction)
            for name, data in content.items():
                if (extraction / prefix / name).read_bytes() != data:
                    raise ArithmeticError("Extracted bytes differ: " + name)


def write_archive(target, prefix, content):
    target.parent.mkdir(parents=True, exist_ok=True)
    with tempfile.TemporaryDirectory(prefix=".crest-package-", dir=target.parent) as directory:
        temporary = Path(directory) / target.name
        with ZipFile(temporary, "w", compression=ZIP_DEFLATED, compresslevel=9) as archive:
            for name, data in sorted(content.items()):
                info = ZipInfo(prefix + name, STAMP)
                info.create_system = 3
                info.compress_type = ZIP_DEFLATED
                info.external_attr = 0o100644 << 16
                archive.writestr(info, data, compress_type=ZIP_DEFLATED, compresslevel=9)
        verify_archive(temporary, prefix, content)
        temporary.replace(target)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--only", choices=("full", "periodic", "both"), default="both")
    parser.add_argument("--output-dir", type=Path, default=ROOT / "downloads")
    mode = parser.add_mutually_exclusive_group()
    mode.add_argument("--list", action="store_true", help="Inspect selected paths and readiness without writing archives")
    mode.add_argument("--check", action="store_true", help="Verify existing archives against current source without replacing them")
    options = parser.parse_args()
    kinds = list(PACKETS) if options.only == "both" else [options.only]
    groups, problems = readiness()
    if options.list:
        print(json.dumps({"verification_groups": groups, "ready": not problems,
                          "readiness_problems": problems,
                          "packets": {kind: {"source_files": len(source_paths(kind)),
                                             "paths": [p.relative_to(ROOT).as_posix() for p in source_paths(kind)]}
                                      for kind in kinds}}, indent=2))
        return
    if problems:
        raise SystemExit("Packet is not ready:\n" + "\n".join(problems))
    rows = re.findall(r"^\| ([A-Z]\d{2}) \| (P[0-3]) \|",
                      (ROOT / "research/backlog-2026-11-23.md").read_text(encoding="utf-8"), re.M)
    if (len(rows) != 57 or len({key for key, _ in rows}) != 57
            or Counter(priority for _, priority in rows) != {"P0": 8, "P1": 25, "P2": 22, "P3": 2}):
        raise ArithmeticError("The original 57-question register changed its identifiers or priority counts.")
    report = []
    for kind in kinds:
        filename, prefix = PACKETS[kind]
        target = options.output_dir / filename
        content = entries(kind, groups)
        preservation = {}
        if kind == "full":
            preservation = check_previous_full_archive(target, content)
            preservation["historical_baseline"] = check_previous_full_archive(ARCHIVED_FULL, content)
        if options.check:
            verify_archive(target, prefix, content)
        else:
            write_archive(target, prefix, content)
        report.append({"packet": kind, "path": str(target), "entries": len(content),
                       "bytes": target.stat().st_size, "sha256": digest(target.read_bytes()),
                       "verified_source_and_extracted_bytes": True, **preservation})
    print(json.dumps({"verification_groups": groups, "archives": report}, indent=2))


if __name__ == "__main__":
    main()
