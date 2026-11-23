"""Build the deterministic, self-contained crossing feasibility packet."""

import argparse
import hashlib
import json
from pathlib import Path
import zipfile


HERE = Path(__file__).resolve().parent
FILES = (
    "README.md", "proof.md", "sources.md", "model.py", "verify.py", "study.py",
    "example.json", "report.json", "verification.json", "verification-symbolic.json",
    "build_packet.py",
)


def build(target):
    contents = {name: (HERE / name).read_bytes() for name in FILES}
    manifest = {name: {"bytes": len(data), "sha256": hashlib.sha256(data).hexdigest()}
                for name, data in sorted(contents.items())}
    contents["MANIFEST.sha256.json"] = (
        json.dumps({"scope": "Archive content identity, not engineering validation",
                    "files": manifest}, indent=2, sort_keys=True) + "\n").encode("utf-8")
    target.parent.mkdir(parents=True, exist_ok=True)
    with zipfile.ZipFile(target, "w", compression=zipfile.ZIP_DEFLATED, compresslevel=9) as archive:
        for name, data in sorted(contents.items()):
            entry = zipfile.ZipInfo(name, (2000, 1, 1, 0, 0, 0))
            entry.compress_type = zipfile.ZIP_DEFLATED
            entry.external_attr = 0o100644 << 16
            archive.writestr(entry, data, compress_type=zipfile.ZIP_DEFLATED, compresslevel=9)
    return {"file": target.name, "members": len(contents), "bytes": target.stat().st_size,
            "sha256": hashlib.sha256(target.read_bytes()).hexdigest()}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output", type=Path)
    args = parser.parse_args()
    print(json.dumps(build(args.output), indent=2, sort_keys=True))


if __name__ == "__main__":
    main()
