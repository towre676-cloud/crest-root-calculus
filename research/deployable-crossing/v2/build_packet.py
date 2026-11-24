#!/usr/bin/env python3
"""Build a deterministic offline research packet, excluding caches/scratch."""
import argparse
import hashlib
import json
from pathlib import Path
import zipfile

ROOT=Path(__file__).resolve().parent
SUFFIXES={".py",".json",".md",".html",".png",".geojson",".cmd",".svg"}


def build(target):
    contents={}
    for path in sorted(ROOT.rglob("*")):
        relative=path.relative_to(ROOT)
        if path.is_file() and path.suffix in SUFFIXES and not any(p.startswith(".") or p=="__pycache__" for p in relative.parts):
            contents[relative.as_posix()]=path.read_bytes()
    required={"README.md","design-review.md","workbench.html","site-comparison.png",
              "verify_release.py","data/atlas-subset.json","site-report.json","hazard-report.json",
              "engineering-report.json","verification-release.json"}
    if not required<=contents.keys():
        raise ValueError("Incomplete release: "+str(sorted(required-contents.keys())))
    manifest={name:{"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()} for name,data in contents.items()}
    contents["MANIFEST.sha256.json"]=(json.dumps({"format":"crest-root-tillamook-v2-packet/1",
         "scope":"File identity only; not engineering approval","files":manifest},indent=2,sort_keys=True)+"\n").encode()
    target=Path(target)
    target.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(target,"w",zipfile.ZIP_DEFLATED,compresslevel=9) as z:
        for name,data in sorted(contents.items()):
            info=zipfile.ZipInfo(name,(2000,1,1,0,0,0))
            info.external_attr=0o100644<<16
            z.writestr(info,data,compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
    return {"file":target.name,"members":len(contents),"bytes":target.stat().st_size,
            "sha256":hashlib.sha256(target.read_bytes()).hexdigest()}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("output",type=Path)
    print(json.dumps(build(parser.parse_args().output),indent=2))
