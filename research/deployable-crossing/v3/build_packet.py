"""Build a deterministic offline packet with the unchanged v2 reference study."""
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
            contents["v3/"+relative.as_posix()]=path.read_bytes()
    v2=ROOT.parent/"v2"
    for path in sorted(v2.rglob("*")):
        if path.is_file() and path.suffix in SUFFIXES and "__pycache__" not in path.parts:
            contents["v2/"+path.relative_to(v2).as_posix()]=path.read_bytes()
    required={"v3/workbench.html","v3/integrated-report.json","v3/verification-release.json","v3/verification-motion.json",
      "v3/verification-assembly.json","v3/verification-site.json","v3/verification-viewer.json","v3/outcome-report.json"}
    if not required<=contents.keys(): raise ValueError("Incomplete packet: "+str(sorted(required-contents.keys())))
    release=json.loads(contents["v3/verification-release.json"])
    if release["checks_passed"] is not True: raise ValueError("Release verification is not successful")
    for name,digest in release["file_sha256"].items():
        if hashlib.sha256(contents["v3/"+name]).hexdigest()!=digest: raise ValueError("Changed file after release verification: "+name)
    contents["README.txt"]=("TILLAMOOK CROSSING V3\n\nExtract all files. On Windows double-click Open-Workbench.cmd.\n"
      "On any desktop open v3/workbench.html in a browser; no server/network is needed.\n\n"
      "The same reference assembly drives every study. Read v3/README.md and v3/integration-review.md.\n"
      "This is research software, not a bridge design approval, controller or public-use authorization.\n\n"
      "With Python 3.11+: python -B v3/verify_release.py\n"
      "Exact core checks use the standard library. Optional GIS replay needs GDAL/GEOS;\n"
      "JavaScript execution review needs Node and does not verify real browser/mobile rendering.\n\n"
      "The manifest checks file identity only. Keep v2 and v3 together; v3 imports pinned v2 dependencies.\n").encode()
    contents["Open-Workbench.cmd"]=b'@echo off\r\nstart "" "%~dp0v3\\workbench.html"\r\n'
    ledger={name:{"bytes":len(data),"sha256":hashlib.sha256(data).hexdigest()} for name,data in contents.items()}
    contents["MANIFEST.sha256.json"]=(json.dumps({"format":"crest-root-tillamook-v3-packet/1","design_sha256":release["design_sha256"],
      "scope":"File identity only; no physical approval","files":ledger},indent=2,sort_keys=True)+"\n").encode()
    target=Path(target);target.parent.mkdir(parents=True,exist_ok=True)
    with zipfile.ZipFile(target,"w",zipfile.ZIP_DEFLATED,compresslevel=9) as archive:
        for name,data in sorted(contents.items()):
            info=zipfile.ZipInfo(name,(2000,1,1,0,0,0));info.external_attr=0o100644<<16
            archive.writestr(info,data,compress_type=zipfile.ZIP_DEFLATED,compresslevel=9)
    return {"file":target.name,"members":len(contents),"bytes":target.stat().st_size,"sha256":hashlib.sha256(target.read_bytes()).hexdigest()}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__);parser.add_argument("output",type=Path)
    print(json.dumps(build(parser.parse_args().output),indent=2))
