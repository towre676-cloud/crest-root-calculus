"""Verify the complete stored research packet and its independent certificates.

Python 3.11+ standard library suffices for the exact core. Optional GIS and
JavaScript replays require GDAL/GEOS and Node respectively. No physical
capacity, road opening or real-world benefit follows from this verifier.
"""
import argparse
from concurrent.futures import ThreadPoolExecutor
from copy import deepcopy
import hashlib
import json
from pathlib import Path
import subprocess
import sys

import integrate

ROOT=Path(__file__).resolve().parent


def execute(script,args=()):
    result=subprocess.run([sys.executable,"-B",str(ROOT/script),*args],cwd=ROOT,capture_output=True,text=True,check=True)
    receipt=json.loads(result.stdout)
    if not (receipt.get("checks_passed") is True or receipt.get("status")=="passed"):
        raise AssertionError("Verifier did not report success: "+script)
    return receipt


def validate_saved():
    saved=integrate.read(ROOT/"integrated-report.json")
    if saved!=integrate.build(): raise AssertionError("Integrated report or its source/report bindings are stale")
    design=integrate.load_design(); reports=saved["reports"]
    rejected=[]
    edits={
      "different_design":lambda r:r["motion"].update(design_sha256="0"*64),
      "different_geometry":lambda r:r["site"].update(geometry_and_mass_source_sha256="0"*64),
      "different_structural_capture":lambda r:r["motion"].update(loaded_capture_source_sha256="0"*64),
      "different_sweep":lambda r:r["site"].update(motion_report_sha256="0"*64),
      "different_outcome_site_source":lambda r:r["outcome"].update(site_report_sha256="0"*64),
      "promoted_public_opening":lambda r:r["operation"]["public_opening"].update(allowed=True),
    }
    for name,edit in edits.items():
        bad=deepcopy(reports);edit(bad)
        try: integrate.validate(design,bad)
        except ValueError: rejected.append(name)
        else: raise AssertionError("Bad integration was accepted: "+name)
    return saved,rejected


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--replay",action="store_true",help="Also recompute the full exact structural and motion reports")
    parser.add_argument("--with-gis",action="store_true",help="Run independent GDAL/GEOS geometry checks with this Python")
    parser.add_argument("--with-viewer",action="store_true",help="Run static and Node DOM/Canvas checks; not a browser-rendering test")
    args=parser.parse_args(); saved,rejected=validate_saved()
    jobs=[("assembly","verify_assembly.py",("--replay",) if args.replay else ()),
          ("motion","verify_motion.py",("--replay",) if args.replay else ()),
          ("evidence","verify_evidence.py",()),("outcomes","verify_outcomes.py",())]
    with ThreadPoolExecutor(max_workers=4) as pool:
        futures={name:pool.submit(execute,script,options) for name,script,options in jobs}
        receipts={name:future.result() for name,future in futures.items()}
    for name,receipt in receipts.items():
        (ROOT/("verification-"+name+".json")).write_text(json.dumps(receipt,indent=2)+"\n",encoding="utf-8")
    optional={"gis":"not_rerun_by_this_invocation","viewer":"not_rerun_by_this_invocation"}
    for enabled,name,script in ((args.with_gis,"gis","verify_site.py"),(args.with_viewer,"viewer","independent_review.py")):
        if enabled:
            subprocess.run([sys.executable,"-B",str(ROOT/script)],cwd=ROOT,check=True)
            optional[name]="passed"
    files={}
    for path in sorted(ROOT.rglob("*")):
        if path.is_file() and path.suffix in {".py",".json",".md",".html",".cmd"} and path.name!="verification-release.json" and "__pycache__" not in path.parts:
            files[path.relative_to(ROOT).as_posix()]=hashlib.sha256(path.read_bytes()).hexdigest()
    result={"schema_version":1,"format":"crest-root-integrated-release-verification/3","checks_passed":True,
      "design_sha256":saved["design_sha256"],"integration_checks":saved["integration_checks"],
      "integration_negative_controls":rejected,"core_verifiers":{name:{"passed":True,"receipt":"verification-"+name+".json"} for name in receipts},
      "full_solver_replay":args.replay,"optional_replays":optional,"python":sys.version.split()[0],"file_sha256":files,
      "physical_approval":False,"real_browser_rendering_verified":False,"mobile_touch_verified":False,
      "scope":"Reproducible software-model identities and independent checks; no physical bridge, site, operation or casualty benefit is approved"}
    (ROOT/"verification-release.json").write_text(json.dumps(result,indent=2)+"\n",encoding="utf-8")
    print(json.dumps({key:result[key] for key in ("checks_passed","integration_checks","core_verifiers","full_solver_replay","optional_replays","physical_approval")},indent=2))


if __name__=="__main__": main()
