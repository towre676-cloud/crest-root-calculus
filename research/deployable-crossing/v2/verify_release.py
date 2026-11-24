#!/usr/bin/env python3
"""Offline focused v2 verification. Standard-library default; no network."""
import argparse
import hashlib
import json
from pathlib import Path
import subprocess
import sys

from atlas_adapter import load_subset
from bridge_examples import build_examples
from engineering_study import run as engineering_report
from verify_engineering import verify as engineering_checks
from verify_truss import verify_all

ROOT=Path(__file__).resolve().parent


def read(name): return json.loads((ROOT/name).read_text(encoding="utf-8"))
def sha(path): return hashlib.sha256(path.read_bytes()).hexdigest()


def verify(recompute_sites=False):
    groups=[]
    atlas=load_subset()
    groups.append({"group":"Atlas subset content/terrain integrity","status":"passed","counts":atlas["counts"]})
    truss=verify_all()
    assert truss==read("verification-truss.json"),"Stored truss receipt differs"
    groups.append({"group":"Exact truss independent controls","status":"passed"})
    bridges=build_examples()
    assert bridges==read("examples/bridge-example-report.json"),"Stored bridge examples differ"
    groups.append({"group":"Bridge sections and virtual-work controls","status":"passed","checks":bridges["independent_checks"]})
    engineering=engineering_checks()
    assert engineering==read("verification-engineering.json"),"Stored engineering checks differ"
    assert engineering_report()==read("engineering-report.json"),"Stored support/join/access report differs"
    groups.append({"group":"Exact support paths, capture transforms, access and geometry","status":"passed","checks":engineering["check_count"]})
    independent=read("verification-engineering-independent.json")
    for name,digest in independent["reviewed_sha256"].items():
        assert sha(ROOT/name)==digest,"Independent engineering review is stale: "+name
    site_review=read("verification-site-independent.json")
    for name,digest in site_review["associated_source_code_sha256"].items():
        assert sha(ROOT/name)==digest,"Independent site review is stale: "+name
    assert sha(ROOT/"site-report.json")==site_review["input_hashes"]["site_report_sha256"]
    assert sha(ROOT/"data/atlas-subset.json")==site_review["input_hashes"]["atlas_subset_sha256"]
    groups.append({"group":"Independent engineering and GIS review byte bindings","status":"passed"})
    site=read("site-report.json")
    assert site["atlas_content_sha256"]==atlas["content_sha256"],"Site report uses a different Atlas"
    if recompute_sites:
        from site_study import run_study
        regenerated=json.loads(json.dumps(run_study(atlas),ensure_ascii=False,allow_nan=False))
        assert regenerated==site,"Placement scenario replay differs"
    groups.append({"group":"Placement scenarios","status":"recomputed" if recompute_sites else "bound saved report",
                   "scenario_count":len(site["scenarios"]),"summary":site["summary"]})
    checked=subprocess.run([sys.executable,"-S",str(ROOT/"hazard_screen.py"),"--check"],
                           check=True,capture_output=True,text=True)
    groups.append({"group":"Hazard source/report integrity and query coverage","status":"passed",
                   "receipt":json.loads(checked.stdout)})
    # Every shipped example output must still solve its own input exactly.
    from truss import solve_truss
    for name in ("service-24m-plane","prejoin-12m-plane"):
        assert solve_truss(read("examples/"+name+".json"))==read("examples/"+name+"-result.json"),name
    groups.append({"group":"Stored complete truss outputs","status":"passed","case_count":2})
    if (ROOT/"MANIFEST.sha256.json").exists():
        manifest=read("MANIFEST.sha256.json")["files"]
        for name,entry in manifest.items():
            path=ROOT/name
            assert path.resolve().is_relative_to(ROOT),name
            assert path.stat().st_size==entry["bytes"] and sha(path)==entry["sha256"],name
        groups.append({"group":"Extracted archive file identity","status":"passed","file_count":len(manifest)})
    return {"format":"crest-root-tillamook-release-verification/1","status":"passed",
            "python":sys.version.split()[0],"stdlib_default":True,"network_required":False,
            "recomputed_sites":recompute_sites,"groups":groups,
            "scope":"Reproducibility and stated model correctness only; no physical tests, site approval or capacity rating"}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--recompute-sites",action="store_true")
    parser.add_argument("--output",type=Path,help="Optional new receipt path; default prints without modifying the packet")
    args=parser.parse_args()
    result=verify(args.recompute_sites)
    text=json.dumps(result,indent=2)+"\n"
    if args.output: args.output.write_text(text,encoding="utf-8")
    print(text,end="")
