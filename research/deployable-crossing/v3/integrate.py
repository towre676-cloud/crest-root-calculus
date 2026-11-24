"""Compose the shared design reports without promoting model results to approval.

Use --refresh to regenerate components in dependency order. The default binds
existing component reports and builds the review packet's integrated record.
"""
import argparse
from fractions import Fraction as F
import hashlib
import json
import os
from pathlib import Path
import subprocess
import sys

from assembly import design_digest, load_design
from outcomes import OUTCOMES, solve_outcome

ROOT=Path(__file__).resolve().parent
COMPONENTS=("assembly","motion","site","evidence","operation","outcome")


def read(path): return json.loads(Path(path).read_text(encoding="utf-8"))
def file_sha(path): return hashlib.sha256(Path(path).read_bytes()).hexdigest()
def require(condition,message):
    if not condition: raise ValueError(message)


def validate(design,reports):
    digest=design_digest(design)
    require(set(reports)==set(COMPONENTS),"All six reports are required")
    for key,report in reports.items():
        require(report.get("design_sha256")==digest,key+" has a different design identity")
    a,m,s,e,o,h=(reports[k] for k in COMPONENTS)
    geometry=design_digest(a["assembly"])
    require(a["assembly"]["design"]==design,"Assembly embeds different design values")
    for key in ("assembly","motion","site"):
        require(reports[key].get("geometry_and_mass_source_sha256")==geometry,key+" uses different geometry/mass")
    require(m.get("structural_report_sha256")==design_digest(a),"Motion uses another structural report")
    require(m.get("loaded_capture_source_sha256")==design_digest(a["capture_reference"]),"Capture uses another structural state")
    require(s.get("motion_report_sha256")==design_digest(m),"Site uses another continuous sweep")
    for key,source in (("assembly_receipt_content_sha256",a),("motion_receipt_content_sha256",m),("evidence_report_content_sha256",e)):
        require(o["receipt_bindings"].get(key)==design_digest(source),"Operation receipt binding changed: "+key)
    for key in ("assembly","motion","site","operation"):
        require(reports[key].get("physical_approval") is False,key+" cannot establish physical approval")
    require(o["public_opening"]["allowed"] is False,"Physical opening cannot be enabled in this research packet")
    require(e["physical_measurement_count"]==0 and e["source_kind"]=="synthetic_fixture","Evidence fixture provenance changed")
    nodes=set(a["assembly"]["physical_nodes"])
    mass=F(a["assembly"]["mass_ledger"]["bridge_mass_kg"])
    for frame in m["frames"]:
        require(set(frame["world_nodes_m"])==nodes,"A motion frame changes physical node identities")
        require(sum(F(half["mass_kg"]) for half in frame["halves"])==mass,"A motion frame changes moving mass")
    require(len(s["sites"])==6 and len(s["placements"])==18 and len(s["scenarios"])==216,"Reference study coverage changed")
    require(s["summary"]["approved_sites"]==0,"Map screens cannot approve a site")
    require(type(h.get("schema_version")) is int and h["schema_version"]==1 and h.get("physical_approval") is False,"Unsupported outcome receipt")
    require(h.get("site_report_sha256")==design_digest(s),"Outcome ledger uses another site report")
    ids=[row["site_id"] for row in h["site_results"]]
    require(len(ids)==len(set(ids)) and set(ids)=={row["id"] for row in s["sites"]},"Outcome sites differ from mapped references")
    for row in h["site_results"]:
        require(len(row["outcomes"])==len(OUTCOMES) and {v["outcome"] for v in row["outcomes"]}==set(OUTCOMES),"Outcomes must remain separate at every site")
        for entry in row["outcomes"]:
            case,result=entry["input"],entry["result"]
            require(case["site_id"]==row["site_id"] and case["outcome"]==entry["outcome"] and case["design_sha256"]==digest,"Outcome input has another design/site/outcome")
            require(solve_outcome(case)==result,"Outcome result differs from its evidence")
            require(result["status"]=="unknown" and result["net_harm_avoided"] is None,"Reference packet has no measured local benefit inputs")
    return {"same_design":True,"same_geometry_and_mass":True,"same_structural_capture_source":True,
      "same_continuous_sweep":True,"operation_consumes_actual_reports":True,"outcomes_bound_to_same_sites":True,"physical_opening_remains_unresolved":True}


def refresh():
    commands=[("assembly.py","--output",str(ROOT/"assembly-report.json")),("motion.py",),("site.py",),
              ("evidence.py","--write"),("operation.py","--write"),("outcomes.py","--write")]
    for script,*args in commands:
        print("Rebuilding "+script,flush=True)
        subprocess.run([sys.executable,"-B",str(ROOT/script),*args],check=True,cwd=ROOT)


def build():
    design=load_design(); reports={key:read(ROOT/(key+"-report.json")) for key in COMPONENTS}
    checks=validate(design,reports)
    a,m,s,e,o,h=(reports[k] for k in COMPONENTS)
    source_paths=[ROOT/"design.json",*sorted(ROOT.glob("*.py")),ROOT.parent/"v2/truss.py",ROOT.parent/"v2/engineering.py",
      ROOT.parent/"v2/atlas_adapter.py",ROOT.parent/"v2/site_geometry.py",ROOT.parent/"v2/data/atlas-subset.json",
      ROOT.parent/"v2/oregon-design-basis.md",ROOT.parent/"v2/oregon-sources.json"]
    source_paths += sorted((ROOT.parent/"v2/data/hazards").glob("*.json"))+sorted((ROOT.parent/"v2/data/hazards").glob("*.geojson"))
    return {"schema_version":1,"format":"crest-root-tillamook-integrated/3","design_sha256":design_digest(design),
      "physical_approval":False,"integration_checks":checks,
      "report_content_sha256":{key:design_digest(value) for key,value in reports.items()},
      "report_file_sha256":{key:file_sha(ROOT/(key+"-report.json")) for key in COMPONENTS},
      "source_file_sha256":{os.path.relpath(path,ROOT).replace('\\','/'):file_sha(path) for path in source_paths},
      "summary":{"span_m":a["assembly"]["dimensions"]["span_m"],"physical_nodes":a["assembly"]["physical_node_count"],
        "physical_members":a["assembly"]["physical_member_count"],"moving_mass_kg":a["assembly"]["mass_ledger"]["bridge_mass_kg"],
        "walking_surface_elevation_m":a["assembly"]["dimensions"]["installed_walking_surface_elevation_m"],
        "structural_case_count":len(a["cases"]),"continuous_support_paths":len(m["support_certificates"]),
        "assembled_base_wrench_cases":len(m["assembled_support"]["case_support"]),
        "continuous_release_support_paths":len(m["assembled_support"]["release_support"]),
        "continuous_sweep_boxes":len(m["swept_boxes_model"]),"site":s["summary"],
        "synthetic_evidence_cases":e["case_count"],"physical_measurements":e["physical_measurement_count"],
        "model_service_reachable":"service" in o["reachable_modes"],"public_opening_allowed":False,
        "measured_deployment_time":None,"estimated_real_world_casualty_reduction":None},
      "reports":reports,
      "scope":["One declared reference geometry and mass ledger; no hardware selection or capacity approval",
        "Lift-slew-short-slide is a prescribed kinematic variant, not an engineered telescoping mechanism",
        "Exact mechanics and source GIS use different arithmetic and evidence standards",
        "Synthetic evidence fixtures establish conditional software behavior, not actual lock or support state",
        "No actual public-benefit, deployment-time or comparative PARI/GP performance claim is made"]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--refresh",action="store_true")
    args=parser.parse_args()
    if args.refresh: refresh()
    result=build(); (ROOT/"integrated-report.json").write_text(json.dumps(result,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"integration_checks":result["integration_checks"],"summary":result["summary"]},indent=2))
