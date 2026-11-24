"""Place the shared assembly, complete motion envelopes and access hypotheses.

Mapped/GIS arithmetic is binary64, separate from exact mechanics. A swept-box
intersection is a possible obstruction, not a demonstrated physical impact.
No unobserved curb, road width, ground capacity or permission is invented.
"""
import argparse
from collections import Counter
from fractions import Fraction as F
import hashlib
from itertools import product
import json
import math
from pathlib import Path
import sys

from assembly import build_assembly, design_digest, load_design
from motion import sha

HERE=Path(__file__).resolve().parent
V2=HERE.parent/"v2"
sys.path.insert(0,str(V2))
from atlas_adapter import load_subset, sample_terrain
from site_geometry import bbox, overlap, polygon_hit, line_hit, rectangle, polygon_samples

VEHICLE_KINDS={"trunk","primary","secondary","tertiary","residential","unclassified","service","track"}
ARCHITECTURES=("straight-ramp","switchback-ramp","elevator-and-stair-allowance")


def access_geometry(rise):
    rise=F(rise)
    if rise<=0: raise ValueError("A positive approach rise is required")
    slope=F(3,40); max_rise=F(381,500); landing=F(381,250); lane=F(9,5)
    runs=math.ceil(rise/max_rise); run=rise/slope
    return {"rise_m":str(rise),"running_slope":str(slope),"horizontal_ramp_run_m":str(run),"run_count":runs,
      "landing_count":runs+1,"straight_length_m":str(run+(runs+1)*landing),"straight_width_m":str(lane),
      "switchback_length_m":str(run/runs+2*landing),"switchback_width_m":str(runs*lane),
      "scope":"Same v2 planning slope/landing hypotheses; unstacked geometry, no approved route or elevator design"}


def scene_mapping(site,offset,span):
    t=site["road_tangent"]; n=[t[1],-t[0]]
    origin=[site["point"][i]+offset*t[i] for i in (0,1)]
    return {"origin_at_span_midpoint":origin,"cross_road_unit":n,"along_road_unit":t,"span_m":float(span),"model_midpoint_x_m":float(span)/2}


def project(mapping,x,y):
    x=float(F(x))-mapping["model_midpoint_x_m"]; y=float(F(y))
    return [mapping["origin_at_span_midpoint"][i]+x*mapping["cross_road_unit"][i]+y*mapping["along_road_unit"][i] for i in (0,1)]


def model_rectangle(mapping,x0,x1,y0,y1):
    return [project(mapping,x,y) for x,y in ((x0,y0),(x1,y0),(x1,y1),(x0,y1))]


def _index(data):
    return {key:[(row,bbox([p for ring in row["rings"] for p in ring])) for row in data[key]] for key in ("buildings","parcels")}


def _hits(data,index,footprints,road_id):
    rings=[f["ring"] for f in footprints]; box=bbox([p for r in rings for p in r])
    out={"building_ids":[],"building_ids_at_3m_sensitivity":[],"vehicle_centerline_ids":[],"other_vehicle_centerline_ids":[],"parcel_ids":[]}
    for feature,bounds in index["buildings"]:
        if overlap(box,bounds,3):
            if any(polygon_hit([ring],feature["rings"]) for ring in rings): out["building_ids"].append(feature["id"])
            if any(polygon_hit([ring],feature["rings"],3) for ring in rings): out["building_ids_at_3m_sensitivity"].append(feature["id"])
    for road in data["streets"]:
        if road["kind"] in VEHICLE_KINDS and any(line_hit(road["p"],[ring]) for ring in rings):
            out["vehicle_centerline_ids"].append(road["id"])
            if road["id"]!=road_id: out["other_vehicle_centerline_ids"].append(road["id"])
    for feature,bounds in index["parcels"]:
        if overlap(box,bounds) and any(polygon_hit([ring],feature["rings"]) for ring in rings): out["parcel_ids"].append(feature["id"])
    return {key:sorted(set(value)) for key,value in out.items()}


def _hazard_inputs():
    folder=V2/"data/hazards"; raw=(folder/"manifest.json").read_bytes(); manifest=json.loads(raw)
    layers=[]; hashes={"manifest.json":hashlib.sha256(raw).hexdigest()}
    if manifest["schema_version"]!=2 or manifest["study_envelope"]["srid"]!=26910: raise ValueError("Unsupported hazard source")
    for entry in manifest["layers"]:
        path=folder/entry["clipped_file"]; blob=path.read_bytes()
        if hashlib.sha256(blob).hexdigest()!=entry["clipped_sha256"]: raise ValueError("Altered hazard source")
        collection=json.loads(blob)
        if len(collection["features"])!=entry["clipped_feature_count"]: raise ValueError("Hazard feature count mismatch")
        summary_path=folder/entry["raw_summary"]; summary=json.loads(summary_path.read_text(encoding="utf-8"))
        if summary["query_srid"]!=26910: raise ValueError("Wrong source frame")
        clip=manifest["processing"]["clip_bounds"]; query=summary["query_envelope"]
        available=[max(clip[0],query[0]),max(clip[1],query[1]),min(clip[2],query[2]),min(clip[3],query[3])]
        features=[]
        for feature in collection["features"]:
            geometry=feature["geometry"]
            parts=[geometry["coordinates"]] if geometry["type"]=="Polygon" else geometry["coordinates"] if geometry["type"]=="MultiPolygon" else None
            if parts is None: raise ValueError("Expected polygon hazard")
            prop=feature["properties"]; key=entry["key"]
            if key=="fema_nfhl_28": label=" | ".join(str(prop[k]) for k in ("FLD_ZONE","ZONE_SUBTY") if prop.get(k) is not None)
            elif key=="dogami_tsunami_12": label=str(prop.get("Tsunami_Region","unlabelled"))
            elif key=="dogami_tsunami_78": label="XXL1 modeled inundation extent"
            elif key=="dogami_ls_coverage_0": label=str(prop.get("REF_ID_COD","unlabelled"))
            elif key=="dogami_ls_deep_1": label=str(prop.get("Suscept","unlabelled"))
            else: label="mapped landslide deposit"
            features.append({"id":str(prop[entry["object_id_field"]]),"label":label,"parts":parts})
        layers.append({"id":entry["key"],"features":features,"available_bounds":available,
          "source_complete":entry["clip_status"]=="complete" and entry["query_status"] in ("features_returned","no_features_returned"),"endpoint":entry["endpoint"]})
        hashes[entry["clipped_file"]]=entry["clipped_sha256"]
        hashes[entry["raw_summary"]]=hashlib.sha256(summary_path.read_bytes()).hexdigest()
    return layers,hashes


def hazard_hits(layers,footprints,origin):
    rings=[[[point[0]+origin[0],point[1]+origin[1]] for point in f["ring"]] for f in footprints]
    out=[]
    for layer in layers:
        bounds=layer["available_bounds"]
        covered=all(bounds[0]<=p[0]<=bounds[2] and bounds[1]<=p[1]<=bounds[3] for ring in rings for p in ring)
        features=[]
        if covered and layer["source_complete"]:
            for feature in layer["features"]:
                if any(polygon_hit([ring],part) for ring in rings for part in feature["parts"]):
                    features.append({"id":feature["id"],"label":feature["label"]})
        out.append({"layer":layer["id"],"all_footprints_within_query":covered,
          "status":"intersections_evaluated" if covered and layer["source_complete"] else "unknown_source_coverage",
          "features":features})
    return out


def _access_footprints(design,model,mapping,architecture,directions,rises):
    base=design["base_interfaces"]; pad={k:F(v) for k,v in base["pad_rectangle_m"].items()}
    pivot=list(map(F,base["pivot_position_in_base_m"])); span=F(model["dimensions"]["span_m"])
    out=[]
    for side,x,direction,rise in zip(("L","R"),(F(0),span),directions,rises):
        bx=x-pivot[0]; by=-pivot[1]
        out.append({"id":"base-"+side,"kind":"base","side":side,"ring":model_rectangle(mapping,bx+pad["x_min"],bx+pad["x_max"],by+pad["y_min"],by+pad["y_max"])})
        access=access_geometry(rise)
        if architecture=="straight-ramp": length=F(access["straight_length_m"]); width=F(access["straight_width_m"])
        elif architecture=="switchback-ramp": length=F(access["switchback_length_m"]); width=F(access["switchback_width_m"])
        else: length=F(20); width=F(6)
        edge=by+(pad["y_max"] if direction>0 else pad["y_min"])
        y0,y1=sorted((edge,edge+direction*length))
        # The unstacked ramp occupies the outer roadside strip.
        if architecture=="switchback-ramp":
            x0,x1=(bx-width+F(9,10),bx+F(9,10)) if side=="L" else (bx-F(9,10),bx+width-F(9,10))
        else: x0,x1=bx-width/2,bx+width/2
        out.append({"id":"access-"+side,"kind":"access","side":side,"ring":model_rectangle(mapping,x0,x1,y0,y1),"access":access})
    return out


def run_site(design,model,motion_report):
    if model["design_sha256"]!=design_digest(design) or motion_report["design_sha256"]!=design_digest(design) or motion_report["geometry_and_mass_source_sha256"]!=sha(model):
        raise ValueError("Site inputs disagree about the shared design")
    data=load_subset(V2/"data/atlas-subset.json"); index=_index(data); layers,hazard_hashes=_hazard_inputs()
    span=F(model["dimensions"]["span_m"]); walk=F(model["dimensions"]["installed_walking_surface_elevation_m"])
    pad=design["base_interfaces"]["pad_rectangle_m"]
    half_base=(F(pad["x_max"])-F(pad["x_min"]))/2
    road_width=span-2*(F(3)+half_base)
    if road_width<=0: raise ValueError("No positive assumed road envelope fits the declared setback")
    placements=[]; scenarios=[]
    for site,offset in product(data["candidates"],(-30,0,30)):
        mapping=scene_mapping(site,offset,span)
        sweep=[{"id":box["id"],"kind":"sweep_bound","stage_id":box["stage_id"],"side":box["side"],
                "ring":model_rectangle(mapping,box["minimum_m"][0],box["maximum_m"][0],box["minimum_m"][1],box["maximum_m"][1])} for box in motion_report["swept_boxes_model"]]
        placement_id=f'{site["id"]}-o{offset}'
        placement={"id":placement_id,"site_id":site["id"],"label":site["label"],"station_offset_m":offset,"scene_mapping":mapping,"sweep_footprints":sweep,
          "sweep_possible_hits":_hits(data,index,sweep,site["road_way_id"]),"sweep_hazards":hazard_hits(layers,sweep,data["origin"]),
          "meaning":"Intersections with continuous conservative sweep bounds; a hit is a possible obstruction, not proof of collision"}
        placements.append(placement)
        road_samples=[sample_terrain(data,*project(mapping,span/2+road_width*F(k,16)-road_width/2,0)) for k in range(17)]
        road_high=max(road_samples)
        for architecture,directions in product(ARCHITECTURES,product((-1,1),repeat=2)):
            rises=[walk+F(str(road_high-sample_terrain(data,*project(mapping,x,-5)))) for x in (0,span)]
            feet=_access_footprints(design,model,mapping,architecture,directions,rises)
            for _ in range(2):
                for i,side in enumerate(("L","R")):
                    ring=next(f["ring"] for f in feet if f["id"]=="access-"+side)
                    low=min(sample_terrain(data,*point) for point in polygon_samples(ring,4))
                    rises[i]=max(rises[i],walk+F(str(road_high-low)))
                feet=_access_footprints(design,model,mapping,architecture,directions,rises)
            hits=_hits(data,index,feet,site["road_way_id"])
            scenarios.append({"id":f'{placement_id}-{architecture}-{directions[0]}-{directions[1]}',"placement_id":placement_id,"site_id":site["id"],
              "architecture":architecture,"directions":list(directions),"road_width_hypothesis_m":str(road_width),"walking_height_above_reference_road_m":str(walk),
              "source_road_elevation_sample_range_m":[min(road_samples),road_high],"footprints":feet,"direct_hits":hits,
              "hazards":hazard_hits(layers,feet,data["origin"]),"site_approval":"unresolved"})
    by_placement={p["id"]:p for p in placements}
    def rank(row):
        h=row["direct_hits"]; s=by_placement[row["placement_id"]]["sweep_possible_hits"]
        return len(h["building_ids"]),len(h["vehicle_centerline_ids"]),len(s["building_ids"]),len(h["building_ids_at_3m_sensitivity"]),row["id"]
    selected=min(scenarios,key=rank)
    clear=[r for r in scenarios if not r["direct_hits"]["building_ids"] and not r["direct_hits"]["vehicle_centerline_ids"]]
    all_clear=[r for r in clear if not by_placement[r["placement_id"]]["sweep_possible_hits"]["building_ids"]]
    return {"schema_version":1,"format":"crest-root-shared-site/3","design_sha256":design_digest(design),"physical_approval":False,
      "geometry_and_mass_source_sha256":sha(model),"motion_report_sha256":sha(motion_report),"atlas_content_sha256":data["content_sha256"],
      "hazard_source_sha256":hazard_hashes,"sites":data["candidates"],"placements":placements,"scenarios":scenarios,"selected_scenario_id":selected["id"],
      "flat_ground_access":access_geometry(walk),"summary":{"placement_count":len(placements),"scenario_count":len(scenarios),
        "without_direct_building_or_road_conflict":len(clear),"also_without_sweep_bound_building_intersection":len(all_clear),"approved_sites":0,
        "possible_sweep_building_conflicts":sum(bool(p["sweep_possible_hits"]["building_ids"]) for p in placements)},
      "basis":{"assumed_road_width_m":str(road_width),"span_m":str(span),"outer_pad_to_assumed_road_edge_setback_m":"3",
        "elevator_and_stair_plan_allowance_m":[20,6],"widths_are_authored":True,"sensitivity_m":3,
        "clearance_and_ramp_reference":"Same captured Oregon references as v2/oregon-design-basis.md; applicability remains unresolved"},
      "limitations":["No survey accuracy follows from coordinate encoding or DEM interpolation",
        "Two terrain expansion passes and source samples are not a certified terrain envelope",
        "Site checks use horizontal geometry; no building-height or utility-clearance survey is supplied",
        "Sweep rectangles are conservative enclosures and can produce false positive obstruction flags",
        "Target road sweep is intended overhead movement; traffic control and physical clearance remain separate obligations",
        "Railing/actuator/transport/queue geometry and complete accessible connections are not designed",
        "Hazard data are the captured v2 source snapshot; missing/out-of-query information stays unknown",
        "No ranking establishes safety, permission, constructability, uptake or public benefit"]}


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__); parser.add_argument("--output",type=Path,default=HERE/"site-report.json")
    args=parser.parse_args(); d=load_design(); model=build_assembly(d)
    motion=json.loads((HERE/"motion-report.json").read_text(encoding="utf-8"))
    report=run_site(d,model,motion)
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(report["summary"],indent=2))
