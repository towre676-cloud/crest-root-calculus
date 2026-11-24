#!/usr/bin/env python3
"""Replay mapped Tillamook placement scenarios, never approve a site.

The input is a pinned Atlas subset. Real geometry chooses the test locations
and identifies possible conflicts. Road widths, vehicles and access layouts
remain explicit design hypotheses. Passing a 2D overlap test is not a finding
of available land, reliable soil, accessibility or structural adequacy.
"""
import argparse
from collections import Counter
from fractions import Fraction
import itertools
import json
import math
from pathlib import Path

from atlas_adapter import load_subset, sample_terrain
from engineering import access_envelope, enc
from site_geometry import bbox, overlap, polygon_hit, line_hit, rectangle, polygon_samples

ROOT=Path(__file__).parent
VEHICLE_KINDS={"trunk","primary","secondary","tertiary","residential","unclassified","service","track"}
ASSUMPTIONS={
    "road_envelope_widths_m":[8,12,16],
    "station_offsets_m":[-30,0,30],
    "road_clearance_m":"5.2832",
    "road_clearance_basis":"ODOT HDM 317.3 pedestrian reference; applicability not approved",
    "future_overlay_m":"0.1",
    "structure_below_walking_surface_m":"0.6",
    "clearance_reserve_m":"0.1",
    "accessible_ramp_running_slope":"3/40",
    "ramp_slope_basis":"7.5% ODOT temporary route planning reference; applicability open",
    "max_ramp_run_rise_m":"0.762",
    "landing_length_m":"1.524",
    "outside_ramp_lane_width_m":"1.8",
    "concept_clear_route_width_m":"1.524",
    "road_edge_to_nearest_support_envelope_m":3,
    "base_and_outrigger_envelope_along_road_m":12,
    "base_and_outrigger_envelope_across_road_m":5,
    "elevator_stair_allowance_along_road_m":20,
    "elevator_stair_allowance_across_road_m":6,
    "building_separation_sensitivity_m":3,
    "sensitivity_basis":"User-independent scenario allowance, not a measured error bound",
    "position_note":"Station offsets are straight offsets along local source secants, not surveyed staking points",
    "access_note":"Bounding layouts only; no stair/elevator design, land rights, headroom or egress approval",
    "clearance_note":"DEM is not a road-surface survey; allowances are assumed, not load/settlement predictions"
}


def _ramp(rise):
    # Decimal conversion records GIS approximation, not false exactness of GIS.
    # Access arithmetic after that conversion is rational.
    return access_envelope(ASSUMPTIONS["road_clearance_m"],ASSUMPTIONS["future_overlay_m"],
        ASSUMPTIONS["structure_below_walking_surface_m"],ASSUMPTIONS["clearance_reserve_m"],
        str(rise),ASSUMPTIONS["accessible_ramp_running_slope"],
        ASSUMPTIONS["max_ramp_run_rise_m"],ASSUMPTIONS["landing_length_m"],
        ASSUMPTIONS["outside_ramp_lane_width_m"])


def _templates(origin,t,n,width,architecture,directions,rises):
    out=[]
    for side,direction,offset in zip((-1,1),directions,rises):
        base=side*(width/2+3+2.5)
        out.append({"id":f"base-{side}","kind":"base","side":side,
                    "ring":rectangle(origin,t,n,-6,6,base-2.5,base+2.5)})
        access=_ramp(offset)
        if architecture=="straight-ramp":
            length=float(access["straight_envelope_length"])
            n0,n1=base-.9,base+.9
        elif architecture=="switchback-ramp":
            length=float(access["unstacked_switchback_envelope_length"])
            wide=float(access["unstacked_switchback_envelope_width"])
            n0,n1=sorted((base-.9*side,base+(wide-.9)*side))
        elif architecture=="elevator-and-stair-allowance":
            length=20
            n0,n1=base-3,base+3
        else:
            raise ValueError("Unknown architecture")
        s0,s1=sorted((6*direction,(6+length)*direction))
        out.append({"id":f"access-{side}","kind":"access","side":side,
                    "ring":rectangle(origin,t,n,s0,s1,n0,n1),"access":enc(access)})
    return out


def _terrain_range(data,ring):
    values=[sample_terrain(data,*p) for p in polygon_samples(ring,3)]
    return {"sample_min_m":min(values),"sample_max_m":max(values),"sample_count":len(values),
            "meaning":"Unblended source terrain samples/interpolants; not a certified elevation envelope"}


def _feature_bounds(row,key="rings"):
    return bbox([p for ring in row[key] for p in ring])


def scenario(data,site,width,offset,architecture,directions):
    t=site["road_tangent"]; n=[t[1],-t[0]]
    origin=[site["point"][i]+offset*t[i] for i in (0,1)]
    road_samples=[sample_terrain(data,origin[0]+u*n[0],origin[1]+u*n[1])
                  for u in [-width/2+i*width/16 for i in range(17)]]
    road_high=max(road_samples)
    base_offsets=[]
    for side in (-1,1):
        centre=[origin[i]+side*(width/2+5.5)*n[i] for i in (0,1)]
        base_offsets.append(road_high-sample_terrain(data,*centre))
    footprints=_templates(origin,t,n,width,architecture,directions,base_offsets)
    # Two monotonic expansion passes use the lowest sampled access ground.
    # This is explicitly a sensitivity layout, not convergence or an extrema proof.
    for _ in range(2):
        for i,side in enumerate((-1,1)):
            ring=next(f["ring"] for f in footprints if f["kind"]=="access" and f["side"]==side)
            base_offsets[i]=max(base_offsets[i],road_high-_terrain_range(data,ring)["sample_min_m"])
        footprints=_templates(origin,t,n,width,architecture,directions,base_offsets)
    rings=[f["ring"] for f in footprints]
    # Every separate footprint is filled independently (not an even-odd union).
    box=bbox([p for r in rings for p in r])
    buildings=[b for b in data["buildings"] if overlap(box,_feature_bounds(b),3)]
    hits={"building_ids":[],"building_ids_at_3m_sensitivity":[],
          "vehicle_centerline_ids":[],"vehicle_width_proxy_ids":[],
          "mapped_footway_ids":[],"parcel_ids":[],"environment_ids":[],"canopy_proxy_ids":[]}
    for b in buildings:
        if any(polygon_hit([r],b["rings"]) for r in rings): hits["building_ids"].append(b["id"])
        if any(polygon_hit([r],b["rings"],3) for r in rings): hits["building_ids_at_3m_sensitivity"].append(b["id"])
    for road in data["streets"]:
        if road["kind"] in VEHICLE_KINDS:
            if any(line_hit(road["p"],[r]) for r in rings): hits["vehicle_centerline_ids"].append(road["id"])
            # The selected road uses the scenario width. Other roads retain a
            # labelled display-width sensitivity, never a measured carriageway.
            road_width=width if road["id"]==site["road_way_id"] else road["width"]
            if any(line_hit(road["p"],[r],road_width/2) for r in rings):
                hits["vehicle_width_proxy_ids"].append(road["id"])
        elif road["kind"] in {"footway","path","steps","pedestrian","cycleway"}:
            if any(line_hit(road["p"],[r]) for r in rings): hits["mapped_footway_ids"].append(road["id"])
    for parcel in data["parcels"]:
        if overlap(box,_feature_bounds(parcel)) and any(polygon_hit([r],parcel["rings"]) for r in rings):
            hits["parcel_ids"].append(parcel["id"])
    for obj in data["environment"]:
        kind=obj["geometry"]
        hit=any(polygon_hit([r],[obj["p"]]) if kind=="area" else
                line_hit([obj["p"]] if kind=="point" else obj["p"],[r]) for r in rings)
        if hit: hits["environment_ids"].append(obj["id"])
    for tree in data["canopy"]:
        # Adapter retains source proxy attributes rather than surveyed trunks.
        p=tree["p"]
        if any(line_hit([p],[r],tree["radius_estimate_m"]) for r in rings): hits["canopy_proxy_ids"].append(tree["id"])
    for f in footprints:
        f["terrain"]=_terrain_range(data,f["ring"])
    conflict=bool(hits["building_ids"] or hits["vehicle_centerline_ids"])
    return {"id":f'{site["id"]}-w{width}-o{offset}-{architecture}-{directions[0]}-{directions[1]}',
        "site_id":site["id"],"width_hypothesis_m":width,"station_offset_m":offset,
        "architecture":architecture,"directions":directions,"origin":origin,
        "base_centre_span_m":width+11,"required_half_reach_m":(width+11)/2,
        "road_source_elevation_sample_range_m":[min(road_samples),max(road_samples)],
        "road_to_approach_elevation_offsets_m":base_offsets,
        "footprints":footprints,"hits":{k:sorted(set(v)) for k,v in hits.items()},
        "mapped_obstacle_screen":"conflict_in_declared_layout" if conflict else "no_tested_building_or_road_centerline_conflict",
        "site_approval":"unresolved",
        "scope":"Static 2D footprint hypotheses; no vehicle approach/swing, utilities completeness or permissions"}


def rank(row):
    h=row["hits"]
    return (len(h["building_ids"]),len(h["vehicle_centerline_ids"]),
            len(h["building_ids_at_3m_sensitivity"]),len(h["vehicle_width_proxy_ids"]),
            len(h["environment_ids"]),len(h["canopy_proxy_ids"]),row["id"])


def run_study(data):
    rows=[]
    architectures=["straight-ramp","switchback-ramp","elevator-and-stair-allowance"]
    for site,width,offset,arch,directions in itertools.product(data["candidates"],
            ASSUMPTIONS["road_envelope_widths_m"],ASSUMPTIONS["station_offsets_m"],
            architectures,list(itertools.product((-1,1),repeat=2))):
        rows.append(scenario(data,site,width,offset,arch,directions))
    sites=[]
    for site in data["candidates"]:
        options=[r for r in rows if r["site_id"]==site["id"]]
        selected=[min((r for r in options if r["architecture"]==arch and r["width_hypothesis_m"]==w),key=rank)["id"]
                  for arch,w in itertools.product(architectures,ASSUMPTIONS["road_envelope_widths_m"])]
        sites.append(dict(site,selected_scenario_ids=selected,
                     count_without_tested_building_or_road_centerline_conflict=sum(not(r["hits"]["building_ids"] or r["hits"]["vehicle_centerline_ids"]) for r in options)))
    return {"format":"crest-root-tillamook-study/2","atlas_content_sha256":data["content_sha256"],
            "assumptions":ASSUMPTIONS,"sites":sites,"scenarios":rows,
            "summary":{"scenario_count":len(rows),"counts_by_screen":dict(Counter(r["mapped_obstacle_screen"] for r in rows)),
                "approved_sites":0,"architecture_count":len(architectures),
                "minimum_assumed_base_span_m":19,"maximum_assumed_base_span_m":27},
            "interpretation":"The ranking orders declared map conflicts only. It is not safety, desirability, constructability or legal feasibility.",
            "open_site_gates":["Road/curb/shoulder survey and actual widths","ROW and property permission",
              "Utility locate and vertical clearance","Geotechnical bearing/settlement/liquefaction",
              "Accessible route design and lift/elevator availability","Vehicle approach, swept paths and traffic control",
              "Current physical changes and mapping completeness","Hazard and evacuation review","Responsible authority design basis"]}


def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas",type=Path,default=ROOT/"data/atlas-subset.json")
    parser.add_argument("--output",type=Path,default=ROOT/"site-report.json")
    args=parser.parse_args()
    report=run_study(load_subset(args.atlas))
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps(report["summary"],indent=2))


if __name__=="__main__": main()
