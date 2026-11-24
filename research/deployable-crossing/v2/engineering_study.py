#!/usr/bin/env python3
"""Reproducible discrete support/path and finite connector-capture challenges."""
from fractions import Fraction as F
from itertools import product
import json
from pathlib import Path

from engineering import enc, support_path, join_capture, access_envelope

ROOT=Path(__file__).parent


def deployment_case(reach, half_across, half_along, ballast, wind_x, wind_y):
    """Whole half-unit weights, unlike the separate per-plane truss examples.

    Constant deck mass at an affine centre of gravity is a prescribed motion
    model, not a validated telescoping mechanism. All load inputs are hypotheses.
    """
    weights=[{"force":100,"x":[0],"y":[0]},
             {"force":30,"x":[0],"y":[0]},
             {"force":4*reach,"x":[-1,1+reach/2],"y":[0]}]
    if ballast: weights.append({"force":ballast,"x":["-5/2"],"y":[0]})
    return {"supports":[[-half_across,-half_along],[half_across,-half_along],
                        [half_across,half_along],[-half_across,half_along]],
            "weights":weights,
            "horizontal_forces":[{"x":wind_x,"y":wind_y,"height":[6]}]}


def run():
    candidates=[]
    for reach,a,b,ballast in product((F(19,2),F(23,2),F(27,2)),
                                    (F(3,2),F(2),F(5,2)),(F(3,2),F(5,2)),range(0,101,10)):
        corners=[]
        for wx,wy in product((-10,10),(-6,6)):
            case=deployment_case(reach,a,b,ballast,wx,wy)
            result=support_path(case)
            corners.append({"wind_x_kN":wx,"wind_y_kN":wy,"result":result})
        feasible=all(c["result"]["status"]=="inside" for c in corners)
        candidates.append({"half_reach_m":reach,"support_half_width_across_m":a,
                           "support_half_length_along_m":b,"ballast_weight_kN":ballast,
                           "support_rectangle_area_m2":4*a*b,
                           "strict_vertical_contact_feasibility":feasible,"wind_corners":corners})
    front=[]
    for reach in (F(19,2),F(23,2),F(27,2)):
        feasible=[r for r in candidates if r["half_reach_m"]==reach and r["strict_vertical_contact_feasibility"]]
        for row in feasible:
            if not any(other["ballast_weight_kN"]<=row["ballast_weight_kN"]
                       and other["support_rectangle_area_m2"]<=row["support_rectangle_area_m2"]
                       and (other["ballast_weight_kN"]<row["ballast_weight_kN"]
                            or other["support_rectangle_area_m2"]<row["support_rectangle_area_m2"])
                       for other in feasible):
                front.append({k:v for k,v in row.items() if k!="wind_corners"})
    connectors=[{"id":f"{y}:{z}","point":[0,y,z]} for y,z in product((-1,1),(-2,2))]
    aligned={"relative_rotation_quaternion":[1,0,0,0],"relative_translation":["1/100",0,0],
             "capture_half_widths":["1/40","1/200","1/200"],"connectors":connectors}
    rotated=dict(aligned,relative_rotation_quaternion=["39999/40001",0,"400/40001",0])
    missed_midpath={"supports":[[-1,-1],[1,-1],[1,1],[-1,1]],
                    "weights":[{"force":1,"x":["3/4",2,-2],"y":[0]}]}
    report={"format":"crest-root-engineering-challenge/2",
            "units":{"length":"m","force":"kN"},
            "candidate_count":len(candidates),
            "strict_feasible_count":sum(r["strict_vertical_contact_feasibility"] for r in candidates),
            "parameter_domain":"Every s in [0,1]; horizontal forces Fx in [-10,10], Fy in [-6,6] kN at 6m",
            "wind_enclosure_reason":"Each edge margin is affine in the two independent wind forces, so minima over their rectangle occur at its four corners. Each corner is minimized analytically over the whole quadratic path.",
            "assumptions":{"chassis_weight_kN":100,"mast_weight_kN":30,
               "entire_half_deck_weight_kN_per_m_reach":4,"ballast_position_m":["-5/2",0],
               "deck_CG_x_path":"-1 + (1 + half_reach/2) * s", "public_loading_during_extension":False,
               "span_link":"Half-reaches correspond to the GIS scenario base spans19/23/27m, not measured roads"},
            "not_resolved":["Actual weights and CG trajectories","Wind area, drag, gusts and operating threshold",
              "Dynamic inertia and actuator failures","Footpad stiffness, contact area and bearing",
              "Horizontal friction, sliding, torsion and anchorage","Occupied/service loads",
              "Mechanical base/mast strength","Transportation weight and axle limits"],
            "pareto_front_under_declared_model":front,"candidates":candidates,
            "endpoint_only_counterexample":{"input":missed_midpath,"result":support_path(missed_midpath),
                "endpoint_CoP_x":"3/4","interior_CoP_x_at_half":"5/4"},
            "join_examples":[{"id":"aligned-box-capture","input":aligned,"result":join_capture(aligned)},
                             {"id":"pitch-causes-connector-miss","input":rotated,"result":join_capture(rotated)}],
            "flat_ground_access_example":access_envelope("5.2832","0.1","0.6","0.1",0,"3/40"),
            "scope":"Research screens only. No computed capacity, safe wind speed, approved configuration or machine interlock."}
    return enc(report)


if __name__=="__main__":
    report=run()
    (ROOT/"engineering-report.json").write_text(json.dumps(report,indent=2)+"\n")
    print(json.dumps({k:report[k] for k in ("candidate_count","strict_feasible_count","pareto_front_under_declared_model","flat_ground_access_example")},indent=2))
