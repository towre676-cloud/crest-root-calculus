#!/usr/bin/env python3
"""Independent analytic/adverse checks of bounded screens and GIS predicates."""
from fractions import Fraction as F
from copy import deepcopy
import json
from pathlib import Path

from engineering import support_polygon,support_path,minimum_quadratic,join_capture,access_envelope,exact
from engineering_study import run
from site_geometry import polygon_hit,line_hit,inside,segment_hit


def verify():
    checks=[]
    def check(name,condition):
        assert condition,name
        checks.append(name)
    pad=[[-3,-2],[3,-2],[3,2],[-3,2]]
    for label,W,Qx,Qy,status in [("interior",10,5,4,"inside"),("edge",10,30,0,"boundary"),
                                 ("outside",10,31,0,"outside"),("lateral",10,0,21,"outside")]:
        r=support_polygon(pad,W,Qx,Qy)
        check("support-"+label,r["status"]==status)
        if r["feasible_reaction_witness"]:
            rows=r["feasible_reaction_witness"]
            check("witness-"+label,sum(v["reaction"] for v in rows)==W
                  and sum(v["reaction"]*v["point"][0] for v in rows)==Qx
                  and sum(v["reaction"]*v["point"][1] for v in rows)==Qy
                  and all(v["reaction"]>=0 for v in rows))
    check("quadratic-interior",minimum_quadratic([1,-4,4])["minimum"]==0)
    check("concave-endpoint",minimum_quadratic([1,0,-2])["minimum"]==-1)
    check("constant",minimum_quadratic([4])["minimum"]==4)
    check("stationary-outside",minimum_quadratic([0,1,1])["minimum"]==0)
    report=run()
    counter=report["endpoint_only_counterexample"]
    result=support_path(counter["input"])
    check("endpoint-only-misses-interior",result["status"]=="outside" and min(e["minimum"] for e in result["edges"])==F(-1,2))
    for row in report["join_examples"]:
        result=join_capture(row["input"])
        check(row["id"],result["status"]==("inside" if row["id"].startswith("aligned") else "outside"))
    q=deepcopy(report["join_examples"][0]["input"])
    q["relative_rotation_quaternion"]=[0,0,0,1]
    r=join_capture(q)
    check("finite-rotation-not-small-angle",r["rotation_matrix"]==[[-1,0,0],[0,-1,0],[0,0,1]])
    access=access_envelope("5.2832","0.1","0.6","0.1",0,"3/40")
    check("Oregon-clearance-conversion",F(17*12+4)*F("0.0254")==F("5.2832"))
    check("access-eight-runs",access["rise"]==F("6.0832") and access["runs"]==8)
    check("access-horizontal-run",access["sloped_horizontal_run"]==F(121664,1500))
    check("access-landings",access["straight_envelope_length"]==F(121664,1500)+9*F("1.524"))
    check("198-discrete-candidates",report["candidate_count"]==198)
    check("pareto-nondominated",all(not any(
        o["half_reach_m"]==r["half_reach_m"] and o["strict_vertical_contact_feasibility"]
        and o["ballast_weight_kN"]<=r["ballast_weight_kN"]
        and F(o["support_rectangle_area_m2"])<=F(r["support_rectangle_area_m2"])
        and (o["ballast_weight_kN"]<r["ballast_weight_kN"] or F(o["support_rectangle_area_m2"])<F(r["support_rectangle_area_m2"]))
        for o in report["candidates"]) for r in report["pareto_front_under_declared_model"]))
    # A hole and a disconnected component test distinguish filled footprints
    # from careless bounding-box, exterior-only, or even-odd-union shortcuts.
    outer=[[0,0],[10,0],[10,10],[0,10]]
    hole=[[3,3],[7,3],[7,7],[3,7]]
    inner=[[4,4],[6,4],[6,6],[4,6]]
    touch=[[10,2],[11,2],[11,3],[10,3]]
    separate=[[12,2],[13,2],[13,3],[12,3]]
    check("point-in-hole",not inside([5,5],[outer,hole]))
    check("polygon-in-hole",not polygon_hit([inner],[outer,hole]))
    check("hole-clearance-sensitivity",polygon_hit([inner],[outer,hole],1))
    check("boundary-contact",polygon_hit([outer],[touch]))
    check("separated-polygons",not polygon_hit([outer],[separate]))
    check("line-through-polygon",line_hit([[-1,1],[11,1]],[outer]))
    check("line-inside-hole",not line_hit([[4,4],[6,6]],[outer,hole]))
    check("width-proxy-buffer",line_hit([[11,0],[11,10]],[outer],1))
    check("collinear-segments",segment_hit([0,0],[2,0],[1,0],[3,0]))
    check("disjoint-collinear",not segment_hit([0,0],[1,0],[2,0],[3,0]))
    invalid=[lambda:exact(0.1),lambda:exact(True),lambda:support_polygon([[0,0],[1,0],[2,0]],1,0,0),
             lambda:support_polygon(pad,0,0,0),lambda:minimum_quadratic([1,2,3,4]),
             lambda:minimum_quadratic([1],(1,0)),
             lambda:join_capture(dict(q,relative_rotation_quaternion=[1,0,0,1])),
             lambda:access_envelope(5,0,0,0,0,0)]
    for i,fn in enumerate(invalid):
        try: fn()
        except (TypeError,ValueError): check(f"invalid-{i}",True)
        else: raise AssertionError(f"invalid-{i} was accepted")
    return {"status":"passed","check_count":len(checks),"checks":checks,
            "scope":"Analytic bounded engineering models and planar geometry controls, no physical validation"}


if __name__=="__main__":
    result=verify()
    (Path(__file__).parent/"verification-engineering.json").write_text(json.dumps(result,indent=2)+"\n")
    print(json.dumps(result,indent=2))
