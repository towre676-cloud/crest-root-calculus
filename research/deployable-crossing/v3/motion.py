"""Exact prescribed lift/slew/short-slide paths for the shared v3 assembly.

This is an ideal fixture experiment, not a designed actuator mechanism.
Rational rotation matrices preserve each half's member lengths. Continuous
support signs and separation of the declared centreline prisms are checked
analytically. Interval boxes enclose the complete nominal geometric sweep;
animation frames are never used as a clearance certificate.
"""
from fractions import Fraction as F
from itertools import product
import argparse
import hashlib
import importlib.util
import json
from pathlib import Path

from assembly import build_assembly, design_digest, load_design

HERE = Path(__file__).resolve().parent
_spec = importlib.util.spec_from_file_location("crest_v2_engineering_for_motion", HERE.parent / "v2" / "engineering.py")
_eng = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_eng)
STAGES = ("lift", "slew", "approach")
LABELS = {"lift": "Lift the independently held halves", "slew": "Slew with the roots retreated", "approach": "Translate into nominal capture"}


def encode(x):
    if isinstance(x, F): return str(x)
    if isinstance(x, dict): return {k: encode(v) for k,v in x.items()}
    if isinstance(x, (list,tuple)): return [encode(v) for v in x]
    return x


def sha(x):
    return hashlib.sha256(json.dumps(x, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def p(x):
    if isinstance(x, (tuple,list)):
        if len(x)>3: raise ValueError("Degree above two")
        return tuple(F(v) for v in x)+(F(0),)*(3-len(x))
    return F(x),F(0),F(0)


def add(a,b): return tuple(x+y for x,y in zip(p(a),p(b)))
def scale(a,c): return tuple(x*F(c) for x in p(a))
def value(a,t):
    a=p(a); return a[0]+F(t)*(a[1]+F(t)*a[2])


def rational_transform(design, side, stage):
    """Return degree<=2 numerator polynomials over a positive denominator.

    Model coordinates x run from left root 0 to right root span; y runs
    along the road, with the base centre behind its pivot; z is up.
    In the slew, u=1-t, c=(1-u^2)/(1+u^2), s=2u/(1+u^2).
    Both halves start pointing along -y. Right placement is a rotation,
    never a reflection. The two endpoint transforms join exactly.
    """
    if side not in ("L","R") or stage not in STAGES: raise ValueError("Unknown side or stage")
    span=F(design["geometry"]["span_m"])
    r=F(design["motion"]["retreat_m"])
    base=design["base_interfaces"]
    low,high=map(F,(base["transport_truss_bottom_m"],base["truss_bottom_elevation_m"]))
    zero=p(0); one=p(1)
    if stage=="slew":
        den=p([2,-2,1]); c=p([0,2,-1]); s=p([2,-2])
        R=[[c,s,zero],[scale(s,-1),c,zero],[zero,zero,den]] if side=="L" else [[scale(c,-1),s,zero],[scale(s,-1),scale(c,-1),zero],[zero,zero,den]]
        T=[scale(den,-r if side=="L" else span+r),zero,scale(den,high)]
    else:
        den=one
        if stage=="lift":
            R=[[zero,one,zero],[p(-1),zero,zero],[zero,zero,one]]
            T=[p(-r if side=="L" else span+r),zero,p([low,high-low])]
        else:
            sign=1 if side=="L" else -1
            R=[[p(sign),zero,zero],[zero,p(sign),zero],[zero,zero,one]]
            T=[p([-r,r]) if side=="L" else p([span+r,-r]),zero,p(high)]
    return {"denominator":den,"rotation_numerator":R,"translation_numerator":T}


def transform_at(path,t):
    d=value(path["denominator"],t)
    if d<=0: raise ValueError("Nonpositive transform denominator")
    return [[value(v,t)/d for v in row] for row in path["rotation_numerator"]],[value(v,t)/d for v in path["translation_numerator"]]


def point_numerator(path, point):
    return [add(path["translation_numerator"][i],sum_poly(scale(path["rotation_numerator"][i][j],F(point[j])) for j in range(3))) for i in range(3)]


def sum_poly(items):
    out=p(0)
    for item in items: out=add(out,item)
    return out


def moved_point(path, point, t):
    den=value(path["denominator"],t)
    return [value(v,t)/den for v in point_numerator(path,point)]


def poly_bounds(poly,lo=0,hi=1):
    low=_eng.minimum_quadratic(poly,(lo,hi))["minimum"]
    high=-_eng.minimum_quadratic(scale(poly,-1),(lo,hi))["minimum"]
    return low,high


def rational_bounds(num,den,lo=0,hi=1):
    n0,n1=poly_bounds(num,lo,hi); d0,d1=poly_bounds(den,lo,hi)
    if d0<=0: raise ValueError("Interval denominator crosses zero")
    values=[a/b for a,b in product((n0,n1),(d0,d1))]
    return min(values),max(values)


def support_certificates(design, model):
    base=design["base_interfaces"]; box={k:F(v) for k,v in base["pad_rectangle_m"].items()}
    pads=[(box["x_min"],box["y_min"]),(box["x_max"],box["y_min"]),(box["x_max"],box["y_max"]),(box["x_min"],box["y_max"])]
    pivot=list(map(F,base["pivot_position_in_base_m"]))
    gravity=F(design["gravity_m_s2"])/1000
    fixed=base["fixed_mass_items_per_base"]
    mass_fixed=sum(F(v["mass_kg"]) for v in fixed)
    fixed_q=[sum(F(v["mass_kg"])*F(v["center_in_base_m"][i])*gravity for v in fixed) for i in (0,1)]
    bounds=design["motion"]["horizontal_fixture_force_bounds_kN"]
    force_height=F(design["motion"]["horizontal_fixture_force_height_m"])
    records=[]
    for side,stage,fx,fy in product(("L","R"),STAGES,map(F,bounds["x"]),map(F,bounds["y"])):
        half=model["halves"][side]; path=rational_transform(design,side,stage); den=path["denominator"]
        mass=F(half["mass_kg"]); weight=mass*gravity; total=(mass+mass_fixed)*gravity
        cg=point_numerator(path,half["local_center_of_mass_m"])
        base_x=F(0) if side=="L" else F(design["geometry"]["span_m"])
        # The physical base centre is pivot_position behind the model root.
        Qx=add(scale(den,fixed_q[0]+fx*force_height+weight*(pivot[0]-base_x)),scale(cg[0],weight))
        Qy=add(scale(den,fixed_q[1]+fy*force_height+weight*pivot[1]),scale(cg[1],weight))
        edges=[]
        for a,b in zip(pads,pads[1:]+pads[:1]):
            ex,ey=b[0]-a[0],b[1]-a[1]
            num=add(add(scale(Qy,ex),scale(Qx,-ey)),scale(den,total*(ey*a[0]-ex*a[1])))
            low=_eng.minimum_quadratic(num)
            d0,d1=poly_bounds(den)
            # The sign certificate uses the numerator. This quotient bound is
            # conservative and is not advertised as the exact minimum ratio.
            bound=low["minimum"]/(d1 if low["minimum"]>=0 else d0)
            edges.append({"from":a,"to":b,"signed_margin_numerator":num,**low,"conservative_margin_lower_bound_kN_m2":bound})
        minimum=min(e["minimum"] for e in edges)
        records.append({"side":side,"stage_id":stage,"fixture_force_kN":[fx,fy],"fixture_force_height_m":force_height,
          "half_mass_kg":mass,"fixed_base_mass_kg":mass_fixed,"total_vertical_load_kN":total,"denominator":den,
          "first_moment_numerators_kN_m":[Qx,Qy],"edges":edges,
          "status":"inside" if minimum>0 else "boundary" if minimum==0 else "outside"})
    return records


def separation_certificates(design,model):
    span=F(design["geometry"]["span_m"]); w=F(model["dimensions"]["truss_width_m"]); r=F(design["motion"]["retreat_m"])
    radius=F(design["geometry"]["member_envelope_radius_m"])
    contact=F(design["geometry"]["joint_contact_zone_length_m"])
    # X-projection separation proves separation of the declared full prisms.
    # Zero at approach t=1 is intended coincident ideal seam-plane contact.
    rows=[]
    for stage in STAGES:
        path=rational_transform(design,"L",stage); den=path["denominator"]
        if stage=="lift": num=p(span+2*r-w)
        elif stage=="approach": num=p([2*r,-2*r])
        else:
            c=path["rotation_numerator"][0][0]; s=path["rotation_numerator"][0][1]
            num=add(scale(den,span+2*r),add(scale(c,-span),scale(s,-w)))
        # In slew use u=1-t: 2r +(2span+2r)u^2 -2w*u.
        low=_eng.minimum_quadratic(num)
        c=path["rotation_numerator"][0][0]
        noncontact=add(add(num,scale(c,2*contact)),scale(den,-2*radius))
        outside_low=_eng.minimum_quadratic(noncontact)
        rows.append({"stage_id":stage,"denominator":den,"x_gap_numerator_m":num,**low,
          "status":"separated" if low["minimum"]>0 else "intended_endpoint_contact" if low["minimum"]==0 and stage=="approach" else "overlap",
          "outside_joint_zone_numerator_m":noncontact,"outside_joint_zone_minimum":outside_low,
          "member_envelope_radius_m":radius,"joint_contact_zone_length_per_half_m":contact,
          "outside_joint_zone_separated":outside_low["minimum"]>0,
          "scope":"Centreline prisms touch at final seam; inflated member envelopes are separated outside the declared joint zone. Actual engagement/contact inside that zone and omitted hardware remain unresolved"})
    return rows


def sweep_boxes(design,model):
    rows=[]
    radius=F(design["geometry"]["member_envelope_radius_m"])
    for stage,side in product(STAGES,("L","R")):
        path=rational_transform(design,side,stage)
        n=design["motion"]["yaw_interval_count"] if stage=="slew" else 1
        points=list(model["halves"][side]["local_nodes"].values())
        # Include the declared floor/rail mass locations for a conservative
        # centreline/mass-point envelope; no unprovided hardware is invented.
        points += [v["position_m"] for v in model["halves"][side]["local_mass_items"]]
        length=F(model["dimensions"]["half_length_m"])
        fw=F(model["dimensions"]["floor_width_m"])
        floor_top=F(model["dimensions"]["walking_surface_offset_m"])
        points += [list(v) for v in product((F(0),length),(-fw/2,fw/2),(floor_top,))]
        for k in range(n):
            lo,hi=F(k,n),F(k+1,n)
            all_bounds=[[rational_bounds(num,path["denominator"],lo,hi) for num in point_numerator(path,point)] for point in points]
            minima=[min(b[i][0] for b in all_bounds)-radius for i in range(3)]
            maxima=[max(b[i][1] for b in all_bounds)+radius for i in range(3)]
            rows.append({"id":f"{stage}-{side}-{k}","side":side,"stage_id":stage,"parameter_interval":[lo,hi],"minimum_m":minima,"maximum_m":maxima})
    return rows


def make_frames(design,model,count=9):
    frames=[]
    for stage in STAGES:
        for i in range(count):
            t=F(i,count-1); nodes={}; halves=[]
            for side in ("L","R"):
                half=model["halves"][side]; path=rational_transform(design,side,stage); R,T=transform_at(path,t)
                for key,point in half["local_nodes"].items(): nodes[key]=moved_point(path,point,t)
                halves.append({"id":side,"transform":{"rotation":R,"translation":T},"mass_kg":half["mass_kg"],"center_of_mass_m":moved_point(path,half["local_center_of_mass_m"],t)})
            frames.append({"id":f"{stage}-{i}","stage_id":stage,"label":LABELS[stage],"progress":t,"world_nodes_m":nodes,"halves":halves})
    return frames


def check_endpoints(design):
    return all(transform_at(rational_transform(design,side,a),1)==transform_at(rational_transform(design,side,b),0)
               for side in ("L","R") for a,b in zip(STAGES,STAGES[1:]))


def nominal_capture(design,model,positions=None):
    limits=list(map(F,design["motion"]["nominal_capture_half_widths_m"]))
    rows=[]
    for pair in model["connector_pairs"]:
        # Accept the assembly's explicit left/right key contract only.
        left,right=pair["left"],pair["right"]
        lp=model["halves"]["L"]["local_nodes"][left] if positions is None else positions["L"][left]
        rp=model["halves"]["R"]["local_nodes"][right] if positions is None else positions["R"][right]
        a=moved_point(rational_transform(design,"L","approach"),lp,1)
        b=moved_point(rational_transform(design,"R","approach"),rp,1)
        mismatch=[y-x for x,y in zip(a,b)]
        rows.append({"left":left,"right":right,"nominal_difference_m":mismatch,"margins_m":[limit-abs(v) for limit,v in zip(limits,mismatch)]})
    return {"connectors":rows,"capture_half_widths_m":limits,"within_envelope":all(min(row["margins_m"])>=0 for row in rows),
       "meaning":"Rigid nominal geometry" if positions is None else "Independent-half elastic dead-load endpoint from the bound structural case; no foundation error, joint compliance or physical engagement"}


def assembled_support(design,assembly_report):
    """Transport the complete FE root wrench to the same base pad plane.

    Reaction-on-bridge is R,M at the structural root. The bridge applies
    -R,-M to the base. For vertical pad forces N, sum(x*N)=external My and
    sum(y*N)=-external Mx. Horizontal/yaw restraint is a separate obligation.
    """
    base=design["base_interfaces"]; box={k:F(v) for k,v in base["pad_rectangle_m"].items()}
    pads=[(box["x_min"],box["y_min"]),(box["x_max"],box["y_min"]),(box["x_max"],box["y_max"]),(box["x_min"],box["y_max"])]
    pivot=list(map(F,base["pivot_position_in_base_m"]))
    height=F(base["truss_bottom_elevation_m"])
    g=F(design["gravity_m_s2"])/1000
    fixed=base["fixed_mass_items_per_base"]
    fixed_w=sum(F(v["mass_kg"])*g for v in fixed)
    fixed_q=[sum(F(v["mass_kg"])*F(v["center_in_base_m"][i])*g for v in fixed) for i in (0,1)]
    bounds=design["motion"]["horizontal_fixture_force_bounds_kN"]
    fixture_h=F(design["motion"]["horizontal_fixture_force_height_m"])
    records=[]
    for case in assembly_report["cases"]:
        if case["status"]!="unique_linear_solution": continue
        for side,fx,fy in product(("L","R"),map(F,bounds["x"]),map(F,bounds["y"])):
            source=case["bases"][side]; wrench=source["reaction_on_bridge"]
            R=list(map(F,wrench["force_kN"])); M=list(map(F,wrench["moment_kN_m"]))
            W=fixed_w+R[2]
            qx=fixed_q[0]-M[1]-height*R[0]+pivot[0]*R[2]+fixture_h*fx
            qy=fixed_q[1]+M[0]+pivot[1]*R[2]-height*R[1]+fixture_h*fy
            edges=[]
            for a,b in zip(pads,pads[1:]+pads[:1]):
                ex,ey=b[0]-a[0],b[1]-a[1]
                edges.append({"from":a,"to":b,"signed_margin_kN_m2":ex*(qy-a[1]*W)-ey*(qx-a[0]*W)})
            least=min(v["signed_margin_kN_m2"] for v in edges)
            records.append({"case_id":case["id"],"side":side,"fixture_force_kN":[fx,fy],
              "source_reference_origin_m":source["reference_origin_m"],"source_reaction_on_bridge":wrench,
              "root_position_in_base_m":[pivot[0],pivot[1],height],
              "total_vertical_load_kN":W,"required_first_moments_kN_m":[qx,qy],
              "resultant_position_in_base_m":[qx/W,qy/W] if W>0 else None,"edges":edges,
              "status":"outside" if W<=0 or least<0 else "boundary" if least==0 else "inside"})
    clamped="joined_clamped/dead"; service="joined_service/dead"; midpoint="joined_release_midpoint/dead"
    lookup={(r["case_id"],r["side"],*r["fixture_force_kN"]):r for r in records}
    release=[]
    for side,fx,fy in product(("L","R"),map(F,bounds["x"]),map(F,bounds["y"])):
        a,b,m=[lookup[(case,side,fx,fy)] for case in (clamped,service,midpoint)]
        w_identity=2*m["total_vertical_load_kN"]==a["total_vertical_load_kN"]+b["total_vertical_load_kN"]
        q_identity=all(2*z==x+y for x,y,z in zip(a["required_first_moments_kN_m"],b["required_first_moments_kN_m"],m["required_first_moments_kN_m"]))
        margins=[{"from":x["from"],"to":x["to"],"affine_margin_coefficients_kN_m2":[x["signed_margin_kN_m2"],y["signed_margin_kN_m2"]-x["signed_margin_kN_m2"]],
                  "minimum_kN_m2":min(x["signed_margin_kN_m2"],y["signed_margin_kN_m2"])} for x,y in zip(a["edges"],b["edges"])]
        release.append({"side":side,"fixture_force_kN":[fx,fy],"parameter_domain":[0,1],"source_cases":[clamped,service],
          "actual_midpoint_wrench_matches":w_identity and q_identity,"edges":margins,
          "inside_for_complete_release":min(a["total_vertical_load_kN"],b["total_vertical_load_kN"])>0 and all(e["minimum_kN_m2"]>0 for e in margins) and w_identity and q_identity})
    return {"case_support":records,"release_support":release,
      "service_support_inside":all(r["status"]=="inside" for r in records if r["case_id"].startswith("joined_service/")),
      "release_support_inside":all(r["inside_for_complete_release"] for r in release),
      "all_regular_case_support_inside":all(r["status"]=="inside" for r in records),
      "scope":"Complete FE reaction force and moment at each root, translated to the pad plane, plus fixed base weights and declared external fixture loads. This proves existence of coplanar nonnegative vertical pad forces only; no ground, sliding, torsional restraint or hardware capacity."}


def reference_clearance(design,assembly_report):
    """Linear element endpoint bounds over a hypothetical flat road plane."""
    bottom=F(design["base_interfaces"]["truss_bottom_elevation_m"])
    radius=F(design["geometry"]["member_envelope_radius_m"])
    reference=F("5.2832"); overlay=F("0.1"); reserve=F("0.1")
    rows=[]
    for case in assembly_report["cases"]:
        if case["status"]!="unique_linear_solution": continue
        coords=case["physical_deformed_coordinates_m"]
        least=min(F(point[2]) for point in coords.values())+bottom-radius
        rows.append({"case_id":case["id"],"lowest_declared_member_surface_m":least,
          "margin_above_reference_plus_overlay_m":least-reference-overlay,
          "remaining_nominal_reserve_m":least-reference-overlay,
          "satisfies_reference_plus_overlay":least>=reference+overlay})
    return {"road_plane_m":"0","source_reference_m":reference,"assumed_overlay_m":overlay,"nominal_reserve_m":reserve,
      "source":"../v2/oregon-design-basis.md; captured ODOT HDM 0300 reference, applicability unresolved",
      "cases":rows,"service_reference_clearance_satisfied":all(r["satisfies_reference_plus_overlay"] for r in rows if r["case_id"].startswith("joined_service/")),
      "scope":"Lowest endpoint of every linear axial member, minus declared radius, in the small-displacement model over one flat reference road plane. No real road crown, differential settlement, dynamic deflection, vehicle envelope or as-built survey is certified."}


def run_motion(design, model=None, assembly_report=None):
    model=build_assembly(design) if model is None else model
    if model["design_sha256"]!=design_digest(design): raise ValueError("Assembly/design mismatch")
    if assembly_report is None:
        assembly_report=json.loads((HERE/"assembly-report.json").read_text(encoding="utf-8"))
    if assembly_report["design_sha256"]!=design_digest(design) or sha(assembly_report["assembly"])!=sha(model):
        raise ValueError("Structural capture source/model mismatch")
    support=support_certificates(design,model); separation=separation_certificates(design,model); capture=nominal_capture(design,model)
    loaded_capture=nominal_capture(design,model,assembly_report["capture_reference"]["deformed_local_nodes_by_half"])
    installed_support=assembled_support(design,assembly_report)
    clearance=reference_clearance(design,assembly_report)
    report={"schema_version":1,"format":"crest-root-shared-motion/3","design_sha256":design_digest(design),
      "geometry_and_mass_source_sha256":sha(model),"physical_approval":False,
      "architecture":"Prescribed ideal lift, yaw and short-slide fixtures carrying two unchanged rigid deck-truss halves",
      "coordinate_convention":"x=0..span across road; y along road; z up; base centre y=-pivot_y; exact proper rotations",
      "stages":[{"id":stage,"label":LABELS[stage],"parameter_domain":[0,1]} for stage in STAGES],
      "rational_paths":{side:{stage:encode(rational_transform(design,side,stage)) for stage in STAGES} for side in ("L","R")},
      "support_certificates":support,"separation_certificates":separation,"capture":capture,"loaded_capture":loaded_capture,
      "loaded_capture_source_sha256":sha(assembly_report["capture_reference"]),
      "structural_report_sha256":sha(assembly_report),"assembled_support":installed_support,"reference_clearance":clearance,
      "swept_boxes_model":sweep_boxes(design,model),"frames":make_frames(design,model),
      "checks":{"continuous_support_inside":all(v["status"]=="inside" for v in support),
        "continuous_half_separation":all(v["status"]!="overlap" and v["outside_joint_zone_separated"] for v in separation),
        "nominal_capture_within_envelope":capture["within_envelope"],"loaded_capture_within_envelope":loaded_capture["within_envelope"],"stage_endpoints_match":check_endpoints(design),
        "service_support_inside":installed_support["service_support_inside"],"release_support_inside":installed_support["release_support_inside"],
        "service_reference_clearance_satisfied":clearance["service_reference_clearance_satisfied"]},
      "scope":["The same physical member and true-mass ledger is used for every pose",
        "Support proof covers all path parameters and the four corners of the declared fixture-force box",
        "Fixed-height horizontal fixture loads act on the base interface; they are not a bridge wind model or wind-speed limit",
        "Coplanar compressive pad reaction existence only; no individual pad load, sliding, ground or anchorage capacity",
        "The prescribed fixture is not a derived actuator mechanism, transport qualification or dynamic model",
        "Sweep boxes enclose nodes, floor corners and mass points with a declared member-radius allowance; complete railing and hardware geometry is absent",
        "Approach endpoint allows ideal common seam contact; connection forces, deformation compatibility and engagement are separate",
        "No public operation follows from any computed pass"]}
    return encode(report)


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design",type=Path,default=HERE/"design.json")
    parser.add_argument("--output",type=Path,default=HERE/"motion-report.json")
    args=parser.parse_args(); report=run_motion(load_design(args.design))
    args.output.write_text(json.dumps(report,indent=2,ensure_ascii=False,allow_nan=False)+"\n",encoding="utf-8")
    print(json.dumps({"checks":report["checks"],"paths":len(report["support_certificates"]),"sweep_boxes":len(report["swept_boxes_model"]),"frames":len(report["frames"])}))
