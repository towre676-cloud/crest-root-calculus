"""Exact, deliberately bounded engineering screens; not a bridge design code.

Coordinate convention: x/y horizontal, z upward. A downward force W at (x,y)
contributes (W*x,W*y) to the vertical resultant's first moments. A horizontal
force (Fx,Fy) at height h contributes (Fx*h,Fy*h). These first moments define
the centre of pressure on a horizontal support plane; torsion/friction and
individual pad stiffness/capacity are outside this model.
"""
from fractions import Fraction as F
from itertools import combinations


def exact(value):
    if isinstance(value, bool) or isinstance(value, float):
        raise TypeError("Use integers or exact decimal/fraction strings, never floats/bools")
    if not isinstance(value, (int, str, F)):
        raise TypeError("Unsupported exact scalar")
    return F(value)


def enc(value):
    if isinstance(value, F):
        return str(value)
    if isinstance(value, dict):
        return {k: enc(v) for k, v in value.items()}
    if isinstance(value, (list, tuple)):
        return [enc(v) for v in value]
    return value


def cross(a, b):
    return a[0] * b[1] - a[1] * b[0]


def sub(a, b):
    return tuple(x-y for x, y in zip(a, b))


def hull(points):
    points = sorted(set(tuple(map(exact, p)) for p in points))
    if any(len(p) != 2 for p in points):
        raise ValueError("Supports must be horizontal 2D points")
    def half(rows):
        out = []
        for p in rows:
            while len(out) >= 2 and cross(sub(out[-1], out[-2]), sub(p, out[-1])) <= 0:
                out.pop()
            out.append(p)
        return out
    result = half(points)[:-1] + half(points[::-1])[:-1]
    if len(result) < 3:
        raise ValueError("At least three non-collinear support locations are required")
    return result


def support_polygon(supports, downward_load, first_moment_x, first_moment_y):
    """Exact existence of nonnegative vertical reactions at coplanar pads.

    An interior resultant does not establish the actual four-pad reactions.
    A witness using at most three pads demonstrates feasibility only.
    """
    polygon = hull(supports)
    W, Qx, Qy = map(exact, (downward_load, first_moment_x, first_moment_y))
    if W <= 0:
        raise ValueError("Total downward load must be positive")
    c = Qx/W, Qy/W
    margins = [cross(sub(b, a), (Qx-W*a[0], Qy-W*a[1]))
               for a, b in zip(polygon, polygon[1:]+polygon[:1])]
    status = "outside" if min(margins) < 0 else "boundary" if min(margins) == 0 else "inside"
    witness = None
    if status != "outside":
        for a, b, d in combinations(polygon, 3):
            area = cross(sub(b, a), sub(d, a))
            if not area:
                continue
            wb = cross(sub(c, a), sub(d, a))/area
            wd = cross(sub(b, a), sub(c, a))/area
            wa = 1-wb-wd
            if min(wa, wb, wd) >= 0:
                witness = [{"point": p, "reaction": W*w} for p, w in zip((a,b,d), (wa,wb,wd))]
                break
        assert witness is not None
    return {"status": status, "centre_of_pressure": c, "hull": polygon,
            "edge_margins_force_length_squared": margins,
            "feasible_reaction_witness": witness,
            "meaning": "Existence of vertical compressive reactions only; not actual pad forces or capacity"}


def quadratic(coefficients):
    if not isinstance(coefficients, (list, tuple)) or not 1 <= len(coefficients) <= 3:
        raise ValueError("A polynomial needs one to three coefficients, ascending powers")
    return tuple(map(exact, coefficients)) + (F(0),)*(3-len(coefficients))


def poly_value(p, s):
    return p[0]+s*(p[1]+s*p[2])


def minimum_quadratic(coefficients, interval=(0,1)):
    p = quadratic(coefficients)
    lo, hi = map(exact, interval)
    if lo >= hi:
        raise ValueError("Nonempty ordered deployment interval required")
    candidates = [lo, hi]
    if p[2] > 0:
        stationary = -p[1]/(2*p[2])
        if lo < stationary < hi:
            candidates.append(stationary)
    value, argmin = min((poly_value(p,s),s) for s in candidates)
    return {"minimum": value, "at": argmin, "tested_exact_parameters": sorted(candidates)}


def support_path(case):
    """Certify every parameter of a declared quadratic load-position path.

    Constant positive weights, quadratic positions and constant horizontal
    forces with quadratic heights. No hidden interpolation of sampled motion.
    This is a quasi-static path; accelerations and mode transitions are absent.
    """
    polygon = hull(case["supports"])
    if not case.get("weights"):
        raise ValueError("At least one constant weight is required")
    W, Qx, Qy = F(0), [F(0)]*3, [F(0)]*3
    for load in case["weights"]:
        w = exact(load["force"])
        if w <= 0:
            raise ValueError("Weights must be positive")
        W += w
        for target, key in ((Qx,"x"),(Qy,"y")):
            for i, v in enumerate(quadratic(load[key])):
                target[i] += w*v
    for load in case.get("horizontal_forces", []):
        h = quadratic(load["height"])
        if minimum_quadratic(h)["minimum"] < 0:
            raise ValueError("Height must remain above the support plane")
        for target, key in ((Qx,"x"),(Qy,"y")):
            f = exact(load[key])
            for i,v in enumerate(h):
                target[i] += f*v
    edges = []
    for a,b in zip(polygon,polygon[1:]+polygon[:1]):
        e = sub(b,a)
        p = [e[0]*Qy[i]-e[1]*Qx[i] for i in range(3)]
        p[0] += W*(e[1]*a[0]-e[0]*a[1])
        edges.append({"from":a,"to":b,"polynomial":p,**minimum_quadratic(p)})
    margin = min(e["minimum"] for e in edges)
    return {"status": "outside" if margin<0 else "boundary" if margin==0 else "inside",
            "domain":[F(0),F(1)],"total_downward_load":W,
            "first_moment_x_polynomial":Qx,"first_moment_y_polynomial":Qy,
            "edges":edges,
            "meaning":"All s in [0,1] for the declared quasi-static polynomial model only"}


def join_capture(case):
    """Finite rigid transform of every connector; exact unit quaternion.

    A box is a declared geometric capture envelope, not a pin/lock capacity.
    Quaternion order [w,x,y,z]; translation and connector points share units.
    """
    q = list(map(exact,case["relative_rotation_quaternion"]))
    t = list(map(exact,case["relative_translation"]))
    limits = list(map(exact,case["capture_half_widths"]))
    if len(q)!=4 or len(t)!=3 or len(limits)!=3 or any(v<0 for v in limits):
        raise ValueError("Quaternion, translation and nonnegative capture dimensions required")
    if sum(v*v for v in q)!=1:
        raise ValueError("Quaternion must be exactly unit length; no silent normalization")
    w,x,y,z=q
    R=[[1-2*(y*y+z*z),2*(x*y-z*w),2*(x*z+y*w)],
       [2*(x*y+z*w),1-2*(x*x+z*z),2*(y*z-x*w)],
       [2*(x*z-y*w),2*(y*z+x*w),1-2*(x*x+y*y)]]
    rows=[]
    ids=set()
    if not case.get("connectors"):
        raise ValueError("At least one connector required")
    for con in case["connectors"]:
        if con["id"] in ids:
            raise ValueError("Connector identifiers must be unique")
        ids.add(con["id"])
        p=list(map(exact,con["point"]))
        if len(p)!=3:
            raise ValueError("Connector point requires three coordinates")
        delta=[t[i]+sum(R[i][j]*p[j] for j in range(3))-p[i] for i in range(3)]
        margins=[limits[i]-abs(delta[i]) for i in range(3)]
        rows.append({"id":con["id"],"displacement":delta,"margins":margins})
    worst=min(v for row in rows for v in row["margins"])
    return {"status":"outside" if worst<0 else "boundary" if worst==0 else "inside",
            "rotation_matrix":R,"connectors":rows,
            "meaning":"Geometric capture only; engagement, lock verification and load capacity remain separate"}


def access_envelope(clearance, overlay, underside_depth, clearance_reserve,
                    lower_approach_offset, slope, max_run_rise="0.762",
                    landing_length="1.524", outside_lane_width="1.8"):
    """Concept ramp dimensions from explicitly declared inputs, exact in metres.

    lower_approach_offset is road high point minus approach elevation. Negative
    offsets may reduce rise. No surface accessibility or approval is inferred.
    """
    c,o,d,r,h,s,m,l,b=map(exact,(clearance,overlay,underside_depth,clearance_reserve,
                               lower_approach_offset,slope,max_run_rise,landing_length,outside_lane_width))
    if min(c,s,m,l,b)<=0 or min(o,d,r)<0:
        raise ValueError("Positive geometric limits and nonnegative allowances required")
    rise=c+o+d+r+h
    if rise<=0:
        raise ValueError("Positive approach rise required")
    runs=(rise/m).__ceil__()
    run=rise/s
    longest=min(m,rise)/s
    return {"rise":rise,"sloped_horizontal_run":run,"runs":runs,
            "intermediate_landings":runs-1,
            "straight_envelope_length":run+(runs+1)*l,
            "straight_envelope_width":b,
            "unstacked_switchback_envelope_length":longest+2*l,
            "unstacked_switchback_envelope_width":runs*b,
            "meaning":"Bounding footprints including end landings; no structural details, headroom or route certification"}
