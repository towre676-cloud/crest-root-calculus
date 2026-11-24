"""Independent exact audit of the v3 prescribed-motion report.

No engineering capacity or operation approval is implied. This checks the
declared rigid nominal geometry, exact polynomial bounds, vertical-contact
existence and linear-model capture/clearance provenance. Default checks do
not reuse motion.py's transform, extrema or wrench formulas. Negative controls
exercise its public computations; --replay also compares a complete rebuild.
"""

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
from itertools import product
import json
from pathlib import Path

import assembly
import motion
import verify_assembly

HERE = Path(__file__).resolve().parent
STAGES = ("lift", "slew", "approach")
CONTROLS = (F(0), F(1, 4), F(1, 2), F(3, 4), F(1))


def require(condition, message):
    if not condition:
        raise ValueError(message)


def vec(values):
    return [F(value) for value in values]


def canonical(value):
    return hashlib.sha256(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode()).hexdigest()


def poly(values):
    values = list(map(F, values))
    while len(values) > 1 and values[-1] == 0:
        values.pop()
    return tuple(values) if values else (F(0),)


def plus(a, b):
    size = max(len(a), len(b))
    return poly([(a[i] if i < len(a) else 0)+(b[i] if i < len(b) else 0) for i in range(size)])


def times(a, b):
    result = [F(0)]*(len(a)+len(b)-1)
    for i, x in enumerate(a):
        for j, y in enumerate(b):
            result[i+j] += x*y
    return poly(result)


def scale(a, value):
    return poly([F(value)*x for x in a])


def evaluate(a, t):
    result = F(0)
    for coefficient in reversed(a):
        result = result*F(t)+coefficient
    return result


def minimum(coefficients, lo=F(0), hi=F(1)):
    """Independent endpoint/stationary-point enumeration, including maxima."""
    coefficients = poly(coefficients)
    require(len(coefficients) <= 3 and lo <= hi, "Expected quadratic on an ordered interval")
    candidates = [F(lo), F(hi)]
    if len(coefficients) == 3 and coefficients[2]:
        stationary = -coefficients[1]/(2*coefficients[2])
        if lo < stationary < hi:
            candidates.append(stationary)
    return min((evaluate(coefficients, point), point) for point in candidates)


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def add_vectors(a, b):
    return [x+y for x, y in zip(a, b)]


def path_parts(path):
    return (poly(path["denominator"]), [[poly(entry) for entry in row] for row in path["rotation_numerator"]],
            [poly(entry) for entry in path["translation_numerator"]])


def pose_from_path(path, t):
    den, rotation, translation = path_parts(path)
    divisor = evaluate(den, t)
    require(divisor > 0, "Nonpositive path denominator")
    return [[evaluate(entry, t)/divisor for entry in row] for row in rotation], [evaluate(entry, t)/divisor for entry in translation]


def analytic_pose(design, side, stage, t):
    """Direct half-angle and affine-coordinate formulas from the specification."""
    t = F(t)
    retreat, span = F(design["motion"]["retreat_m"]), F(design["geometry"]["span_m"])
    low, high = (F(design["base_interfaces"][key]) for key in ("transport_truss_bottom_m", "truss_bottom_elevation_m"))
    zero = F(0)
    if stage == "lift":
        return [[zero, F(1), zero], [F(-1), zero, zero], [zero, zero, F(1)]], [-retreat if side == "L" else span+retreat, zero, low+t*(high-low)]
    if stage == "approach":
        sign = F(1) if side == "L" else F(-1)
        return [[sign, zero, zero], [zero, sign, zero], [zero, zero, F(1)]], [-retreat*(1-t) if side == "L" else span+retreat*(1-t), zero, high]
    require(stage == "slew", "Unknown stage")
    u = 1-t
    c, s = (1-u*u)/(1+u*u), 2*u/(1+u*u)
    sign = F(1) if side == "L" else F(-1)
    return [[sign*c, s, zero], [-s, sign*c, zero], [zero, zero, F(1)]], [-retreat if side == "L" else span+retreat, zero, high]


def apply_pose(pose, point):
    rotation, translation = pose
    return [dot(row, vec(point))+shift for row, shift in zip(rotation, translation)]


def point_polynomials(path, point):
    den, rotation, translation = path_parts(path)
    answer = []
    for row, shift in zip(rotation, translation):
        component = shift
        for entry, coordinate in zip(row, vec(point)):
            component = plus(component, scale(entry, coordinate))
        answer.append(component)
    return answer


def check_paths(report, design, model):
    identities, controls, endpoint_checks = 0, 0, 0
    require(set(report["rational_paths"]) == {"L", "R"}, "Missing half path")
    for side in ("L", "R"):
        require(set(report["rational_paths"][side]) == set(STAGES), "Missing stage path")
        for stage in STAGES:
            path = report["rational_paths"][side][stage]
            den, rotation, translation = path_parts(path)
            require(len(den) <= 3 and all(len(entry) <= 3 for row in rotation for entry in row)
                    and all(len(entry) <= 3 for entry in translation), "Path numerator degree exceeds its exact identity controls")
            require(minimum(den)[0] > 0, "Transform denominator is not positive over the entire stage")
            den2, den3 = times(den, den), times(times(den, den), den)
            for i in range(3):
                for j in range(3):
                    value = (F(0),)
                    for k in range(3):
                        value = plus(value, times(rotation[i][k], rotation[j][k]))
                    require(value == (den2 if i == j else (F(0),)), "Polynomial rotation orthogonality failed")
                    identities += 1
            determinant = (F(0),)
            for columns, sign in (((0, 1, 2), 1), ((1, 2, 0), 1), ((2, 0, 1), 1),
                                  ((0, 2, 1), -1), ((2, 1, 0), -1), ((1, 0, 2), -1)):
                term = times(times(rotation[0][columns[0]], rotation[1][columns[1]]), rotation[2][columns[2]])
                determinant = plus(determinant, scale(term, sign))
            require(determinant == den3, "Transform is not a proper rotation")
            identities += 1
            for t in CONTROLS:
                require(pose_from_path(path, t) == analytic_pose(design, side, stage, t), "Rational path disagrees with independent analytic pose")
                controls += 1
        for before, after in zip(STAGES, STAGES[1:]):
            require(pose_from_path(report["rational_paths"][side][before], 1) == pose_from_path(report["rational_paths"][side][after], 0), "Stage endpoints do not join")
            endpoint_checks += 1
        pose = analytic_pose(design, side, "approach", 1)
        bottom = F(design["base_interfaces"]["truss_bottom_elevation_m"])
        for node, local in model["halves"][side]["local_nodes"].items():
            expected = vec(model["physical_nodes"][node])
            expected[2] += bottom
            require(apply_pose(pose, local) == expected, "Final pose disagrees with installed structural assembly")
            endpoint_checks += 1
    return {"polynomial_rotation_identities": identities, "analytic_pose_controls": controls, "endpoint_checks": endpoint_checks}


def pad_vertices(design):
    b = {key: F(value) for key, value in design["base_interfaces"]["pad_rectangle_m"].items()}
    return [(b["x_min"], b["y_min"]), (b["x_max"], b["y_min"]), (b["x_max"], b["y_max"]), (b["x_min"], b["y_max"])]


def fixture_corners(design):
    bounds = design["motion"]["horizontal_fixture_force_bounds_kN"]
    return list(product(map(F, bounds["x"]), map(F, bounds["y"])))


def direct_deployment_wrench(design, model, side, stage, t, fixture):
    """Sum every true mass point and r cross F directly in the base frame."""
    base = design["base_interfaces"]
    pivot, gravity = vec(base["pivot_position_in_base_m"]), F(design["gravity_m_s2"])/1000
    span = F(design["geometry"]["span_m"])
    base_origin_world = [(-pivot[0] if side == "L" else span-pivot[0]), -pivot[1], F(0)]
    total_force, total_moment = [F(0)]*3, [F(0)]*3
    pose = analytic_pose(design, side, stage, t)
    for item in model["halves"][side]["local_mass_items"]:
        world = apply_pose(pose, item["position_m"])
        point = [a-b for a, b in zip(world, base_origin_world)]
        force = [F(0), F(0), -F(item["mass_kg"])*gravity]
        total_force = add_vectors(total_force, force)
        total_moment = add_vectors(total_moment, cross(point, force))
    for item in base["fixed_mass_items_per_base"]:
        force = [F(0), F(0), -F(item["mass_kg"])*gravity]
        total_force = add_vectors(total_force, force)
        total_moment = add_vectors(total_moment, cross(vec(item["center_in_base_m"]), force))
    horizontal = [*fixture, F(0)]
    total_force = add_vectors(total_force, horizontal)
    total_moment = add_vectors(total_moment, cross([F(0), F(0), F(design["motion"]["horizontal_fixture_force_height_m"])], horizontal))
    return -total_force[2], [total_moment[1], -total_moment[0]]


def check_support(report, design, model):
    rows = report["support_certificates"]
    expected = {(side, stage, fx, fy) for side, stage, (fx, fy) in product(("L", "R"), STAGES, fixture_corners(design))}
    require(len(rows) == len(expected), "Deployment support record count mismatch")
    seen, polynomial_checks, wrench_checks = set(), 0, 0
    pads = pad_vertices(design)
    for row in rows:
        side, stage = row["side"], row["stage_id"]
        fixture = tuple(vec(row["fixture_force_kN"]))
        key = (side, stage, *fixture)
        require(key in expected and key not in seen, "Duplicate or unexpected deployment support record")
        seen.add(key)
        path = report["rational_paths"][side][stage]
        den = poly(row["denominator"])
        require(den == poly(path["denominator"]), "Support uses another path denominator")
        qx, qy = map(poly, row["first_moment_numerators_kN_m"])
        W = F(row["total_vertical_load_kN"])
        require(W > 0, "Nonpositive vertical load")
        require(F(row["fixture_force_height_m"]) == F(design["motion"]["horizontal_fixture_force_height_m"]), "Fixture-force height label mismatch")
        require(F(row["half_mass_kg"]) == F(model["halves"][side]["mass_kg"])
                and F(row["fixed_base_mass_kg"]) == sum((F(item["mass_kg"]) for item in design["base_interfaces"]["fixed_mass_items_per_base"]), F(0)), "Support mass labels mismatch")
        for t in CONTROLS:
            weight, moments = direct_deployment_wrench(design, model, side, stage, t, fixture)
            require(W == weight, "True mass ledger and support weight disagree")
            require([evaluate(qx, t)/evaluate(den, t), evaluate(qy, t)/evaluate(den, t)] == moments, "Support moment polynomial disagrees with direct mass/fixture wrench")
            wrench_checks += 1
        # Both sides of each equality are degree<=2 numerators. Five distinct
        # exact controls therefore also certify their polynomial identity.
        require(len(qx) <= 3 and len(qy) <= 3, "Moment degree exceeds certified controls")
        require(len(row["edges"]) == 4, "Missing support polygon edge")
        least = None
        for edge, start, end in zip(row["edges"], pads, pads[1:]+pads[:1]):
            require(tuple(vec(edge["from"])) == start and tuple(vec(edge["to"])) == end, "Wrong support edge")
            ex, ey = end[0]-start[0], end[1]-start[1]
            expected_poly = plus(plus(scale(qy, ex), scale(qx, -ey)), scale(den, W*(ey*start[0]-ex*start[1])))
            require(poly(edge["signed_margin_numerator"]) == expected_poly, "Support edge numerator sign/formula error")
            low, at = minimum(expected_poly)
            require(F(edge["minimum"]) == low and F(edge["at"]) == at, "Whole-interval quadratic minimum is incorrect")
            dmin, dmax = minimum(den)[0], -minimum(scale(den, -1))[0]
            require(F(edge["conservative_margin_lower_bound_kN_m2"]) == low/(dmax if low >= 0 else dmin), "Rational conservative lower bound is incorrect")
            least = low if least is None else min(least, low)
            polynomial_checks += 1
        status = "inside" if least > 0 else "boundary" if least == 0 else "outside"
        require(row["status"] == status, "No-lift status ignores exact boundary/outside result")
    require(seen == expected, "A force corner or stage is missing")
    return {"records": len(rows), "direct_mass_wrench_controls": wrench_checks, "whole_interval_edge_minima": polynomial_checks,
            "force_box_coverage": "All four corners; each signed edge margin is affine in the two independent fixture force coordinates"}


def prism_points(model, contact=F(0)):
    length = F(model["dimensions"]["half_length_m"])-contact
    width, height = F(model["dimensions"]["truss_width_m"]), F(model["dimensions"]["truss_depth_m"])
    return [list(point) for point in product((F(0), length), (-width/2, width/2), (F(0), height))]


def check_separation(report, design, model):
    rows = report["separation_certificates"]
    require({row["stage_id"] for row in rows} == set(STAGES) and len(rows) == 3, "Incomplete separation stages")
    radius = F(design["geometry"]["member_envelope_radius_m"])
    contact = F(design["geometry"]["joint_contact_zone_length_m"])
    controls = 0
    for row in rows:
        stage, den = row["stage_id"], poly(row["denominator"])
        require(den == poly(report["rational_paths"]["L"][stage]["denominator"]), "Gap denominator mismatch")
        num, outside = poly(row["x_gap_numerator_m"]), poly(row["outside_joint_zone_numerator_m"])
        require(len(num) <= 3 and len(outside) <= 3, "Gap degree unsupported")
        for t in CONTROLS:
            full, trimmed = [], []
            for side in ("L", "R"):
                pose = analytic_pose(design, side, stage, t)
                full.append([apply_pose(pose, point)[0] for point in prism_points(model)])
                trimmed.append([apply_pose(pose, point)[0] for point in prism_points(model, contact)])
            gap = min(full[1])-max(full[0])
            finite_gap = min(trimmed[1])-max(trimmed[0])-2*radius
            require(evaluate(num, t)/evaluate(den, t) == gap, "Gap polynomial disagrees with transformed full prism extrema")
            require(evaluate(outside, t)/evaluate(den, t) == finite_gap, "Reserved-zone finite envelope gap formula is wrong")
            controls += 1
        low, at = minimum(num)
        require(F(row["minimum"]) == low and F(row["at"]) == at, "Incorrect whole-stage centerline gap minimum")
        low_out, at_out = minimum(outside)
        require(F(row["outside_joint_zone_minimum"]["minimum"]) == low_out and F(row["outside_joint_zone_minimum"]["at"]) == at_out, "Incorrect whole-stage finite envelope gap minimum")
        expected_status = "separated" if low > 0 else "intended_endpoint_contact" if low == 0 and stage == "approach" else "overlap"
        require(row["status"] == expected_status and row["outside_joint_zone_separated"] == (low_out > 0), "Separation status mismatch")
        require(F(row["member_envelope_radius_m"]) == radius and F(row["joint_contact_zone_length_per_half_m"]) == contact, "Gap uses different geometry assumptions")
    return {"transformed_prism_controls": controls, "whole_interval_gap_minima": 6,
            "scope": "Declared centerline prisms and radius-inflated trimmed prisms; no contact/capacity claim inside the reserved seam zone"}


def check_sweep(report, design, model):
    rows = report["swept_boxes_model"]
    expected_count = 2*(2+design["motion"]["yaw_interval_count"])
    require(len(rows) == expected_count and len({row["id"] for row in rows}) == len(rows), "Sweep interval inventory mismatch")
    radius = F(design["geometry"]["member_envelope_radius_m"])
    inequalities = 0
    for stage, side in product(STAGES, ("L", "R")):
        group = [row for row in rows if row["stage_id"] == stage and row["side"] == side]
        n = design["motion"]["yaw_interval_count"] if stage == "slew" else 1
        require(len(group) == n, "Wrong stage interval count")
        group.sort(key=lambda row: F(row["parameter_interval"][0]))
        half = model["halves"][side]
        points = list(half["local_nodes"].values())+[item["position_m"] for item in half["local_mass_items"]]
        length, width, ztop = (F(model["dimensions"][key]) for key in ("half_length_m", "floor_width_m", "walking_surface_offset_m"))
        points += [list(point) for point in product((F(0), length), (-width/2, width/2), (ztop,))]
        path = report["rational_paths"][side][stage]
        den = poly(path["denominator"])
        for i, row in enumerate(group):
            lo, hi = vec(row["parameter_interval"])
            require((lo, hi) == (F(i, n), F(i+1, n)), "Sweep intervals leave a gap or overlap unexpectedly")
            lower, upper = vec(row["minimum_m"]), vec(row["maximum_m"])
            require(all(a < b for a, b in zip(lower, upper)), "Empty sweep box")
            require(minimum(den, lo, hi)[0] > 0, "Sweep denominator may vanish")
            for point in points:
                for axis, num in enumerate(point_polynomials(path, point)):
                    low_poly = plus(num, scale(den, -(lower[axis]+radius)))
                    high_poly = plus(scale(den, upper[axis]-radius), scale(num, -1))
                    require(minimum(low_poly, lo, hi)[0] >= 0 and minimum(high_poly, lo, hi)[0] >= 0, "A complete-path point/member-radius envelope escapes its sweep box")
                    inequalities += 2
    return {"interval_boxes": len(rows), "whole_interval_polynomial_enclosure_inequalities": inequalities,
            "basis": "Point-coordinate numerator inequalities over every full interval; linear bar segments and floor prisms follow by convexity, with the declared radius allowance"}


def check_frames(report, design, model):
    frames = report["frames"]
    require(len(frames) == 27 and len({frame["id"] for frame in frames}) == len(frames), "Animation frame inventory changed")
    node_checks, member_checks = 0, 0
    for frame in frames:
        stage, t = frame["stage_id"], F(frame["progress"])
        require(stage in STAGES and 0 <= t <= 1, "Invalid display frame parameter")
        require(set(frame["world_nodes_m"]) == set(model["physical_nodes"]), "Display frame loses physical nodes")
        halves = {half["id"]: half for half in frame["halves"]}
        for side in ("L", "R"):
            pose = analytic_pose(design, side, stage, t)
            shown = halves[side]
            require(([vec(row) for row in shown["transform"]["rotation"]], vec(shown["transform"]["translation"])) == pose, "Display transform inconsistent")
            require(F(shown["mass_kg"]) == F(model["halves"][side]["mass_kg"]), "Display frame changes mass")
            require(vec(shown["center_of_mass_m"]) == apply_pose(pose, model["halves"][side]["local_center_of_mass_m"]), "Display centroid inconsistent")
            for node, point in model["halves"][side]["local_nodes"].items():
                require(vec(frame["world_nodes_m"][node]) == apply_pose(pose, point), "Display node does not follow the certified transform")
                node_checks += 1
        for bar in model["physical_members"]:
            a, b = vec(frame["world_nodes_m"][bar["start"]]), vec(frame["world_nodes_m"][bar["end"]])
            delta = [y-x for x, y in zip(a, b)]
            require(dot(delta, delta) == F(bar["length_m"])**2, "Display changes a member length")
            member_checks += 1
    return {"frames": len(frames), "node_controls": node_checks, "member_length_controls": member_checks,
            "role": "Display consistency only; continuum certificates are checked independently above"}


def check_capture(report, design, structural):
    model = structural["assembly"]
    require(report["loaded_capture_source_sha256"] == canonical(structural["capture_reference"]), "Loaded capture provenance mismatch")
    source = next(case for case in structural["cases"] if case["id"] == "independent_half_cantilevers/dead")
    limits = vec(design["motion"]["nominal_capture_half_widths_m"])
    for key in ("capture", "loaded_capture"):
        data = report[key]
        require(vec(data["capture_half_widths_m"]) == limits and len(data["connectors"]) == 4, "Capture box geometry changed")
        expected_within = True
        for pair, result in zip(model["connector_pairs"], data["connectors"]):
            left, right = pair["left"], pair["right"]
            require(result["left"] == left and result["right"] == right, "Wrong connector correspondence")
            # In the final installed frame the vertical rigid translation is
            # common and cancels. Read the actual verified structural solve,
            # independently of the convenience local capture export.
            if key == "capture":
                mismatch = [b-a for a, b in zip(vec(model["physical_nodes"][left]), vec(model["physical_nodes"][right]))]
            else:
                mismatch = [b-a for a, b in zip(vec(source["physical_deformed_coordinates_m"][left]), vec(source["physical_deformed_coordinates_m"][right]))]
            margins = [limit-abs(value) for limit, value in zip(limits, mismatch)]
            require(vec(result["nominal_difference_m"]) == mismatch and vec(result["margins_m"]) == margins, "Capture mismatch does not follow structural provenance")
            expected_within = expected_within and min(margins) >= 0
        require(data["within_envelope"] == expected_within, "Incorrect capture status")
    return {"nominal_connector_checks": 4, "loaded_connector_checks": 4,
            "loaded_source": source["id"], "scope": "Declared geometric capture box only; actual engagement, tolerance stack, contact and capacity remain unmodeled"}


def direct_installed_wrench(design, case, side, fixture):
    base = design["base_interfaces"]
    source = case["bases"][side]["reaction_on_bridge"]
    force = [-v for v in vec(source["force_kN"])]
    moment = [-v for v in vec(source["moment_kN_m"])]
    pivot = vec(base["pivot_position_in_base_m"])
    root = [pivot[0], pivot[1], F(base["truss_bottom_elevation_m"])]
    moment = add_vectors(moment, cross(root, force))
    gravity = F(design["gravity_m_s2"])/1000
    for item in base["fixed_mass_items_per_base"]:
        item_force = [F(0), F(0), -F(item["mass_kg"])*gravity]
        force = add_vectors(force, item_force)
        moment = add_vectors(moment, cross(vec(item["center_in_base_m"]), item_force))
    horizontal = [*fixture, F(0)]
    force = add_vectors(force, horizontal)
    moment = add_vectors(moment, cross([F(0), F(0), F(design["motion"]["horizontal_fixture_force_height_m"])], horizontal))
    return -force[2], [moment[1], -moment[0]]


def check_assembled_support(report, design, structural):
    data = report["assembled_support"]
    cases = {case["id"]: case for case in structural["cases"] if case["status"] == "unique_linear_solution"}
    expected = {(name, side, fx, fy) for name, side, (fx, fy) in product(cases, ("L", "R"), fixture_corners(design))}
    rows = data["case_support"]
    require(len(rows) == len(expected), "Installed support case inventory mismatch")
    lookup, pads = {}, pad_vertices(design)
    for row in rows:
        key = row["case_id"], row["side"], *vec(row["fixture_force_kN"])
        require(key in expected and key not in lookup, "Wrong or duplicate installed support case")
        lookup[key] = row
        case, side = cases[key[0]], key[1]
        source = case["bases"][side]
        require(row["source_reaction_on_bridge"] == source["reaction_on_bridge"] and row["source_reference_origin_m"] == source["reference_origin_m"], "Installed support substituted another root wrench")
        pivot = vec(design["base_interfaces"]["pivot_position_in_base_m"])
        require(vec(row["root_position_in_base_m"]) == [pivot[0], pivot[1], F(design["base_interfaces"]["truss_bottom_elevation_m"])], "Reported root lever arm differs from wrench transport")
        W, moments = direct_installed_wrench(design, case, side, key[2:])
        require(F(row["total_vertical_load_kN"]) == W and vec(row["required_first_moments_kN_m"]) == moments, "Full root wrench transport sign/origin error")
        require(row["resultant_position_in_base_m"] is None if W <= 0 else vec(row["resultant_position_in_base_m"]) == [q/W for q in moments], "Incorrect normal resultant position")
        least = None
        require(len(row["edges"]) == 4, "Installed support lost polygon edge")
        for edge, start, end in zip(row["edges"], pads, pads[1:]+pads[:1]):
            ex, ey = end[0]-start[0], end[1]-start[1]
            margin = ex*(moments[1]-start[1]*W)-ey*(moments[0]-start[0]*W)
            require(tuple(vec(edge["from"])) == start and tuple(vec(edge["to"])) == end and F(edge["signed_margin_kN_m2"]) == margin, "Installed support edge margin error")
            least = margin if least is None else min(least, margin)
        status = "outside" if W <= 0 or least < 0 else "boundary" if least == 0 else "inside"
        require(row["status"] == status, "Installed support status mismatch")
    release = data["release_support"]
    require(len(release) == 8, "Missing release side/force-corner certificate")
    seen = set()
    for record in release:
        side, fixture = record["side"], tuple(vec(record["fixture_force_kN"]))
        require((side, fixture) not in seen, "Repeated release certificate")
        seen.add((side, fixture))
        start, end, middle = [lookup[(name+"/dead", side, *fixture)] for name in ("joined_clamped", "joined_service", "joined_release_midpoint")]
        require(2*F(middle["total_vertical_load_kN"]) == F(start["total_vertical_load_kN"])+F(end["total_vertical_load_kN"]), "Release vertical load is not affine")
        require([2*v for v in vec(middle["required_first_moments_kN_m"])] == [a+b for a, b in zip(vec(start["required_first_moments_kN_m"]), vec(end["required_first_moments_kN_m"]))], "Release wrench is not affine")
        inside = min(F(start["total_vertical_load_kN"]), F(end["total_vertical_load_kN"])) > 0
        for actual, a, b in zip(record["edges"], start["edges"], end["edges"]):
            initial, final = F(a["signed_margin_kN_m2"]), F(b["signed_margin_kN_m2"])
            coefficients = [initial, final-initial]
            require(vec(actual["affine_margin_coefficients_kN_m2"]) == coefficients and F(actual["minimum_kN_m2"]) == minimum(coefficients)[0], "Complete release edge minimum is wrong")
            inside = inside and minimum(coefficients)[0] > 0
        require(record["actual_midpoint_wrench_matches"] is True and record["inside_for_complete_release"] == inside, "Release support status mismatch")
    require(data["service_support_inside"] == all(row["status"] == "inside" for row in rows if row["case_id"].startswith("joined_service/")), "Service support summary mismatch")
    require(data["release_support_inside"] == all(row["inside_for_complete_release"] for row in release), "Release summary mismatch")
    require(data["all_regular_case_support_inside"] == all(row["status"] == "inside" for row in rows), "Regular-state support summary mismatch")
    return {"independent_full_wrench_transports": len(rows), "affine_release_certificates": len(release),
            "whole_release_edge_minima": 4*len(release), "scope": "Vertical normal reaction existence; horizontal force, yaw, friction, anchorage and capacity are separate obligations"}


def check_clearance(report, design, structural):
    data = report["reference_clearance"]
    cases = {case["id"]: case for case in structural["cases"] if case["status"] == "unique_linear_solution"}
    require({row["case_id"] for row in data["cases"]} == set(cases) and len(data["cases"]) == len(cases), "Clearance case inventory mismatch")
    bottom = F(design["base_interfaces"]["truss_bottom_elevation_m"])
    radius = F(design["geometry"]["member_envelope_radius_m"])
    require(F(data["road_plane_m"]) == 0 and F(data["source_reference_m"]) == F("5.2832")
            and F(data["assumed_overlay_m"]) == F("0.1") and F(data["nominal_reserve_m"]) == F("0.1"), "Declared reference-plane constants changed")
    for row in data["cases"]:
        positions = cases[row["case_id"]]["physical_deformed_coordinates_m"]
        lows = [min(F(positions[bar["start"]][2]), F(positions[bar["end"]][2]))+bottom-radius
                for bar in structural["assembly"]["physical_members"]]
        lowest = min(lows)
        margin = lowest-F(data["source_reference_m"])-F(data["assumed_overlay_m"])
        require(F(row["lowest_declared_member_surface_m"]) == lowest, "Member-envelope clearance bound is wrong")
        require(F(row["margin_above_reference_plus_overlay_m"]) == margin == F(row["remaining_nominal_reserve_m"]), "Remaining reserve arithmetic mismatch")
        require(row["satisfies_reference_plus_overlay"] == (margin >= 0), "Reference-clearance status mismatch")
    require(data["service_reference_clearance_satisfied"] == all(row["satisfies_reference_plus_overlay"] for row in data["cases"] if row["case_id"].startswith("joined_service/")), "Reference-clearance summary mismatch")
    return {"cases": len(cases), "member_endpoint_envelope_checks": len(cases)*len(structural["assembly"]["physical_members"]),
            "scope": "Arithmetic against declared flat datum/reference only; this audit does not determine regulatory applicability or actual surveyed clearance"}


def check_report(report, design, structural, check_structure=True):
    model = structural["assembly"]
    require(report["schema_version"] == 1 and report["format"] == "crest-root-shared-motion/3" and report["physical_approval"] is False, "Incorrect motion scope/schema")
    require(report["design_sha256"] == assembly.design_digest(design), "Motion design digest mismatch")
    require(report["geometry_and_mass_source_sha256"] == canonical(model) and report["structural_report_sha256"] == canonical(structural), "Motion source bindings mismatch")
    if check_structure:
        verify_assembly.check_report(structural, design)
    certificates = {
        "rigid_paths": check_paths(report, design, model),
        "deployment_support": check_support(report, design, model),
        "separation": check_separation(report, design, model),
        "sweep_enclosures": check_sweep(report, design, model),
        "display_frames": check_frames(report, design, model),
        "capture": check_capture(report, design, structural),
        "assembled_support": check_assembled_support(report, design, structural),
        "reference_clearance": check_clearance(report, design, structural),
    }
    expected_checks = {
        "continuous_support_inside": all(row["status"] == "inside" for row in report["support_certificates"]),
        "continuous_half_separation": all(row["status"] != "overlap" and row["outside_joint_zone_separated"] for row in report["separation_certificates"]),
        "nominal_capture_within_envelope": report["capture"]["within_envelope"],
        "loaded_capture_within_envelope": report["loaded_capture"]["within_envelope"],
        "stage_endpoints_match": True,
        "service_support_inside": report["assembled_support"]["service_support_inside"],
        "release_support_inside": report["assembled_support"]["release_support_inside"],
        "service_reference_clearance_satisfied": report["reference_clearance"]["service_reference_clearance_satisfied"],
    }
    require(report["checks"] == expected_checks, "Top-level motion checks disagree with their certificates")
    return certificates


def negative_controls(report, design, structural):
    rejected, failed_scenarios = [], []
    model = structural["assembly"]

    def rejects(name, function):
        try:
            function()
        except (ValueError, ArithmeticError, KeyError):
            rejected.append(name)
        else:
            raise ValueError("Negative control accepted: "+name)

    bad = deepcopy(design)
    bad["motion"]["retreat_m"] = "0"
    rows = motion.separation_certificates(bad, assembly.build_assembly(bad))
    require(any(row["status"] == "overlap" for row in rows), "Zero retreat should fail centerline separation")
    failed_scenarios.append("zero_retreat_overlap")
    bad = deepcopy(design)
    for item in bad["base_interfaces"]["fixed_mass_items_per_base"]:
        item["mass_kg"] = "0"
    rows = motion.support_certificates(bad, assembly.build_assembly(bad))
    require(any(row["status"] == "outside" for row in rows), "Absent fixed base mass should fail this support screen")
    failed_scenarios.append("absent_fixed_base_mass_support_failure")
    bad = deepcopy(design)
    bad["motion"]["horizontal_fixture_force_bounds_kN"] = {"x": ["-1000", "1000"], "y": ["-600", "600"]}
    rows = motion.support_certificates(bad, assembly.build_assembly(bad))
    require(any(row["status"] == "outside" for row in rows), "Large declared fixture force should fail support")
    failed_scenarios.append("excess_fixture_force_support_failure")
    bad = deepcopy(design)
    bad["motion"]["horizontal_fixture_force_bounds_kN"] = {"x": ["0", "0"], "y": ["0", "0"]}
    weight, moments = direct_deployment_wrench(bad, assembly.build_assembly(bad), "L", "lift", F(0), (F(0), F(0)))
    bad["base_interfaces"]["pad_rectangle_m"]["x_min"] = str(moments[0]/weight)
    rows = motion.support_certificates(bad, assembly.build_assembly(bad))
    require(all(row["status"] == "boundary" for row in rows if row["side"] == "L" and row["stage_id"] == "lift"), "Exact no-lift equality must be marginal, not strict inside")
    failed_scenarios.append("no_lift_equality_is_boundary")
    bad = deepcopy(design)
    bad["motion"]["nominal_capture_half_widths_m"] = ["0.000001"]*3
    loaded = motion.nominal_capture(bad, model, structural["capture_reference"]["deformed_local_nodes_by_half"])
    require(loaded["within_envelope"] is False, "Tightened capture envelope should reject elastic mismatch")
    failed_scenarios.append("tightened_loaded_capture_failure")
    damaged = deepcopy(report)
    damaged["rational_paths"]["L"]["slew"]["rotation_numerator"][0][0][0] = "1"
    rejects("tampered_rotation", lambda: check_paths(damaged, design, model))
    damaged = deepcopy(report)
    damaged["support_certificates"][0]["edges"][0]["minimum"] = "1"
    rejects("tampered_path_minimum", lambda: check_support(damaged, design, model))
    damaged = deepcopy(report)
    damaged["separation_certificates"][1]["x_gap_numerator_m"][0] = "1"
    rejects("tampered_prism_gap", lambda: check_separation(damaged, design, model))
    damaged = deepcopy(report)
    damaged["swept_boxes_model"][0]["maximum_m"][2] = damaged["swept_boxes_model"][0]["minimum_m"][2]
    rejects("clipped_swept_envelope", lambda: check_sweep(damaged, design, model))
    damaged = deepcopy(report)
    damaged["loaded_capture"]["connectors"][0]["nominal_difference_m"] = ["0"]*3
    rejects("substituted_undeformed_capture", lambda: check_capture(damaged, design, structural))
    damaged = deepcopy(report)
    damaged["assembled_support"]["case_support"][0]["required_first_moments_kN_m"][0] = "0"
    rejects("tampered_root_wrench_transport", lambda: check_assembled_support(damaged, design, structural))
    damaged = deepcopy(report)
    damaged["reference_clearance"]["cases"][0]["lowest_declared_member_surface_m"] = "10"
    rejects("tampered_clearance", lambda: check_clearance(damaged, design, structural))
    damaged = deepcopy(report)
    damaged["structural_report_sha256"] = "0"*64
    rejects("mismatched_structural_provenance", lambda: check_report(damaged, design, structural, False))
    return {"scenario_controls": failed_scenarios, "tampering_rejections": rejected,
            "total_controls": len(failed_scenarios)+len(rejected)}


def verify_all(report_path=HERE/"motion-report.json", structural_path=HERE/"assembly-report.json", design_path=HERE/"design.json", replay=False):
    report_path, structural_path, design_path = map(Path, (report_path, structural_path, design_path))
    report, structural = json.loads(report_path.read_text()), json.loads(structural_path.read_text())
    design = assembly.load_design(design_path)
    certificates = check_report(report, design, structural)
    controls = negative_controls(report, design, structural)
    if replay:
        require(motion.run_motion(design, structural["assembly"], structural) == report, "Deterministic motion replay differs from saved report")
    paths = {"motion.py": HERE/"motion.py", "assembly.py": HERE/"assembly.py", "verify_assembly.py": HERE/"verify_assembly.py",
             "verify_motion.py": Path(__file__), "v2/engineering.py": HERE.parent/"v2"/"engineering.py"}
    return {"schema_version": 1, "format": "crest-root-independent-motion-verification/3", "checks_passed": True, "physical_approval": False,
            "design_sha256": assembly.design_digest(design), "geometry_and_mass_source_sha256": report["geometry_and_mass_source_sha256"],
            "motion_report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
            "assembly_report_sha256": hashlib.sha256(structural_path.read_bytes()).hexdigest(),
            "source_sha256": {name: hashlib.sha256(path.read_bytes()).hexdigest() for name, path in paths.items()},
            "certificates": certificates, "controls": controls, "complete_production_replay": replay,
            "review_findings": ["No transform, normal-support moment-sign, interval-minimum, separation or capture-provenance defect found for the declared reference inputs",
                                "Support certificates establish existence of vertical coplanar reactions only, not sliding/yaw resistance, individual pad loads or capacity",
                                "Sweeps concern declared nominal member/floor/mass-point geometry; actual guards, locks, fixture bodies and reserved-zone engagement remain outside the model",
                                "Reference-clearance arithmetic uses a hypothetical flat datum and declared small-displacement axial shapes; legal applicability and as-built clearance are not established"]}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=HERE/"motion-report.json")
    parser.add_argument("--structural-report", type=Path, default=HERE/"assembly-report.json")
    parser.add_argument("--design", type=Path, default=HERE/"design.json")
    parser.add_argument("--replay", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    result = verify_all(args.report, args.structural_report, args.design, args.replay)
    rendered = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
