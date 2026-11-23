"""Independent verification of the inspected LW-159-R1 hand-turn reference."""

import argparse
from fractions import Fraction as F
import itertools
import json
import math
from pathlib import Path

from bounds import Box, S, TWO, circle_arcs, plus_pose, segment_gap

HERE = Path(__file__).resolve().parent


def exact_json(path):
    return json.loads(path.read_text(encoding="utf-8"), parse_float=F, parse_int=F)


def rounded(value):
    return round(float(value), 6)


def down(value):
    """Display at most the proved lower bound in mm."""
    return math.floor(value * 1000) / 1000


def min_angle_bounds(lengths):
    g, a, b, c = lengths
    dmin = max(F(0), g[0] - a[1], a[0] - g[1])
    dmax = g[1] + a[1]
    fg = dmin - max(abs(b[0] - c[1]), abs(b[1] - c[0]))
    rt = b[0] + c[0] - dmax
    if min(fg, rt) <= 0:
        raise ValueError("strict closure fails in the tolerance box")
    # With d at its minimum, maximum cos(mu) occurs at a corner: for fixed c,
    # the only interior stationary value in b is a minimum, and conversely.
    hi = max((x * x + y * y - dmin * dmin) / (2 * x * y) for x in b for y in c)
    numerator = b[0] ** 2 + c[0] ** 2 - dmax**2
    lo = numerator / (2 * (b[0] * c[0] if numerator < 0 else b[1] * c[1]))
    abs_cos = max(abs(hi), abs(lo))
    if abs_cos >= 1:
        raise ValueError("transmission lower bound unavailable")
    return {
        "inner_closure_mm": str(fg),
        "outer_closure_mm": str(rt),
        "max_abs_transmission_cos_upper": str(abs_cos),
        "acute_transmission_lower_deg": down(math.degrees(math.acos(abs_cos))),
    }


def arc_cos(numerator, denominator):
    v = numerator / denominator
    if v.lo <= -S or v.hi >= S:
        raise ValueError("rocker reversal triangle nearly flat; check dimensions")
    return tuple(math.degrees(math.acos(float(y) / S)) for y in (v.hi, v.lo))


def reversals(g, a, b, c):
    # The reference has g>a and b>a. At the rocker extrema O,A,B are collinear.
    if b.lo <= a.hi or g.lo <= a.hi:
        raise ValueError("expected reference crank-rocker arrangement")
    angles = []
    for r in (b + a, b - a):
        phi = arc_cos(r.square() - g.square() - c.square(), TWO * g * c)
        direction = arc_cos(r.square() + g.square() - c.square(), TWO * r * g)
        angles.append((phi, direction))
    return angles


def stroke_report(nominal, tol, subdivisions=4):
    def bands(dim):
        step = 2 * tol / subdivisions
        return [
            Box.enclose(dim - tol + i * step, dim - tol + (i + 1) * step)
            for i in range(subdivisions)
        ]

    swings, ratios = [], []
    for dims in itertools.product(*(bands(v) for v in nominal)):
        (outer, t1), (inner, t2) = reversals(*dims)
        swings.extend([inner[0] - outer[1], inner[1] - outer[0]])
        delta = (180 + t2[0] - t1[1], 180 + t2[1] - t1[0])
        if not 180 < delta[0] <= delta[1] < 360:
            raise ValueError("long/short stroke assumption fails")
        ratios.extend(s / (360 - s) for s in delta)
    g, a, b, c = nominal
    r_outer, r_inner = b + a, b - a

    def in_degrees(cosine):
        return math.degrees(math.acos(cosine))

    p1, p2 = [
        in_degrees((r * r - g * g - c * c) / (2 * g * c)) for r in (r_outer, r_inner)
    ]
    r1, r2 = [
        in_degrees((r * r + g * g - c * c) / (2 * g * r)) for r in (r_outer, r_inner)
    ]
    r2 += 180
    delta = r2 - r1
    return {
        "nominal_rocker_extrema_deg": [rounded(p1), rounded(p2)],
        "nominal_crank_reversal_deg_ccw": [rounded(r1), rounded(r2)],
        "nominal_swing_deg": rounded(p2 - p1),
        "nominal_stroke_ratio": rounded(delta / (360 - delta)),
        "fixed_spacing_box_swing_deg": [
            down(min(swings) - 1e-7),
            -down(-max(swings) - 1e-7),
        ],
        "fixed_spacing_box_stroke_ratio": [
            down(min(ratios) - 1e-7),
            -down(-max(ratios) - 1e-7),
        ],
    }


def axial(spec):
    """Parallel plates/shafts; stack corner bounds include play and measured flatness."""
    board, links, j = spec["base"], spec["links"], spec["joint"]

    def ends(obj, key, tol):
        return obj[key] - obj[tol], obj[key] + obj[tol]

    h0, h1 = ends(j, "collar_height", "collar_height_tol")
    w0, w1 = ends(j, "washer_thickness", "washer_thickness_tol")
    t0, t1 = ends(links, "thickness", "thickness_tol")
    s = j["spacer_length_tol"]
    play = j["endplay"][1]
    stick = j["free_end_projection"][1]
    bow = links["flatness_max"]
    flat = board["flatness_max"]
    lower_min = h0 + j["ground_spacer"] - s + w0 - flat - bow
    lower_max = h1 + j["ground_spacer"] + s + w1 + flat + bow + play
    ground_top = lower_max + t1 + w1 + play + h1 + stick
    upper_floor = lower_min + t0 + 2 * w0 + j["moving_spacer"] - s - bow
    lower_back = lower_min - w1 - play - h1 - stick
    result = {
        "moving_shafts_over_base_mm": lower_back - flat,
        "crank_rocker_over_base_mm": lower_min - flat,
        "coupler_over_ground_pin_stacks_mm": upper_floor - ground_top,
        "coupler_over_lower_links_mm": upper_floor - (lower_max + t1),
        "ground_retention_over_table_mm": board["foot_height"]
        - board["foot_height_tol"]
        - flat
        - h1
        - stick,
    }
    return result


def sweep(spec, nominal, radial, subdivisions=128):
    board, links, j = spec["base"], spec["links"], spec["joint"]
    fixed, free = radial
    tol = spec["spacings"]["acceptance_plus_minus"]
    dims = [
        (v - tol - 2 * r, v + tol + 2 * r)
        for v, r in zip(nominal, (fixed, free, free, free))
    ]
    lengths = [Box.enclose(*ends) for ends in dims]
    closure = min_angle_bounds(dims)
    # An actual drilled line and its two effective shaft centres can differ in
    # direction. Bound both angles using the perpendicular offset <= 2*free.
    eps = sum(
        math.degrees(math.asin(2 * free / min(ends))) for ends in (dims[2], dims[3])
    )
    closure["cut_link_acute_transmission_lower_deg"] = down(
        closure["acute_transmission_lower_deg"] - eps - 1e-7
    )
    axial_gaps = axial(spec)
    # Each link centre is offset from the line of shaft axes by at most free,
    # growing past B on the integral pointer by normalized-axis perturbation.
    # The extension uses actual hole-axis direction, versus shaft-axis direction:
    # |a/|a|-b/|b|| <= 2|a-b|/min(|a|,|b|).
    ext = links["pointer_extension"] + links["pointer_extension_tol"]
    half = (links["width"] + links["width_tol"]) / 2 + links[
        "edge_profile_deviation_max"
    ]
    body_radii = {
        "crank": half + free,
        "coupler": half + free,
        "rocker": half + free + ext * 4 * free / dims[3][0],
    }
    # A loose washer/spacer can shift by its bore-to-shaft radial clearance.
    radial_hardware = max(
        j["collar_od_envelope_max"] / 2,
        (j["washer_od_max"] + j["washer_bore"][1] - j["shaft_diameter"][0]) / 2,
        (j["spacer_od_max"] + j["spacer_bore"][1] - j["shaft_diameter"][0]) / 2,
    )
    grip = (j["grip_diameter"] + j["grip_diameter_tol"]) / 2
    grip += (j["grip_blind_bore_diameter"][1] - j["shaft_diameter"][0]) / 2
    radii = {**body_radii, **{key: radial_hardware for key in ("O", "D", "A", "B")}}
    radii["A"] = max(radial_hardware, grip)
    parts = {
        "crank": ("O", "A"),
        "coupler": ("A", "B"),
        "rocker": ("D", "T"),
        **{key: (key, key) for key in "ODAB"},
    }
    # Each allowed coaxial fastener joint or washer/face contact is intentional.
    # Upper/lower link contact remains checked through their z separation.
    joint_sets = {"crank": "OA", "coupler": "AB", "rocker": "DB"}
    pairs = [("crank", "rocker")]
    pairs += [
        (body, joint)
        for body, ends in joint_sets.items()
        for joint in "ODAB"
        if joint not in ends
    ]
    pairs += list(itertools.combinations("ODAB", 2))
    gap = {p: None for p in pairs}
    boardgap = {side: None for side in ("left", "right", "bottom", "top")}
    angle_sine = 2 * (board["fixed_position_tol_per_axis"] + fixed) / dims[0][0]
    if angle_sine >= F(1, 10):
        raise ValueError("ground datum orientation too uncertain")
    # 1-cos(beta) = sin(beta)^2/(1+cos(beta)) <= sin(beta)^2, |beta|<90deg.
    turn_err = Box.enclose(angle_sine)
    cos_err = Box.enclose(angle_sine**2)
    offset = Box.enclose(board["fixed_position_tol_per_axis"] + fixed)
    ext_ratio = (
        Box.enclose(links["pointer_extension"] - links["pointer_extension_tol"], ext)
        / lengths[3]
    )
    for cs in circle_arcs(subdivisions):
        pos = plus_pose(*lengths, cs)
        pos["T"] = tuple(b + (b - d) * ext_ratio for b, d in zip(pos["B"], pos["D"]))
        # Intersect loose interval expressions with the known rocker disk.
        # Each coordinate differs from the corresponding D coordinate by at
        # most the pin-axis length, plus the physical tip extension for T.
        for key, cap in (
            ("B", lengths[3].hi),
            ("T", lengths[3].hi + Box.enclose(ext).hi),
        ):
            pos[key] = tuple(
                Box(max(coord.lo, d.lo - cap), min(coord.hi, d.hi + cap))
                for coord, d in zip(pos[key], pos["D"])
            )
        for p in pairs:
            lhs, rhs = p
            raw = segment_gap(*(pos[key] for key in (*parts[lhs], *parts[rhs])))
            margin = F(raw, S) - radii[lhs] - radii[rhs]
            if gap[p] is None or margin < gap[p]:
                gap[p] = margin
        for name, keys in parts.items():
            radius = radii[name]
            for key in set(keys):
                x, y = pos[key]
                ax = Box(0, max(abs(x.lo), abs(x.hi)))
                ay = Box(0, max(abs(y.lo), abs(y.hi)))
                ex = (ax * cos_err + ay * turn_err + offset).hi
                ey = (ay * cos_err + ax * turn_err + offset).hi
                edges = [
                    F(x.lo - ex, S) - board["xmin"],
                    board["xmax"] - F(x.hi + ex, S),
                    F(y.lo - ey, S) - board["ymin"],
                    board["ymax"] - F(y.hi + ey, S),
                ]
                for side, edge in zip(boardgap, edges):
                    m = edge - radius - board["outline_tol"]
                    if boardgap[side] is None or m < boardgap[side]:
                        boardgap[side] = m
    # Coupler may pass over a ground joint in plan; all of that joint lies below it.
    # For these two pairs use the stronger valid separating dimension.
    for jname in "OD":
        gap[("coupler", jname)] = max(
            gap[("coupler", jname)], axial_gaps["coupler_over_ground_pin_stacks_mm"]
        )
    # Both ground collars have y=0 nominally. Either bottom or top foot row
    # can only shift toward them by the stated location and board-edge errors.
    foot_y = min(
        abs(board["ymin"] + board["foot_inset"]),
        abs(board["ymax"] - board["foot_inset"]),
    )
    foot_y -= board["outline_tol"] + board["foot_position_tol"]
    # 0.708 > sqrt(2)/2, enclosing an arbitrarily rotated square foot.
    foot_radius = F(708, 1000) * board["foot_width"]
    foot_to_collar = (
        foot_y
        - foot_radius
        - board["fixed_position_tol_per_axis"]
        - fixed
        - radial_hardware
    )
    foot_edge = board["foot_inset"] - board["foot_position_tol"] - foot_radius
    if foot_edge <= 0:
        raise ValueError("corner foot may overhang the deck")
    floor = spec["acceptance"]["minimum_measured_unintended_gap"]
    for label, bound in {
        **gap,
        **axial_gaps,
        **boardgap,
        "foot to lower collar": foot_to_collar,
    }.items():
        if bound < floor:
            raise ValueError(
                f"{label}: clearance lower bound {float(bound):.3f} mm falls below {float(floor)}"
            )
    return {
        "crank_arcs": 4 * subdivisions,
        "effective_pin_axis_lengths_mm": [list(map(str, d)) for d in dims],
        "closure_and_transmission": closure,
        "nonincident_pair_clearance_lower_mm": {
            " / ".join(k): down(v) for k, v in gap.items()
        },
        "parallel_stack_axial_clearance_lower_mm": {
            k: down(v) for k, v in axial_gaps.items()
        },
        "inside_board_edge_lower_mm": {k: down(v) for k, v in boardgap.items()},
        "ground_bottom_collar_to_foot_lower_mm": down(foot_to_collar),
        "foot_inside_its_measured_edges_mm": down(foot_edge),
    }


def registration(nominal, extension):
    dimensions = [Box.enclose(value) for value in nominal]
    rows = []
    for turn, axis in enumerate(((1, 0), (0, 1), (-1, 0), (0, -1))):
        positions = plus_pose(*dimensions, tuple(Box.enclose(x) for x in axis))
        ratio = Box.enclose(extension) / dimensions[3]
        positions["T"] = tuple(
            b + (b - d) * ratio for b, d in zip(positions["B"], positions["D"])
        )
        row = {"crank_deg": 90 * turn}
        for name in ("A", "B", "T"):
            row[name + "_mm"] = [
                round((v.lo + v.hi) / (2 * S), 3) for v in positions[name]
            ]
        rows.append(row)
    return rows


def verify(spec):
    if spec["units"] != "mm":
        raise ValueError("expected millimetres")
    j = spec["joint"]
    free = (j["moving_bore"][1] - j["shaft_diameter"][0]) / 2
    fixed = (j["ground_bore"][1] - j["shaft_diameter"][0]) / 2
    if free < 0 or j["moving_bore"][0] <= j["shaft_diameter"][1]:
        raise ValueError("moving bores require free clearance")
    if j["ground_bore"][0] < j["shaft_diameter"][1] or fixed < 0:
        raise ValueError("ground fit incompatible with shafts")
    for part in ("washer_bore", "spacer_bore", "grip_blind_bore_diameter"):
        if j[part][0] <= j["shaft_diameter"][1]:
            raise ValueError(f"{part}: part may not slide onto the shaft")
    if spec["datum"]["branch"] != "cross(D-A, B-A) > 0":
        raise ValueError("this reference design requires the upper assembly mode")
    washer_radial = (j["washer_od_min"] - j["moving_bore"][1]) / 2
    # The shaft and the two independently loose parts can move to opposite sides.
    washer_radial -= free + (j["washer_bore"][1] - j["shaft_diameter"][0]) / 2
    # Smallest support annulus when loose washers/spacers sit off-centre.
    washer_ecc = (j["washer_bore"][1] - j["shaft_diameter"][0]) / 2
    spacer_ecc = (j["spacer_bore"][1] - j["shaft_diameter"][0]) / 2
    spacer_support = (
        (j["spacer_od_min"] - j["washer_bore"][1]) / 2 - spacer_ecc - washer_ecc
    )
    collar_support = (j["collar_od_min"] - j["washer_bore"][1]) / 2 - washer_ecc
    if min(spacer_support, collar_support) <= 0:
        raise ValueError("spacer or collar may escape the washer")
    if washer_radial <= 0:
        raise ValueError("moving bore may escape washer")
    wall = (
        spec["links"]["width"] - spec["links"]["width_tol"] - j["moving_bore"][1]
    ) / 2
    wall -= spec["links"]["edge_profile_deviation_max"]
    if wall <= 0:
        raise ValueError("moving bore may break through its end profile")
    nominal = [spec["spacings"][k] for k in ("ground", "crank", "coupler", "rocker")]
    if (
        nominal[0] != spec["datum"]["D"][0]
        or spec["datum"]["D"][1] != 0
        or list(spec["datum"]["O"]) != [0, 0]
    ):
        raise ValueError("base datum disagrees with ground spacing")
    return {
        "revision": spec["revision"],
        "reference_units": "mm",
        "model": "rigid, parallel shafts and link planes; upper assembly branch; installed lean and flex tested on bench",
        "minimum_washer_overlap_at_link_mm": str(washer_radial),
        "minimum_spacer_overlap_at_washer_mm": str(spacer_support),
        "minimum_collar_overlap_at_washer_mm": str(collar_support),
        "minimum_radial_material_at_end_bore_mm": str(wall),
        "max_radial_offsets_mm": {"fixed": str(fixed), "moving": str(free)},
        "motion": stroke_report(nominal, spec["spacings"]["acceptance_plus_minus"]),
        "nominal_bench_coordinates": registration(
            nominal, spec["links"]["pointer_extension"]
        ),
        "geometry": sweep(spec, nominal, (fixed, free)),
    }


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--spec", type=Path, default=HERE.parent / "kit" / "design.json")
    ap.add_argument(
        "--check", action="store_true", help="recompute and compare released receipt"
    )
    args = ap.parse_args()
    result = verify(exact_json(args.spec))
    if args.check:
        expected = json.loads((HERE / "receipt.json").read_text())
        if result != expected:
            raise SystemExit("FAIL: recomputed engineering receipt differs")
        print(
            "PASS: reference dimensions, motion, tolerance box, collision pairs, axial stacks and footprint"
        )
    else:
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
