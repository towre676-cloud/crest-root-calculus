"""One exact, explicitly hypothetical 3D deck-truss assembly for v3 research.

Public API: load_design(path=DEFAULT_DESIGN), build_assembly(design),
run_assembly(design). All returned scalars are JSON-compatible exact strings;
counts/indices are integers. The unchanged v2 axial-truss kernel supplies each
linear solve. Physical nodes, bars and true mass items persist across modes;
perfect seam ties are represented by node aliases, never by removing members.

Coordinates: x along the installed span, y across, z up; local roots are at
x=0,y=0,z=0. The right half is placed by a proper 180-degree z rotation and
translation (23-x,-y,z). Forces are kN; moments r cross F are kN m. Stress is
kN/m^2. Floor/rail masses are above the top chord. Their vertical forces are
transferred to top nodes preserving x/y first moments, NOT their true z CG.
The true mass ledger, rather than the force-node CG, is used for rigid motion.
"""

import argparse
from copy import deepcopy
from fractions import Fraction
import hashlib
import importlib.util
import json
from math import isqrt
from pathlib import Path

HERE = Path(__file__).resolve().parent
DEFAULT_DESIGN = HERE / "design.json"
TRUSS_PATH = HERE.parent / "v2" / "truss.py"
_spec = importlib.util.spec_from_file_location("crest_v2_truss_for_v3", TRUSS_PATH)
_truss = importlib.util.module_from_spec(_spec)
_spec.loader.exec_module(_truss)

ROLES = ("longitudinal", "transverse", "vertical", "xy_diagonal", "xz_diagonal", "yz_diagonal")
LIMITATIONS = [
    "Reference geometry, areas, material, masses and crowd pressures are declared hypotheses, not observations or selected capacities",
    "Linear pin-jointed axial bars only: no member bending, torsion, buckling, geometric nonlinearity, yielding or fatigue checks",
    "No deck plate, guards, connections, center locks, mast, bearings, actuators, chassis, pads, soil or anchorage capacity model",
    "No dynamic, vibration, impact, wind-on-truss or moving-load response model",
    "The rigid/deformable mode change uses ideal ties and prescribed translations, not a validated physical control sequence",
    "No accessibility, pedestrian layout, transport envelope, fabrication or operational approval",
]


class DesignError(ValueError):
    pass


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (str, int, Fraction)):
        raise DesignError("Use integers or exact rational/decimal strings, not floats")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as exc:
        raise DesignError("Invalid exact scalar") from exc


def _keys(value, keys, label):
    if not isinstance(value, dict) or set(value) != set(keys):
        raise DesignError(label + " must contain exactly " + ", ".join(keys))


def _positive(value, label, allow_zero=False):
    value = exact(value)
    if value < 0 or (value == 0 and not allow_zero):
        raise DesignError(label + " must be " + ("nonnegative" if allow_zero else "positive"))
    return value


def _v(values):
    return [exact(value) for value in values]


def _s(values):
    return [str(value) for value in values]


def _cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def _sum_vectors(vectors):
    answer = [Fraction(0)]*3
    for vector in vectors:
        answer = [a+b for a, b in zip(answer, vector)]
    return answer


def _rational_sqrt(value):
    a, b = isqrt(value.numerator), isqrt(value.denominator)
    if a*a != value.numerator or b*b != value.denominator or not a:
        raise DesignError("Every bar must have a positive rational length")
    return Fraction(a, b)


def design_digest(design):
    """Digest of sorted compact UTF-8 JSON, independent of file whitespace."""
    return hashlib.sha256(json.dumps(design, sort_keys=True, separators=(",", ":"), ensure_ascii=False).encode()).hexdigest()


def load_design(path=DEFAULT_DESIGN):
    design = json.loads(Path(path).read_text(encoding="utf-8"))
    _validate(design)
    return design


def _validate(design):
    _keys(design, ("schema", "id", "basis", "units", "geometry", "material", "member_areas_m2", "mass", "gravity_m_s2", "crowd", "motion", "base_interfaces"), "design")
    if design["schema"] != "crest-root-deployable-crossing-design-v3":
        raise DesignError("Unsupported design schema")
    for field in ("id", "basis"):
        if not isinstance(design[field], str) or not design[field].strip():
            raise DesignError(field + " must be a nonempty string")
    if design["units"] != {"length": "m", "force": "kN", "mass": "kg"}:
        raise DesignError("This assembly uses m, kN and kg explicitly")
    geometry = design["geometry"]
    _keys(geometry, ("span_m", "panels_per_half", "brick_ratio", "floor_width_m", "floor_thickness_m", "member_envelope_radius_m", "joint_contact_zone_length_m"), "geometry")
    n = geometry["panels_per_half"]
    if type(n) is not int or n < 1:
        raise DesignError("panels_per_half must be a positive integer")
    _keys(geometry["brick_ratio"], ("panel", "width", "depth"), "brick_ratio")
    ratios = [_positive(geometry["brick_ratio"][key], key) for key in ("panel", "width", "depth")]
    span = _positive(geometry["span_m"], "span")
    width = span*ratios[1]/(2*n*ratios[0])
    if _positive(geometry["floor_width_m"], "floor width") > width:
        raise DesignError("Floor width must lie within the reference top-chord width")
    _positive(geometry["floor_thickness_m"], "floor thickness")
    _positive(geometry["member_envelope_radius_m"], "declared member envelope radius")
    if _positive(geometry["joint_contact_zone_length_m"], "reserved joint contact zone") >= span/2:
        raise DesignError("Reserved joint contact zone must be shorter than a half-span")
    _keys(design["material"], ("young_modulus_kN_m2", "density_kg_m3"), "material")
    for key, value in design["material"].items():
        _positive(value, key)
    _keys(design["member_areas_m2"], ROLES, "member areas")
    for key, value in design["member_areas_m2"].items():
        _positive(value, key)
    _keys(design["mass"], ("floor_kg_m2", "rail_kg_m", "rail_center_above_walking_surface_m", "joint_kg_per_physical_node"), "mass")
    for key, value in design["mass"].items():
        _positive(value, key, allow_zero=True)
    _positive(design["gravity_m_s2"], "gravity")
    _keys(design["crowd"], ("symmetric_pressure_kN_m2", "asymmetric_pressure_kN_m2", "asymmetric_region"), "crowd")
    for key in ("symmetric_pressure_kN_m2", "asymmetric_pressure_kN_m2"):
        _positive(design["crowd"][key], key, allow_zero=True)
    if design["crowd"]["asymmetric_region"] != "global_y_positive_half_of_floor":
        raise DesignError("Unsupported asymmetric crowd region")
    motion = design["motion"]
    _keys(motion, ("retreat_m", "yaw_interval_count", "nominal_capture_half_widths_m", "horizontal_fixture_force_bounds_kN", "horizontal_fixture_force_height_m"), "motion")
    _positive(motion["retreat_m"], "retreat", allow_zero=True)
    _positive(motion["horizontal_fixture_force_height_m"], "fixture force height", allow_zero=True)
    if type(motion["yaw_interval_count"]) is not int or motion["yaw_interval_count"] <= 0:
        raise DesignError("Yaw interval count must be a positive integer")
    if not isinstance(motion["nominal_capture_half_widths_m"], list) or len(motion["nominal_capture_half_widths_m"]) != 3:
        raise DesignError("Capture box needs three half-widths")
    for value in motion["nominal_capture_half_widths_m"]:
        _positive(value, "capture half-width")
    _keys(motion["horizontal_fixture_force_bounds_kN"], ("x", "y"), "horizontal force bounds")
    for bounds in motion["horizontal_fixture_force_bounds_kN"].values():
        if not isinstance(bounds, list) or len(bounds) != 2 or exact(bounds[0]) > exact(bounds[1]):
            raise DesignError("Force bounds must be an ordered exact pair")
    base = design["base_interfaces"]
    _keys(base, ("truss_bottom_elevation_m", "transport_truss_bottom_m", "pivot_position_in_base_m", "pad_rectangle_m", "fixed_mass_items_per_base", "interpretation"), "base interfaces")
    for key in ("truss_bottom_elevation_m", "transport_truss_bottom_m"):
        _positive(base[key], key, allow_zero=True)
    if not isinstance(base["pivot_position_in_base_m"], list) or len(base["pivot_position_in_base_m"]) != 3:
        raise DesignError("Base pivot requires three exact coordinates")
    _v(base["pivot_position_in_base_m"])
    _keys(base["pad_rectangle_m"], ("x_min", "x_max", "y_min", "y_max"), "pad rectangle")
    pads = {key: exact(value) for key, value in base["pad_rectangle_m"].items()}
    if pads["x_min"] >= pads["x_max"] or pads["y_min"] >= pads["y_max"]:
        raise DesignError("Pad rectangle must have positive area")
    if not isinstance(base["fixed_mass_items_per_base"], list) or not base["fixed_mass_items_per_base"]:
        raise DesignError("Declare the fixed base mass items")
    seen = set()
    for item in base["fixed_mass_items_per_base"]:
        _keys(item, ("id", "mass_kg", "center_in_base_m"), "fixed mass item")
        if not isinstance(item["id"], str) or not item["id"] or item["id"] in seen:
            raise DesignError("Fixed mass item IDs must be nonempty and unique")
        seen.add(item["id"])
        _positive(item["mass_kg"], "fixed mass", allow_zero=True)
        if not isinstance(item["center_in_base_m"], list) or len(item["center_in_base_m"]) != 3:
            raise DesignError("Fixed mass requires three coordinates")
        _v(item["center_in_base_m"])


def _placed(local, side, span):
    x, y, z = local
    return [x, y, z] if side == "L" else [span-x, -y, z]


def build_assembly(design):
    """Return one physical model/true mass ledger, not a capacity assessment."""
    _validate(design)
    span = exact(design["geometry"]["span_m"])
    n = design["geometry"]["panels_per_half"]
    panel = span/(2*n)
    ratio = design["geometry"]["brick_ratio"]
    width = panel*exact(ratio["width"])/exact(ratio["panel"])
    height = panel*exact(ratio["depth"])/exact(ratio["panel"])
    floor_width = exact(design["geometry"]["floor_width_m"])
    thickness = exact(design["geometry"]["floor_thickness_m"])
    density = exact(design["material"]["density_kg_m3"])
    gravity = exact(design["gravity_m_s2"])
    halves, physical_nodes, physical_members, mass_items = {}, {}, [], []
    dead = {}
    for side in ("L", "R"):
        local = {f"{side}_{i}_{j}_{k}": [i*panel, (j-Fraction(1, 2))*width, k*height]
                 for i in range(n+1) for j in (0, 1) for k in (0, 1)}
        world = {name: _placed(position, side, span) for name, position in local.items()}
        physical_nodes.update({name: _s(position) for name, position in world.items()})
        dead.update({name: [Fraction(0)]*3 for name in local})
        members, items = [], []

        def nid(i, j, k):
            return f"{side}_{i}_{j}_{k}"

        def mass_item(name, category, mass, center, allocations):
            if sum(allocations.values(), Fraction(0)) != mass:
                raise ArithmeticError("Mass-to-load allocation does not conserve weight")
            # Only x/y first moments must agree for vertical gravity. True z is retained.
            for axis in (0, 1):
                if sum((amount*local[node][axis] for node, amount in allocations.items()), Fraction(0)) != mass*center[axis]:
                    raise ArithmeticError("Gravity-load allocation changes its external moment")
            item = {"id": name, "half": side, "category": category, "mass_kg": str(mass),
                    "position_m": _s(center), "local_center_m": _s(center), "installed_center_m": _s(_placed(center, side, span)),
                    "gravity_mass_allocation_kg": {node: str(value) for node, value in allocations.items()}}
            items.append(item)
            mass_items.append(item)
            for node, amount in allocations.items():
                dead[node][2] -= amount*gravity/1000

        def bar(role, a, b):
            start, end = nid(*a), nid(*b)
            delta = [v-u for u, v in zip(local[start], local[end])]
            length = _rational_sqrt(sum((v*v for v in delta), Fraction(0)))
            area = exact(design["member_areas_m2"][role])
            name = f"{side}_{role}_{len(members):03d}"
            mass = density*area*length
            member = {"id": name, "half": side, "start": start, "end": end, "role": role,
                      "length_m": str(length), "area_m2": str(area), "mass_kg": str(mass),
                      "young_modulus_kN_m2": str(exact(design["material"]["young_modulus_kN_m2"]))}
            members.append(member)
            physical_members.append(member)
            center = [(u+v)/2 for u, v in zip(local[start], local[end])]
            mass_item(name+"_mass", "member", mass, center, {start: mass/2, end: mass/2})

        for i in range(n):
            for j in (0, 1):
                for k in (0, 1):
                    bar("longitudinal", (i, j, k), (i+1, j, k))
            for k in (0, 1):
                bar("xy_diagonal", (i, 0, k), (i+1, 1, k))
            for j in (0, 1):
                bar("xz_diagonal", (i, j, 0), (i+1, j, 1))
        for i in range(n+1):
            for k in (0, 1):
                bar("transverse", (i, 0, k), (i, 1, k))
            for j in (0, 1):
                bar("vertical", (i, j, 0), (i, j, 1))
            bar("yz_diagonal", (i, 0, 0), (i, 1, 1))

        def top_allocation(i, mass, center_y):
            positive_fraction = Fraction(1, 2)+center_y/width
            return {nid(i+step, j, 1): mass*(positive_fraction if j else 1-positive_fraction)/2
                    for step in (0, 1) for j in (0, 1)}

        for i in range(n):
            mass = exact(design["mass"]["floor_kg_m2"])*panel*floor_width
            mass_item(f"{side}_floor_{i}", "floor", mass,
                      [(i+Fraction(1, 2))*panel, Fraction(0), height+thickness/2], top_allocation(i, mass, Fraction(0)))
            for edge in (0, 1):
                center_y = (edge-Fraction(1, 2))*floor_width
                mass = exact(design["mass"]["rail_kg_m"])*panel
                center_z = height+thickness+exact(design["mass"]["rail_center_above_walking_surface_m"])
                mass_item(f"{side}_rail_{i}_{edge}", "rail", mass,
                          [(i+Fraction(1, 2))*panel, center_y, center_z], top_allocation(i, mass, center_y))
        for node, position in local.items():
            mass = exact(design["mass"]["joint_kg_per_physical_node"])
            mass_item(node+"_joint", "joint", mass, position, {node: mass})
        total_mass = sum((exact(item["mass_kg"]) for item in items), Fraction(0))
        cg = [sum((exact(item["mass_kg"])*exact(item["local_center_m"][axis]) for item in items), Fraction(0))/total_mass for axis in range(3)]
        halves[side] = {"local_nodes": {name: _s(p) for name, p in local.items()},
                        "global_nodes": {name: _s(p) for name, p in world.items()}, "members": members,
                        "local_mass_items": items, "mass_kg": str(total_mass), "local_center_of_mass_m": _s(cg),
                        "installed_center_of_mass_m": _s(_placed(cg, side, span)),
                        "installed_rotation": [["1", "0", "0"], ["0", "1", "0"], ["0", "0", "1"]] if side == "L" else [["-1", "0", "0"], ["0", "-1", "0"], ["0", "0", "1"]],
                        "installed_translation_m": ["0", "0", "0"] if side == "L" else [str(span), "0", "0"]}
    total_mass = sum((exact(item["mass_kg"]) for item in mass_items), Fraction(0))
    categories = {category: str(sum((exact(item["mass_kg"]) for item in mass_items if item["category"] == category), Fraction(0)))
                  for category in ("member", "floor", "rail", "joint")}
    cg = [sum((exact(item["mass_kg"])*exact(item["installed_center_m"][axis]) for item in mass_items), Fraction(0))/total_mass for axis in range(3)]
    pairs = [{"id": f"center_{j}_{k}", "left": f"L_{n}_{j}_{k}", "right": f"R_{n}_{1-j}_{k}", "level": "lower" if k == 0 else "upper"}
             for j in (0, 1) for k in (0, 1)]
    fixed_mass = sum((exact(item["mass_kg"]) for item in design["base_interfaces"]["fixed_mass_items_per_base"]), Fraction(0))
    return {"schema": "crest-root-assembly-v3", "design_sha256": design_digest(design), "design": deepcopy(design),
            "units": {"length": "m", "force": "kN", "mass": "kg", "moment": "kN m", "stress": "kN/m^2"},
            "dimensions": {"span_m": str(span), "half_length_m": str(span/2), "panel_length_m": str(panel),
                           "truss_width_m": str(width), "truss_depth_m": str(height), "floor_width_m": str(floor_width), "floor_thickness_m": str(thickness),
                           "member_envelope_radius_m": str(exact(design["geometry"]["member_envelope_radius_m"])),
                           "joint_contact_zone_length_m": str(exact(design["geometry"]["joint_contact_zone_length_m"])),
                           "truss_bottom_centerline_elevation_m": str(exact(design["base_interfaces"]["truss_bottom_elevation_m"])),
                           "nominal_lowest_member_surface_elevation_m": str(exact(design["base_interfaces"]["truss_bottom_elevation_m"])-exact(design["geometry"]["member_envelope_radius_m"])),
                           "walking_surface_offset_m": str(height+thickness),
                           "installed_walking_surface_elevation_m": str(exact(design["base_interfaces"]["truss_bottom_elevation_m"])+height+thickness),
                           "transport_walking_surface_elevation_m": str(exact(design["base_interfaces"]["transport_truss_bottom_m"])+height+thickness)},
            "halves": halves, "physical_nodes": physical_nodes, "physical_members": physical_members,
            "physical_node_count": len(physical_nodes), "physical_member_count": len(physical_members),
            "mass_ledger": {"items": mass_items, "bridge_mass_kg": str(total_mass), "bridge_weight_kN": str(total_mass*gravity/1000),
                            "by_category_kg": categories, "installed_bridge_center_of_mass_m": _s(cg),
                            "fixed_mass_per_base_kg": str(fixed_mass), "complete_declared_mass_kg": str(total_mass+2*fixed_mass),
                            "fixed_mass_items_per_base": deepcopy(design["base_interfaces"]["fixed_mass_items_per_base"]),
                            "scope": "All declared bars, floor, rails and node fittings plus two stationary base mass ledgers; omitted hardware is not implicitly massless in reality"},
            "nodal_dead_forces": {name: _s(force) for name, force in dead.items()}, "connector_pairs": pairs,
            "mode_definitions": {
                "independent_half_cantilevers": "No seam ties; all four root-node translations held on each half",
                "joined_clamped": "All four seam-node pairs tied; all root translations still held",
                "joined_service": "All four seam-node pairs tied; four z bearings, two y guides on one global edge, one x restraint",
                "lower_only_join_after_release": "Only the two lower seam pairs tied; permanent service bearings only",
                "premature_release_before_join": "No seam ties; permanent service bearings only",
                "joined_release_midpoint": "All four ties; permanent bearings zero and temporary root translations prescribed at half their service-dead values"},
            "checks_not_performed": list(LIMITATIONS)}


def permanent_restraints(assembly):
    """Seven independent ideal bearing DOFs, all at the outer bottom corners."""
    answer = [{"node": f"{side}_0_{j}_0", "axis": 2, "displacement": "0"} for side in ("L", "R") for j in (0, 1)]
    answer += [{"node": "L_0_0_0", "axis": 1, "displacement": "0"},
               {"node": "R_0_1_0", "axis": 1, "displacement": "0"},
               {"node": "L_0_0_0", "axis": 0, "displacement": "0"}]
    return answer


def _loads(assembly, load_case):
    loads = {name: _v(force) for name, force in assembly["nodal_dead_forces"].items()}
    if load_case == "dead":
        return loads
    if load_case not in ("dead_plus_symmetric_crowd", "dead_plus_asymmetric_crowd"):
        raise DesignError("Unknown load case")
    asymmetric = load_case == "dead_plus_asymmetric_crowd"
    design, dimensions = assembly["design"], assembly["dimensions"]
    n = design["geometry"]["panels_per_half"]
    panel, floor_width, width = (exact(dimensions[key]) for key in ("panel_length_m", "floor_width_m", "truss_width_m"))
    pressure = exact(design["crowd"]["asymmetric_pressure_kN_m2" if asymmetric else "symmetric_pressure_kN_m2"])
    per_panel = pressure*panel*floor_width/(2 if asymmetric else 1)
    center_y = floor_width/4 if asymmetric else Fraction(0)
    for side in ("L", "R"):
        local_center_y = center_y if side == "L" else -center_y
        positive_fraction = Fraction(1, 2)+local_center_y/width
        for i in range(n):
            for step in (0, 1):
                for j in (0, 1):
                    loads[f"{side}_{i+step}_{j}_1"][2] -= per_panel*(positive_fraction if j else 1-positive_fraction)/2
    return loads


def make_case(assembly, mode, load_case="dead", temporary_displacements=None):
    """Build a v2 case and return (case, physical-to-solver alias, physical loads)."""
    if mode not in assembly["mode_definitions"]:
        raise DesignError("Unknown structural mode")
    aliases = {name: name for name in assembly["physical_nodes"]}
    full_join = mode in ("joined_clamped", "joined_service", "joined_release_midpoint")
    for pair in assembly["connector_pairs"]:
        if full_join or (mode == "lower_only_join_after_release" and pair["level"] == "lower"):
            aliases[pair["right"]] = pair["left"]
    nodes = {}
    for physical, name in aliases.items():
        position = assembly["physical_nodes"][physical]
        if name in nodes and nodes[name] != position:
            raise DesignError("A perfect tie was requested between distinct reference points")
        nodes[name] = position
    members = [{"id": member["id"], "start": aliases[member["start"]], "end": aliases[member["end"]],
                "young_modulus": member["young_modulus_kN_m2"], "area": member["area_m2"]}
               for member in assembly["physical_members"]]
    restraints = permanent_restraints(assembly)
    if mode in ("independent_half_cantilevers", "joined_clamped", "joined_release_midpoint"):
        permanent = {(entry["node"], entry["axis"]) for entry in restraints}
        for side in ("L", "R"):
            for j in (0, 1):
                for k in (0, 1):
                    name = f"{side}_0_{j}_{k}"
                    for axis in range(3):
                        if (name, axis) not in permanent:
                            value = "0" if temporary_displacements is None else str(exact(temporary_displacements[name][axis]))
                            restraints.append({"node": name, "axis": axis, "displacement": value})
    if mode == "joined_release_midpoint" and temporary_displacements is None:
        raise DesignError("Release midpoint requires prescribed temporary-root displacements")
    physical_loads = _loads(assembly, load_case)
    loads = [{"node": aliases[node], "components": _s(force)} for node, force in physical_loads.items()]
    return {"id": mode+"/"+load_case, "dimension": 3, "units": {"length": "m", "force": "kN"},
            "nodes": nodes, "members": members, "restraints": restraints, "loads": loads}, aliases, physical_loads


def _resultants(positions, forces, origin):
    total_force = _sum_vectors(forces.values())
    moment = _sum_vectors(_cross([a-b for a, b in zip(_v(positions[node]), origin)], force) for node, force in forces.items())
    return {"force_kN": _s(total_force), "moment_kN_m": _s(moment)}


def _summarize(assembly, mode, load_case, raw, aliases, physical_loads):
    result = {"id": mode+"/"+load_case, "mode": mode, "load_case": load_case,
              "design_sha256": assembly["design_sha256"], "status": raw["status"],
              "physical_member_count": assembly["physical_member_count"], "bridge_mass_kg": assembly["mass_ledger"]["bridge_mass_kg"],
              "solver_node_count": len(raw["node_order"]), "free_dof_count": raw["free_dof_count"],
              "supported_rigidity_rank": raw["supported_rigidity_rank"], "mechanism_dimension": raw["mechanism_dimension"],
              "supported_self_stress_dimension": raw["supported_self_stress_dimension"],
              "physical_to_solver_node": aliases, "physical_applied_nodal_forces_kN": {name: _s(force) for name, force in physical_loads.items()},
              "applied_resultant_about_global_origin": _resultants(assembly["physical_nodes"], physical_loads, [Fraction(0)]*3)}
    if raw["status"] == "mechanism":
        result.update(equilibrium_compatible=raw["equilibrium_compatible"], interpretation=raw["interpretation"],
                      mechanisms=[{"physical_nodal_motion": {node: witness["nodal_motion"][alias] for node, alias in aliases.items()},
                                   "load_virtual_work": witness["load_virtual_work"]} for witness in raw["mechanisms"]])
        return result
    displacements = {node: raw["nodal_displacements"][alias] for node, alias in aliases.items()}
    result.update(physical_nodal_displacements_m=displacements,
                  physical_deformed_coordinates_m={node: _s([a+b for a, b in zip(_v(assembly["physical_nodes"][node]), _v(displacements[node]))]) for node in aliases},
                  physical_support_reactions_kN={node: raw["support_reactions"][alias] for node, alias in aliases.items()},
                  members=raw["members"], strain_energy_kN_m=raw["strain_energy"],
                  total_force_equilibrium_residual_kN=raw["total_force_equilibrium_residual"],
                  reference_moment_equilibrium_residual_kN_m=raw["reference_geometry_moment_equilibrium_residual"],
                  interpretation=raw["interpretation"])
    span = exact(assembly["dimensions"]["span_m"])
    bases = {}
    for side, origin in (("L", [Fraction(0)]*3), ("R", [span, Fraction(0), Fraction(0)])):
        forces = {node: _v(raw["support_reactions"][aliases[node]]) for node in aliases if node.startswith(side+"_0_")}
        resultant = _resultants(assembly["physical_nodes"], forces, origin)
        bases[side] = {"reference_origin_m": _s(origin), "reaction_on_bridge": resultant,
                       "bridge_on_base": {key: _s([-value for value in _v(vector)]) for key, vector in resultant.items()}}
    result["bases"] = bases
    gradients = {node: [Fraction(0)]*3 for node in aliases}
    by_id = {member["id"]: member for member in raw["members"]}
    for member in assembly["physical_members"]:
        force = exact(by_id[member["id"]]["axial_force"])
        delta = [b-a for a, b in zip(_v(assembly["physical_nodes"][member["start"]]), _v(assembly["physical_nodes"][member["end"]]))]
        contribution = [force*d/exact(member["length_m"]) for d in delta]
        for node, sign in ((member["start"], -1), (member["end"], 1)):
            gradients[node] = [a+sign*b for a, b in zip(gradients[node], contribution)]
    connectors = []
    for pair in assembly["connector_pairs"]:
        left, right = pair["left"], pair["right"]
        tied = aliases[left] == aliases[right]
        residuals = {node: [a-b for a, b in zip(gradients[node], physical_loads[node])] for node in (left, right)}
        connectors.append({**pair, "tied": tied, "force_on_left_kN": _s(residuals[left]), "force_on_right_kN": _s(residuals[right]),
                           "force_pair_residual_kN": _s([a+b for a, b in zip(residuals[left], residuals[right])]),
                           "relative_displacement_right_minus_left_m": _s([b-a for a, b in zip(_v(displacements[left]), _v(displacements[right]))])})
        if any(a+b for a, b in zip(residuals[left], residuals[right])) or (not tied and any(residuals[left]+residuals[right])):
            raise ArithmeticError("Seam free-node or tie-force equilibrium failed")
    result["connectors"] = connectors
    max_node = max(displacements, key=lambda node: abs(exact(displacements[node][2])))
    max_member = max(raw["members"], key=lambda member: abs(exact(member["axial_stress"])))
    result["scalar_extrema"] = {"largest_absolute_vertical_displacement": {"node": max_node, "signed_value_m": displacements[max_node][2]},
                                "largest_absolute_axial_stress": {"member": max_member["id"], "signed_value_kN_m2": max_member["axial_stress"]},
                                "interpretation": "Response extrema only; no allowable-stress, deformation or capacity criterion applied"}
    return result


def _solve(assembly, mode, load_case="dead", temporary_displacements=None):
    case, aliases, physical_loads = make_case(assembly, mode, load_case, temporary_displacements)
    raw = _truss.solve_truss(case)
    return _summarize(assembly, mode, load_case, raw, aliases, physical_loads)


def run_assembly(design):
    """Run fixed physical inventory through explicit support/join/load states."""
    assembly = build_assembly(design)
    cases = [_solve(assembly, mode, load_case) for mode, load_case in (
        ("independent_half_cantilevers", "dead"), ("joined_clamped", "dead"), ("joined_service", "dead"),
        ("joined_service", "dead_plus_symmetric_crowd"), ("joined_service", "dead_plus_asymmetric_crowd"),
        ("lower_only_join_after_release", "dead"), ("premature_release_before_join", "dead"))]
    by_id = {case["id"]: case for case in cases}
    clamped, service = by_id["joined_clamped/dead"], by_id["joined_service/dead"]
    if clamped["status"] != "unique_linear_solution" or service["status"] != "unique_linear_solution":
        raise ArithmeticError("Reference fully joined modes must have unique linear solutions")
    midpoint_values = {node: _s([exact(value)/2 for value in vector]) for node, vector in service["physical_nodal_displacements_m"].items()}
    midpoint = _solve(assembly, "joined_release_midpoint", "dead", midpoint_values)
    cases.append(midpoint)
    midpoint_matches = all(exact(value) == (exact(clamped["physical_nodal_displacements_m"][node][axis])+exact(service["physical_nodal_displacements_m"][node][axis]))/2
                           for node, vector in midpoint["physical_nodal_displacements_m"].items() for axis, value in enumerate(vector))
    reaction_matches = all(exact(value) == (exact(clamped["physical_support_reactions_kN"][node][axis])+exact(service["physical_support_reactions_kN"][node][axis]))/2
                          for node, vector in midpoint["physical_support_reactions_kN"].items() for axis, value in enumerate(vector))
    force_matches = all(exact(m["axial_force"]) == (exact(c["axial_force"])+exact(s["axial_force"]))/2
                       for m, c, s in zip(midpoint["members"], clamped["members"], service["members"]))
    permanent = {(r["node"], r["axis"]) for r in permanent_restraints(assembly)}
    temporary = [{"node": f"{side}_0_{j}_{k}", "axis": axis} for side in ("L", "R") for j in (0, 1) for k in (0, 1) for axis in range(3)
                 if (f"{side}_0_{j}_{k}", axis) not in permanent]
    zero_temporary = all(exact(service["physical_support_reactions_kN"][r["node"]][r["axis"]]) == 0 for r in temporary)
    if not all((midpoint_matches, reaction_matches, force_matches, zero_temporary)):
        raise ArithmeticError("Affine commanded-support release identity failed")
    prejoin = by_id["independent_half_cantilevers/dead"]
    deformed_local = {}
    span = exact(assembly["dimensions"]["span_m"])
    for side in ("L", "R"):
        deformed_local[side] = {node: _s(_placed(_v(prejoin["physical_deformed_coordinates_m"][node]), side, span))
                                for node in assembly["halves"][side]["local_nodes"]}
    invalid_modes = {"lower_only_join_after_release", "premature_release_before_join"}
    expected_ids = {member["id"] for member in assembly["physical_members"]}
    checks = {
        "mass_conserved": all(case["physical_member_count"] == assembly["physical_member_count"]
                              and case["bridge_mass_kg"] == assembly["mass_ledger"]["bridge_mass_kg"]
                              and (case["status"] != "unique_linear_solution" or {member["id"] for member in case["members"]} == expected_ids)
                              for case in cases),
        "all_expected_models_solved": all(case["status"] == "unique_linear_solution" for case in cases if case["mode"] not in invalid_modes),
        "invalid_modes_have_mechanism": all(case["status"] == "mechanism" for case in cases if case["mode"] in invalid_modes),
        "controlled_release_verified": all((midpoint_matches, reaction_matches, force_matches, zero_temporary)),
    }
    return {"schema": "crest-root-assembly-report-v3", "schema_version": 1, "physical_approval": False,
            "design_sha256": assembly["design_sha256"], "geometry_and_mass_source_sha256": design_digest(assembly), "checks": checks,
            "source_sha256": {"assembly.py": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
                              "v2/truss.py": hashlib.sha256(TRUSS_PATH.read_bytes()).hexdigest()},
            "assembly": assembly, "cases": cases,
            "release_certificate": {"model": "quasi-static prescribed support translations in one fixed linear stiffness model",
                                    "parameter": "0 <= lambda <= 1", "displacement_identity": "u(lambda)=(1-lambda)*u_clamped+lambda*u_service",
                                    "temporary_root_dofs": temporary, "command": "u_temporary(lambda)=lambda*u_service_temporary; permanent bearing displacements remain zero",
                                    "actual_midpoint_solve_case": midpoint["id"], "midpoint_displacements_match": midpoint_matches,
                                    "midpoint_member_forces_match": force_matches, "midpoint_reactions_match": reaction_matches,
                                    "temporary_reactions_zero_at_lambda_one": zero_temporary,
                                    "endpoint_permanent_reactions_from_case": service["id"],
                                    "scope": "Algebraic equilibrium path for commanded ideal restraints, not instant release or validated actuator/control hardware"},
            "capture_reference": {"source_case": prejoin["id"], "deformed_local_nodes_by_half": deformed_local,
                                  "scope": "Small-displacement dead-load cantilever shapes under ideal root clamps; capture tolerance, contact and actuation are separate"},
            "checks_not_performed": list(LIMITATIONS)}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--design", type=Path, default=DEFAULT_DESIGN)
    parser.add_argument("--geometry-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    design = load_design(args.design)
    result = build_assembly(design) if args.geometry_only else run_assembly(design)
    rendered = json.dumps(result, indent=2)+"\n"
    if args.output:
        args.output.write_text(rendered, encoding="utf-8")
    else:
        print(rendered, end="")
