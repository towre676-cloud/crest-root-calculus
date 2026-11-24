"""Independent arithmetic checks on the saved v3 assembly and linear results.

This verifier does not call the truss solver. It reconstructs external loads,
bar constitutive forces, nodal equilibrium and energy from elementary formulas;
certifies ranks by modular lower bounds plus exact nullspace upper bounds; and
checks a hand-derived rigid-half folding velocity for the invalid lower join.
Default operation needs only Python's standard library. --replay additionally
regenerates the report using the production solver and compares every field.
No check here establishes physical capacity or operation permission.
"""

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import json
from pathlib import Path

import assembly

HERE = Path(__file__).resolve().parent
REPORT = HERE / "assembly-report.json"
PRIME = 1000003


def require(condition, explanation):
    if not condition:
        raise ValueError(explanation)


def vec(values):
    return [F(value) for value in values]


def cross(a, b):
    return [a[1]*b[2]-a[2]*b[1], a[2]*b[0]-a[0]*b[2], a[0]*b[1]-a[1]*b[0]]


def dot(a, b):
    return sum((x*y for x, y in zip(a, b)), F(0))


def add_to(target, name, value):
    target[name] = [a+b for a, b in zip(target[name], value)]


def resultants(nodes, forces, origin=(F(0), F(0), F(0))):
    total, moment = [F(0)]*3, [F(0)]*3
    for node, force in forces.items():
        total = [a+b for a, b in zip(total, force)]
        lever = [a-b for a, b in zip(vec(nodes[node]), origin)]
        moment = [a+b for a, b in zip(moment, cross(lever, force))]
    return total, moment


def modular_rank(matrix):
    """Nonzero minors modulo a prime certify rational rank lower bounds."""
    rows = [[(F(value).numerator % PRIME)*pow(F(value).denominator % PRIME, -1, PRIME) % PRIME for value in row] for row in matrix]
    if not rows:
        return 0
    height, width = len(rows), len(rows[0])
    rank = 0
    for column in range(width):
        found = next((r for r in range(rank, height) if rows[r][column]), None)
        if found is None:
            continue
        rows[rank], rows[found] = rows[found], rows[rank]
        inverse = pow(rows[rank][column], -1, PRIME)
        pivot = rows[rank]
        for r in range(rank+1, height):
            multiplier = rows[r][column]*inverse % PRIME
            if multiplier:
                rows[r] = [(a-multiplier*b) % PRIME for a, b in zip(rows[r], pivot)]
        rank += 1
        if rank == min(height, width):
            break
    return rank


def independent_loads(a, case):
    """Closed-form tributary forces: bar endpoints, joints, floor and two rails."""
    d, dims = a["design"], a["dimensions"]
    nodes = a["physical_nodes"]
    loads = {node: [F(0)]*3 for node in nodes}
    g = F(d["gravity_m_s2"])/1000
    for member in a["physical_members"]:
        force = -F(member["mass_kg"])*g/2
        for node in (member["start"], member["end"]):
            loads[node][2] += force
    n = d["geometry"]["panels_per_half"]
    panel, floor, width = (F(dims[key]) for key in ("panel_length_m", "floor_width_m", "truss_width_m"))
    # Paired rails are symmetric. Each of four panel top nodes gets one fourth
    # of floor mass and one fourth of the two rail masses, independent of rail z.
    nonstructural_mass_per_panel_corner = F(d["mass"]["floor_kg_m2"])*panel*floor/4+F(d["mass"]["rail_kg_m"])*panel/2
    for node in nodes:
        side, station, j, k = node.split("_")
        station, j, k = int(station), int(j), int(k)
        loads[node][2] -= F(d["mass"]["joint_kg_per_physical_node"])*g
        adjacent = 1 if station in (0, n) else 2
        if k:
            loads[node][2] -= adjacent*nonstructural_mass_per_panel_corner*g
            if case["load_case"] == "dead_plus_symmetric_crowd":
                loads[node][2] -= adjacent*F(d["crowd"]["symmetric_pressure_kN_m2"])*panel*floor/4
            elif case["load_case"] == "dead_plus_asymmetric_crowd":
                positive_edge = F(nodes[node][1]) > 0
                edge_share = F(1, 2)+(1 if positive_edge else -1)*floor/(4*width)
                loads[node][2] -= adjacent*F(d["crowd"]["asymmetric_pressure_kN_m2"])*panel*floor*edge_share/4
            else:
                require(case["load_case"] == "dead", "Unknown load case in report")
    return loads


def independent_restraints(a, mode):
    restraints = {(f"{side}_0_{j}_0", 2) for side in ("L", "R") for j in (0, 1)}
    restraints.update({("L_0_0_0", 0), ("L_0_0_0", 1), ("R_0_1_0", 1)})
    if mode in ("independent_half_cantilevers", "joined_clamped", "joined_release_midpoint"):
        restraints = {(f"{side}_0_{j}_{k}", axis) for side in ("L", "R") for j in (0, 1) for k in (0, 1) for axis in range(3)}
    return restraints


def check_geometry(a, design):
    require(a == assembly.build_assembly(design), "Saved assembly differs from current design-to-assembly builder")
    n = design["geometry"]["panels_per_half"]
    nodes, bars = a["physical_nodes"], a["physical_members"]
    require(len(nodes) == 8*(n+1), "Physical node inventory changed")
    require(len(bars) == 2*(8*n+5*(n+1)), "Physical member inventory changed")
    require(len({m["id"] for m in bars}) == len(bars), "Duplicate physical member")
    p, w, h = (F(a["dimensions"][key]) for key in ("panel_length_m", "truss_width_m", "truss_depth_m"))
    span, density = F(design["geometry"]["span_m"]), F(design["material"]["density_kg_m3"])
    require(p == span/(2*n), "Panel span inconsistency")
    role_lengths_squared = {"longitudinal": p*p, "transverse": w*w, "vertical": h*h,
                            "xy_diagonal": p*p+w*w, "xz_diagonal": p*p+h*h, "yz_diagonal": w*w+h*h}
    totals = {key: F(0) for key in ("member", "floor", "rail", "joint")}
    for member in bars:
        delta = [b-a for a, b in zip(vec(nodes[member["start"]]), vec(nodes[member["end"]]))]
        L = F(member["length_m"])
        require(L*L == dot(delta, delta) == role_lengths_squared[member["role"]], "Bar metric/role inconsistency")
        require(F(member["mass_kg"]) == density*F(member["area_m2"])*L, "Member self-weight mass inconsistent")
        totals["member"] += F(member["mass_kg"])
    totals["floor"] = span*F(design["geometry"]["floor_width_m"])*F(design["mass"]["floor_kg_m2"])
    totals["rail"] = 2*span*F(design["mass"]["rail_kg_m"])
    totals["joint"] = len(nodes)*F(design["mass"]["joint_kg_per_physical_node"])
    require({k: F(v) for k, v in a["mass_ledger"]["by_category_kg"].items()} == totals, "Independent category mass mismatch")
    require(sum(totals.values()) == F(a["mass_ledger"]["bridge_mass_kg"]), "Total mass mismatch")
    true_first_moment = [F(0)]*3
    for item in a["mass_ledger"]["items"]:
        mass = F(item["mass_kg"])
        require(sum(map(F, item["gravity_mass_allocation_kg"].values())) == mass, "An item loses gravity load")
        position = vec(item["installed_center_m"])
        for axis in range(3):
            true_first_moment[axis] += mass*position[axis]
        for axis in (0, 1):
            moment = sum((F(value)*F(nodes[node][axis]) for node, value in item["gravity_mass_allocation_kg"].items()), F(0))
            require(moment == mass*position[axis], "Lumped gravity changes an x/y first moment")
    total_mass = sum(totals.values())
    require([value/total_mass for value in true_first_moment] == vec(a["mass_ledger"]["installed_bridge_center_of_mass_m"]), "True mass centroid mismatch")
    require(true_first_moment[0] == total_mass*span/2 and true_first_moment[1] == 0, "Expected whole-bridge symmetry lost")
    for pair in a["connector_pairs"]:
        require(nodes[pair["left"]] == nodes[pair["right"]], "Reference seam nodes do not coincide")
    require(F(a["dimensions"]["walking_surface_offset_m"]) == h+F(design["geometry"]["floor_thickness_m"]), "Access rise omitted truss depth")
    return {"physical_nodes": len(nodes), "physical_members": len(bars), "mass_items": len(a["mass_ledger"]["items"]), "category_mass_checks": 4}


def check_rank_and_witnesses(a, case):
    aliases = case["physical_to_solver_node"]
    require(set(aliases) == set(a["physical_nodes"]), "Missing physical node alias")
    tied_pairs = {(p["left"], p["right"]) for p in a["connector_pairs"]
                  if case["mode"] in ("joined_clamped", "joined_service", "joined_release_midpoint")
                  or (case["mode"] == "lower_only_join_after_release" and p["level"] == "lower")}
    expected = {node: node for node in a["physical_nodes"]}
    for left, right in tied_pairs:
        expected[right] = left
    require(aliases == expected, "Mode seam aliases were changed")
    order = list(dict.fromkeys(aliases.values()))
    constrained = independent_restraints(a, case["mode"])
    free = [(node, axis) for node in order for axis in range(3) if (node, axis) not in constrained]
    indices = {entry: i for i, entry in enumerate(free)}
    rows = []
    for bar in a["physical_members"]:
        delta = [b-a for a, b in zip(vec(a["physical_nodes"][bar["start"]]), vec(a["physical_nodes"][bar["end"]]))]
        row = [F(0)]*len(free)
        for node, sign in ((aliases[bar["start"]], -1), (aliases[bar["end"]], 1)):
            for axis in range(3):
                if (node, axis) in indices:
                    row[indices[node, axis]] += sign*delta[axis]
        rows.append(row)
    lower_rank = modular_rank(rows)
    require(lower_rank == case["supported_rigidity_rank"], "Modular rank certificate differs from reported rank")
    require(len(free) == case["free_dof_count"], "Free DOF count mismatch")
    if case["status"] == "unique_linear_solution":
        require(lower_rank == len(free) and case["mechanism_dimension"] == 0, "Unique solution lacks full-rank certificate")
        return {"rank": lower_rank, "free_dofs": len(free), "witnesses": 0}
    require(case["status"] == "mechanism", "Unexpected case status")
    basis = []
    compatible = True
    loads = independent_loads(a, case)
    for witness in case["mechanisms"]:
        motion = {node: vec(v) for node, v in witness["physical_nodal_motion"].items()}
        require(set(motion) == set(a["physical_nodes"]), "Incomplete mechanism witness")
        require(any(any(vector) for vector in motion.values()), "Zero mechanism witness")
        for left, right in tied_pairs:
            require(motion[left] == motion[right], "Mechanism violates seam tie")
        require(all(motion[node][axis] == 0 for node, axis in constrained), "Mechanism violates a support")
        for bar in a["physical_members"]:
            delta = [b-a for a, b in zip(vec(a["physical_nodes"][bar["start"]]), vec(a["physical_nodes"][bar["end"]]))]
            du = [b-a for a, b in zip(motion[bar["start"]], motion[bar["end"]])]
            require(dot(delta, du) == 0, "Mechanism changes first-order member length")
        work = sum((dot(loads[node], motion[node]) for node in loads), F(0))
        require(work == F(witness["load_virtual_work"]), "Mechanism load virtual work mismatch")
        compatible = compatible and work == 0
        basis.append([motion[node][axis] for node, axis in free])
    require(modular_rank(basis) == len(basis) == case["mechanism_dimension"], "Mechanism witnesses are dependent or miscounted")
    require(lower_rank+len(basis) == len(free), "Rank lower bound and nullity upper bound do not meet")
    require(case["equilibrium_compatible"] == compatible, "Equilibrium compatibility flag mismatch")
    return {"rank": lower_rank, "free_dofs": len(free), "witnesses": len(basis)}


def check_unique_case(a, case, service=None):
    nodes, bars = a["physical_nodes"], a["physical_members"]
    loads = independent_loads(a, case)
    require({node: vec(force) for node, force in case["physical_applied_nodal_forces_kN"].items()} == loads, "Declared physical loading differs from independent tributary loads")
    u = {node: vec(value) for node, value in case["physical_nodal_displacements_m"].items()}
    reactions = {node: vec(value) for node, value in case["physical_support_reactions_kN"].items()}
    require(set(u) == set(reactions) == set(nodes), "Missing displacement/reaction node")
    prescribed = independent_restraints(a, case["mode"])
    permanent = independent_restraints(a, "joined_service")
    for node in nodes:
        require(vec(case["physical_deformed_coordinates_m"][node]) == [a+b for a, b in zip(vec(nodes[node]), u[node])], "Deformed coordinate mismatch")
        for axis in range(3):
            if (node, axis) not in prescribed:
                require(reactions[node][axis] == 0, "Reaction at an unrestrained DOF")
            else:
                target = F(0)
                if case["mode"] == "joined_release_midpoint" and (node, axis) not in permanent:
                    require(service is not None, "Missing service endpoint for release check")
                    target = F(service["physical_nodal_displacements_m"][node][axis])/2
                require(u[node][axis] == target, "Incorrect prescribed support translation")
    for physical, alias in case["physical_to_solver_node"].items():
        require(u[physical] == u[alias], "Tied nodes have different displacements")
    results = {member["id"]: member for member in case["members"]}
    require(len(results) == len(bars) == len(case["members"]) and set(results) == {bar["id"] for bar in bars}, "Physical bar inventory changed in solve")
    gradient, energy = {node: [F(0)]*3 for node in nodes}, F(0)
    for bar in bars:
        answer = results[bar["id"]]
        delta = [b-a for a, b in zip(vec(nodes[bar["start"]]), vec(nodes[bar["end"]]))]
        du = [b-a for a, b in zip(u[bar["start"]], u[bar["end"]])]
        length, EA = F(bar["length_m"]), F(bar["young_modulus_kN_m2"])*F(bar["area_m2"])
        extension = dot(delta, du)/length
        force = EA*extension/length
        member_energy = EA*extension*extension/(2*length)
        require(F(answer["extension"]) == extension and F(answer["axial_force"]) == force, "Independent bar constitutive equation failed")
        require(F(answer["axial_stress"]) == force/F(bar["area_m2"]), "Stress conversion failed")
        require(F(answer["strain_energy"]) == member_energy, "Bar energy mismatch")
        require(answer["force_sign"] == ("tension" if force > 0 else "compression" if force < 0 else "zero"), "Member sign label mismatch")
        energy += member_energy
        add_to(gradient, bar["start"], [-force*d/length for d in delta])
        add_to(gradient, bar["end"], [force*d/length for d in delta])
    require(energy == F(case["strain_energy_kN_m"]), "Total strain energy mismatch")
    residual = {node: [g-f-r for g, f, r in zip(gradient[node], loads[node], reactions[node])] for node in nodes}
    aggregated = {alias: [F(0)]*3 for alias in case["physical_to_solver_node"].values()}
    for node, alias in case["physical_to_solver_node"].items():
        add_to(aggregated, alias, residual[node])
    require(not any(any(value) for value in aggregated.values()), "Independent nodal equilibrium failed")
    work = sum((dot(u[node], [f+r for f, r in zip(loads[node], reactions[node])]) for node in nodes), F(0))
    require(work == 2*energy, "Applied/reaction work does not equal twice energy")
    total, moment = resultants(nodes, {node: [f+r for f, r in zip(loads[node], reactions[node])] for node in nodes})
    require(total == moment == [F(0)]*3, "Whole-bridge force/moment equilibrium failed")
    require(vec(case["total_force_equilibrium_residual_kN"]) == vec(case["reference_moment_equilibrium_residual_kN_m"]) == [F(0)]*3, "Equilibrium status falsely reported")
    applied_force, applied_moment = resultants(nodes, loads)
    require(vec(case["applied_resultant_about_global_origin"]["force_kN"]) == applied_force
            and vec(case["applied_resultant_about_global_origin"]["moment_kN_m"]) == applied_moment, "Applied resultant mismatch")
    span = F(a["dimensions"]["span_m"])
    for side in ("L", "R"):
        origin = [F(0) if side == "L" else span, F(0), F(0)]
        roots = {node: force for node, force in reactions.items() if node.startswith(side+"_0_")}
        force, moment = resultants(nodes, roots, origin)
        base = case["bases"][side]
        require(vec(base["reference_origin_m"]) == origin, "Base resultant uses wrong origin")
        require(vec(base["reaction_on_bridge"]["force_kN"]) == force and vec(base["reaction_on_bridge"]["moment_kN_m"]) == moment, "Base reaction resultant mismatch")
        require(vec(base["bridge_on_base"]["force_kN"]) == [-v for v in force]
                and vec(base["bridge_on_base"]["moment_kN_m"]) == [-v for v in moment], "Base/bridge action-reaction sign mismatch")
        if case["mode"] == "independent_half_cantilevers":
            half_weight = F(a["mass_ledger"]["bridge_weight_kN"])/2
            require(force == [F(0), F(0), half_weight], "Analytic half weight reaction failed")
            require(moment == [F(0), (-1 if side == "L" else 1)*half_weight*span/4, F(0)], "Analytic cantilever gravity moment failed")
    require(len(case["connectors"]) == 4, "Four seam pairs must be reported")
    for pair in case["connectors"]:
        left, right = pair["left"], pair["right"]
        require(vec(pair["force_on_left_kN"]) == residual[left] and vec(pair["force_on_right_kN"]) == residual[right], "Connector demand does not balance its physical half")
        require(vec(pair["force_pair_residual_kN"]) == [F(0)]*3, "Connector action-reaction does not close")
        require(vec(pair["relative_displacement_right_minus_left_m"]) == [b-a for a, b in zip(u[left], u[right])], "Seam displacement mismatch")
        require(pair["tied"] == (case["physical_to_solver_node"][left] == case["physical_to_solver_node"][right]), "False connector status")
    return {"constitutive_bars": len(bars), "physical_node_equilibrium_checks": len(nodes), "root_resultants": 2, "connector_pairs": 4}


def check_release(report):
    cases = {case["id"]: case for case in report["cases"]}
    c, s, m = (cases[name+"/dead"] for name in ("joined_clamped", "joined_service", "joined_release_midpoint"))
    for field in ("physical_nodal_displacements_m", "physical_support_reactions_kN"):
        for node in m[field]:
            require(vec(m[field][node]) == [(a+b)/2 for a, b in zip(vec(c[field][node]), vec(s[field][node]))], "Release midpoint is not affine")
    by_case = [{bar["id"]: F(bar["axial_force"]) for bar in case["members"]} for case in (c, s, m)]
    require(all(by_case[2][key] == (by_case[0][key]+by_case[1][key])/2 for key in by_case[0]), "Release member-force identity failed")
    permanent = independent_restraints(report["assembly"], "joined_service")
    temporary = independent_restraints(report["assembly"], "joined_clamped")-permanent
    certificate = report["release_certificate"]
    require({(entry["node"], entry["axis"]) for entry in certificate["temporary_root_dofs"]} == temporary, "Incorrect temporary support set")
    require(all(F(s["physical_support_reactions_kN"][node][axis]) == 0 for node, axis in temporary), "Final temporary reactions are nonzero")
    for flag in ("midpoint_displacements_match", "midpoint_member_forces_match", "midpoint_reactions_match", "temporary_reactions_zero_at_lambda_one"):
        require(certificate[flag] is True, "Incorrect release status flag")
    return {"temporary_root_dofs": len(temporary), "actual_solved_midpoint": True, "member_affinity_checks": len(by_case[0])}


def check_analytic_invalid_mode(a, case):
    """Independent rigid-half pitch field: finite folding is not asserted."""
    span = F(a["dimensions"]["span_m"])
    motion = {}
    for node, xyz in a["physical_nodes"].items():
        x, y, z = vec(xyz)
        motion[node] = [z, F(0), -x] if node.startswith("L_") else [-z, F(0), x-span]
    for node, axis in independent_restraints(a, "lower_only_join_after_release"):
        require(motion[node][axis] == 0, "Analytic pitch violates bearing")
    for pair in a["connector_pairs"]:
        if pair["level"] == "lower":
            require(motion[pair["left"]] == motion[pair["right"]], "Analytic pitch violates lower tie")
    for bar in a["physical_members"]:
        delta = [b-a for a, b in zip(vec(a["physical_nodes"][bar["start"]]), vec(a["physical_nodes"][bar["end"]]))]
        du = [b-a for a, b in zip(motion[bar["start"]], motion[bar["end"]])]
        require(dot(delta, du) == 0, "Analytic pitch changes bar length")
    loads = independent_loads(a, case)
    work = sum((dot(loads[node], motion[node]) for node in loads), F(0))
    require(work == F(a["mass_ledger"]["bridge_weight_kN"])*span/4 > 0, "Analytic load incompatibility proof failed")
    return {"independent_rigid_half_pitch_witness": True, "exact_virtual_work": str(work), "meaning": "First-order mechanism with nonzero gravity virtual work; no finite folding path claim"}


def check_report(report, design):
    require(report["schema"] == "crest-root-assembly-report-v3" and report["schema_version"] == 1 and report["physical_approval"] is False, "Incorrect scope/schema flags")
    require(report["design_sha256"] == assembly.design_digest(design), "Design digest mismatch")
    for name, path in (("assembly.py", HERE/"assembly.py"), ("v2/truss.py", assembly.TRUSS_PATH)):
        require(report["source_sha256"][name] == hashlib.sha256(path.read_bytes()).hexdigest(), "Source digest mismatch: "+name)
    a = report["assembly"]
    require(report["geometry_and_mass_source_sha256"] == assembly.design_digest(a), "Geometry/mass digest mismatch")
    geometry = check_geometry(a, design)
    cases = {case["id"]: case for case in report["cases"]}
    expected = {"independent_half_cantilevers/dead", "joined_clamped/dead", "joined_service/dead", "joined_service/dead_plus_symmetric_crowd", "joined_service/dead_plus_asymmetric_crowd", "lower_only_join_after_release/dead", "premature_release_before_join/dead", "joined_release_midpoint/dead"}
    require(set(cases) == expected and len(cases) == len(report["cases"]), "A required case is missing, duplicated or unexpected")
    rank_checks, solution_checks = {}, {}
    for case in cases.values():
        require(case["design_sha256"] == report["design_sha256"], "A case uses another design")
        require(case["physical_member_count"] == len(a["physical_members"]) and F(case["bridge_mass_kg"]) == F(a["mass_ledger"]["bridge_mass_kg"]), "Mode change lost physical mass/members")
        require({node: vec(force) for node, force in case["physical_applied_nodal_forces_kN"].items()} == independent_loads(a, case), "Case load prescription changed")
        rank_checks[case["id"]] = check_rank_and_witnesses(a, case)
        if case["status"] == "unique_linear_solution":
            solution_checks[case["id"]] = check_unique_case(a, case, cases["joined_service/dead"])
    require(cases["lower_only_join_after_release/dead"]["mechanism_dimension"] == 1, "Lower-only reference should have one mechanism")
    require(cases["premature_release_before_join/dead"]["mechanism_dimension"] == 5, "Unjoined reference should have five mechanisms")
    require(not cases["lower_only_join_after_release/dead"]["equilibrium_compatible"] and not cases["premature_release_before_join/dead"]["equilibrium_compatible"], "Invalid gravity states falsely compatible")
    release = check_release(report)
    invalid = check_analytic_invalid_mode(a, cases["lower_only_join_after_release/dead"])
    capture = report["capture_reference"]
    require(capture["source_case"] == "independent_half_cantilevers/dead", "Capture uses wrong structural state")
    span = F(a["dimensions"]["span_m"])
    for side in ("L", "R"):
        for node, local in capture["deformed_local_nodes_by_half"][side].items():
            x, y, z = vec(local)
            world = [x, y, z] if side == "L" else [span-x, -y, z]
            require(world == vec(cases[capture["source_case"]]["physical_deformed_coordinates_m"][node]), "Capture local/world mapping error")
    require(report["checks"] == {"mass_conserved": True, "all_expected_models_solved": True, "invalid_modes_have_mechanism": True, "controlled_release_verified": True}, "Top-level check flags disagree")
    return {"geometry": geometry, "rank_certificates": rank_checks, "linear_solution_certificates": solution_checks,
            "release": release, "independent_invalid_mode": invalid, "capture_mapping": True}


def negative_controls(report, design):
    rejected = []

    def rejects(name, operation):
        try:
            operation()
        except (ValueError, ArithmeticError, KeyError):
            rejected.append(name)
        else:
            raise ValueError("Negative control unexpectedly accepted: "+name)

    for name, field, value in (("float_span", ("geometry", "span_m"), 23.0),
                               ("negative_area", ("member_areas_m2", "longitudinal"), "-0.003"),
                               ("zero_gravity", ("gravity_m_s2",), "0"),
                               ("floor_outside_truss", ("geometry", "floor_width_m"), "3"),
                               ("invalid_force_bounds", ("motion", "horizontal_fixture_force_bounds_kN", "x"), ["2", "1"])):
        bad = deepcopy(design)
        target = bad
        for key in field[:-1]:
            target = target[key]
        target[field[-1]] = value
        rejects(name, lambda bad=bad: assembly.build_assembly(bad))
    bad = deepcopy(design)
    bad["geometry"]["brick_ratio"]["width"] = 253
    rejects("irrational_bar_metric", lambda: assembly.build_assembly(bad))
    cases = {case["id"]: case for case in report["cases"]}
    bad_case = deepcopy(cases["joined_service/dead"])
    bad_case["members"][0]["axial_force"] = str(F(bad_case["members"][0]["axial_force"])+1)
    rejects("tampered_member_force", lambda: check_unique_case(report["assembly"], bad_case))
    bad_case = deepcopy(cases["premature_release_before_join/dead"])
    bad_case["mechanisms"][0]["physical_nodal_motion"]["L_0_0_0"][2] = "1"
    rejects("tampered_mechanism_support", lambda: check_rank_and_witnesses(report["assembly"], bad_case))
    bad_case = deepcopy(cases["joined_service/dead"])
    bad_case["physical_applied_nodal_forces_kN"]["L_1_0_1"][2] = "0"
    rejects("tampered_gravity_load", lambda: check_unique_case(report["assembly"], bad_case))
    bad_report = deepcopy(report)
    next(c for c in bad_report["cases"] if c["mode"] == "joined_release_midpoint")["physical_nodal_displacements_m"]["L_1_0_0"][2] = "0"
    rejects("tampered_release_midpoint", lambda: check_release(bad_report))
    bad_report = deepcopy(report)
    bad_report["design_sha256"] = "0"*64
    rejects("tampered_design_digest", lambda: check_report(bad_report, design))
    bad_report = deepcopy(report)
    bad_report["geometry_and_mass_source_sha256"] = "0"*64
    rejects("tampered_geometry_digest", lambda: check_report(bad_report, design))
    return rejected


def verify_all(report_path=REPORT, design_path=assembly.DEFAULT_DESIGN, replay=False):
    report_path, design_path = Path(report_path), Path(design_path)
    report, design = json.loads(report_path.read_text()), assembly.load_design(design_path)
    certificates = check_report(report, design)
    controls = negative_controls(report, design)
    if replay:
        require(assembly.run_assembly(design) == report, "Full deterministic replay differs from saved report")
    return {"schema": "crest-root-assembly-verification-v3", "checks_passed": True, "physical_approval": False,
            "report_sha256": hashlib.sha256(report_path.read_bytes()).hexdigest(),
            "design_sha256": assembly.design_digest(design),
            "verifier_sha256": hashlib.sha256(Path(__file__).read_bytes()).hexdigest(),
            "method": "Independent exact mechanics identities, modular rank lower bounds plus exact nullspace witnesses, analytic half-cantilever and invalid-join controls",
            "certificates": certificates, "negative_controls_rejected": controls, "negative_control_count": len(controls),
            "full_solver_replay": replay,
            "scope": "Arithmetic consistency and completeness for this declared linear reference assembly; no capacity, stability beyond the linear model, or operation approval"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--report", type=Path, default=REPORT)
    parser.add_argument("--design", type=Path, default=assembly.DEFAULT_DESIGN)
    parser.add_argument("--replay", action="store_true")
    args = parser.parse_args()
    print(json.dumps(verify_all(args.report, args.design, args.replay), indent=2))
