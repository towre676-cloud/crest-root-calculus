"""Generate two declared, synthetic through-truss-plane studies.

These 24 m and 12 m lengths are illustrative geometry, NOT dimensions measured
from the Tillamook atlas. E and area are arbitrary elastic-analysis inputs,
not a member selection. The two planes are assumed to share vertical deck
loading equally; their spatial connections, lateral behavior and torsion are
absent. This script makes no strength, buckling, joint, mast, vibration or
operational-capacity claim. Run with --write to save cases and exact outputs
under this module's examples directory.
"""

import argparse
from fractions import Fraction as F
import json
from pathlib import Path

from truss import solve_truss


HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def generate_case(*, case_id, panels, panel_width, height, young_modulus, area,
                  line_load, configuration):
    """One vertical plane, lower-deck nodal loads, and explicit ideal restraints.

    configuration='service': even panel count, pin/vertical roller, diagonals
    mirrored about midspan. configuration='prejoin': both translational DOFs
    held at both root chord nodes; diagonals descend toward the free end.
    Every dimensional or load input is explicit. Unsupported lengths are
    rejected by the exact truss solver.
    """
    if type(panels) is not int or panels < 1:
        raise ValueError("panels must be a positive integer")
    if configuration not in {"service", "prejoin"}:
        raise ValueError("Unknown configuration")
    if configuration == "service" and panels % 2:
        raise ValueError("Symmetric service example requires an even panel count")
    for value in (panel_width, height, young_modulus, area, line_load):
        if isinstance(value, (float, bool)) or not isinstance(value, (int, str, F)):
            raise ValueError("Use exact rational input values")
    p, h, E, A, w = map(F, (panel_width, height, young_modulus, area, line_load))
    if min(p, h, E, A) <= 0 or w < 0:
        raise ValueError("Positive geometry/material values and nonnegative line load are required")
    nodes = {}
    for i in range(panels+1):
        nodes["B"+str(i)] = [str(i*p), "0"]
        nodes["T"+str(i)] = [str(i*p), str(h)]
    members = []

    def add(name, start, end):
        members.append({"id": name, "start": start, "end": end,
                        "young_modulus": str(E), "area": str(A)})

    for i in range(panels):
        add("lower_"+str(i), "B"+str(i), "B"+str(i+1))
        add("upper_"+str(i), "T"+str(i), "T"+str(i+1))
    for i in range(panels+1):
        add("vertical_"+str(i), "B"+str(i), "T"+str(i))
    for i in range(panels):
        if configuration == "prejoin" or i < panels//2:
            add("diagonal_"+str(i), "T"+str(i), "B"+str(i+1))
        else:
            add("diagonal_"+str(i), "B"+str(i), "T"+str(i+1))
    held = [("B0", 0), ("B0", 1), ("B"+str(panels), 1)] if configuration == "service" else [
        ("B0", 0), ("B0", 1), ("T0", 0), ("T0", 1)]
    return {
        "id": case_id, "dimension": 2, "units": {"length": "m", "force": "kN"},
        "nodes": nodes, "members": members,
        "restraints": [{"node": node, "axis": axis, "displacement": "0"} for node, axis in held],
        "loads": [{"node": "B"+str(i), "components": ["0", str(-w*p/(2 if i in (0, panels) else 1))]}
                  for i in range(panels+1)],
    }


def _section_forces(*, panels, p, h, w, diagonal_length, configuration, unit_load=False):
    """Independent method-of-sections forces, tension positive.

    A separate unit downward nodal load is placed at the service midspan or
    prejoin free tip for the virtual-work displacement check. The root
    vertical in the held-root case has zero extension and zero force.
    """
    if configuration == "service":
        if unit_load:
            moment = [p*min(j, panels-j)/2 for j in range(panels+1)]
            shear = [F(1, 2) if i < panels//2 else -F(1, 2) for i in range(panels)]
        else:
            moment = [w*p*p*j*(panels-j)/2 for j in range(panels+1)]
            shear = [w*p*(F(panels-1, 2)-i) for i in range(panels)]
        result = {}
        for i in range(panels):
            result["lower_"+str(i)] = moment[i if i < panels//2 else i+1]/h
            result["upper_"+str(i)] = -moment[i+1 if i < panels//2 else i]/h
            result["diagonal_"+str(i)] = diagonal_length*abs(shear[i])/h
        for j in range(panels+1):
            result["vertical_"+str(j)] = (-shear[j] if j < panels//2 else
                                               shear[j-1] if j > panels//2 else F(0))
        return result
    if unit_load:
        moment = [p*(panels-j) for j in range(panels+1)]
        shear = [F(1)]*panels
    else:
        moment = [w*p*p*(panels-j)**2/2 for j in range(panels+1)]
        shear = [w*p*(panels-i-F(1, 2)) for i in range(panels)]
    result = {}
    for i in range(panels):
        result["lower_"+str(i)] = -moment[i]/h
        result["upper_"+str(i)] = moment[i+1]/h
        result["diagonal_"+str(i)] = diagonal_length*shear[i]/h
    for j in range(panels+1):
        result["vertical_"+str(j)] = -shear[j] if 0 < j < panels else F(0)
    return result


def verify_case(case, result, *, panels, p, h, w, configuration):
    """Exact totals, analytic sections, root moments and unit-load deflection."""
    length = p*panels
    total = -sum(F(load["components"][1]) for load in case["loads"])
    applied_moment = sum(F(case["nodes"][load["node"]][0])*F(load["components"][1])
                         for load in case["loads"])
    require(total == w*length, "Tributary nodal loads do not sum to the declared uniform load")
    require(applied_moment == -w*length*length/2,
            "Nodal load moment differs from the uniform-load analytic moment")
    require(result["status"] == "unique_linear_solution", "Illustrative truss has a mechanism")
    require(all(F(v) == 0 for v in result["total_force_equilibrium_residual"])
            and all(F(v) == 0 for v in result["reference_geometry_moment_equilibrium_residual"]),
            "Global equilibrium residual is nonzero")
    outputs = {item["id"]: item for item in result["members"]}
    diagonal = F(outputs["diagonal_0"]["length"])
    expected = _section_forces(panels=panels, p=p, h=h, w=w, diagonal_length=diagonal,
                              configuration=configuration)
    require(set(outputs) == set(expected), "Section-force coverage differs")
    require(all(F(outputs[key]["axial_force"]) == value for key, value in expected.items()),
            "Independent method-of-sections force differs from the assembled solver")
    unit = _section_forces(panels=panels, p=p, h=h, w=w, diagonal_length=diagonal,
                          configuration=configuration, unit_load=True)
    downward_deflection = sum(expected[key]*unit[key]*F(outputs[key]["length"])
                              / F(outputs[key]["axial_rigidity"]) for key in outputs)
    observed_node = "B"+str(panels//2 if configuration == "service" else panels)
    require(-F(result["nodal_displacements"][observed_node][1]) == downward_deflection,
            "Independent unit-load virtual-work deflection differs")
    reactions = result["support_reactions"]
    reaction_moment = sum(F(case["nodes"][node][0])*F(force[1])
                          - F(case["nodes"][node][1])*F(force[0])
                          for node, force in reactions.items())
    require(reaction_moment == -applied_moment, "Support reactions do not balance the external moment")
    if configuration == "service":
        require(reactions["B0"] == ["0", str(total/2)]
                and reactions["B"+str(panels)] == ["0", str(total/2)],
                "Analytic simply-supported reactions differ")
    else:
        require(F(reactions["T0"][0]) == -w*length*length/(2*h)
                and F(reactions["B0"][0]) == w*length*length/(2*h),
                "Ideal mast's horizontal reaction couple differs")
        require(F(reactions["T0"][1])+F(reactions["B0"][1]) == total,
                "Ideal mast's vertical resultant differs")
    return {
        "case_id": case["id"], "configuration": configuration,
        "span_or_reach_m": str(length), "truss_depth_m": str(h), "panels": panels,
        "line_load_per_plane_kN_m": str(w), "total_downward_load_per_plane_kN": str(total),
        "end_nodal_load_kN": str(w*p/2), "interior_nodal_load_kN": str(w*p),
        "clockwise_external_moment_about_left_root_kN_m": str(-applied_moment),
        "counterclockwise_reaction_moment_about_left_root_kN_m": str(reaction_moment),
        "observed_lower_chord_node": observed_node,
        "observed_node_displacement_m": list(result["nodal_displacements"][observed_node]),
        "downward_vertical_displacement_m": str(downward_deflection),
        "independent_member_force_checks": len(expected),
        "independent_virtual_work_displacement_check": True,
        "support_reactions_kN": {node: values for node, values in reactions.items()
                                 if any(F(value) for value in values)},
        "supported_rigidity_rank": result["supported_rigidity_rank"],
        "free_dof_count": result["free_dof_count"],
        "status": result["status"],
    }


def build_examples(write=False):
    settings = [("service-24m-plane", 8, F(6), "service"),
                ("prejoin-12m-plane", 4, F(1), "prejoin")]
    cases, results, summaries = {}, {}, []
    for name, panels, w, configuration in settings:
        case = generate_case(case_id=name, panels=panels, panel_width=3, height=4,
                             young_modulus=200000000, area="3/1000", line_load=w,
                             configuration=configuration)
        result = solve_truss(case)
        summary = verify_case(case, result, panels=panels, p=F(3), h=F(4), w=w,
                              configuration=configuration)
        cases[name], results[name] = case, result
        summaries.append(summary)
    require(F(5+1)*F(2)/2 == 6 and F(1)*F(2)/2 == 1, "Declared pressure-to-plane conversion failed")
    report = {
        "checks_passed": True, "study_type": "synthetic_declared_geometry_and_loading",
        "site_derived_span": False, "physical_approval_or_load_rating": False,
        "declared_inputs": {
            "panel_width_m": "3", "truss_depth_m": "4", "deck_width_m": "2", "plane_count": 2,
            "crowd_pressure_kPa_service_only": "5", "dead_pressure_kPa": "1",
            "young_modulus_kN_m2": "200000000", "member_area_m2": "3/1000",
            "material_and_area_basis": "Arbitrary illustrative elastic inputs, not selected or capacity-verified members",
            "dead_load_basis": "Prescribed uniform input, not calculated member/equipment mass",
        },
        "idealizations": [
            "Vertical load is shared equally by two independent truss planes; no spatial deck-load distribution is solved",
            "Service top and lower midspan nodes are each shared: ideal continuity at both chord levels, not a single center hinge",
            "The shared nodes do not verify a physical connection's strength, stiffness, alignment, fatigue or reliability",
            "Prejoin root upper and lower nodes have both translations held by an ideal mast; no mast or chassis capacity is calculated",
            "Prejoin contains declared dead load only and no public crowd load",
            "Only two idealized endpoint configurations are analyzed; intermediate extension, lifting and joining states are absent",
            "No buckling, lateral stability, torsion, second-order response, vibration, connection or mast capacity is established",
            "The 24 m span and 12 m reach are illustrative and are not derived from any Tillamook road or site measurement",
        ],
        "cases": summaries,
        "independent_checks": {"member_section_force_checks": sum(s["independent_member_force_checks"] for s in summaries),
                               "unit_load_virtual_work_deflections": 2, "load_and_reaction_resultants": 2},
    }
    if write:
        destination = HERE/"examples"
        destination.mkdir(exist_ok=True)
        for name in cases:
            (destination/(name+".json")).write_text(json.dumps(cases[name], indent=2)+"\n")
            (destination/(name+"-result.json")).write_text(json.dumps(results[name], indent=2)+"\n")
        (destination/"bridge-example-report.json").write_text(json.dumps(report, indent=2)+"\n")
    return report


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true", help="save generated examples and exact outputs under v2/examples")
    args = parser.parse_args()
    print(json.dumps(build_examples(write=args.write), indent=2))
