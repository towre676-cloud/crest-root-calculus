"""Independent small exact controls for truss.py; standard library by default.

Run this file directly. --symbolic adds independent SymPy matrix solves and
rank checks. Synthetic algebra fixtures are not a proposed bridge or assigned
material specification. The default verifier uses hand-derived expectations,
normalized-direction assembly, determinant expansion/Cramer's rule, member
compatibility, nodal equilibrium, force/moment balance and energy identities.
"""

import argparse
from copy import deepcopy
from fractions import Fraction as F
import importlib
from itertools import permutations
from math import isqrt
import json

from truss import ModelError, UnsupportedGeometryError, solve_truss


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def member(name, start, end, E, area):
    return {"id": name, "start": start, "end": end, "young_modulus": E, "area": area}


def restraint(node, axis, displacement=0):
    return {"node": node, "axis": axis, "displacement": displacement}


def example_cases():
    """Synthetic declared inputs, ready for JSON serialization."""
    units = {"length": "synthetic_length_unit", "force": "synthetic_force_unit"}
    triangle = {
        "id": "synthetic_345_triangle", "dimension": 2, "units": dict(units),
        "nodes": {"A": [0, 0], "B": [8, 0], "C": [4, 3]},
        "members": [member("AB", "A", "B", 1000, 2), member("AC", "A", "C", 1000, 2),
                    member("BC", "B", "C", 1000, 2)],
        "restraints": [restraint("A", 0), restraint("A", 1), restraint("B", 1)],
        "loads": [{"node": "C", "components": [0, -6]}],
    }
    rectangle = {
        "id": "synthetic_redundant_34_rectangle", "dimension": 2, "units": dict(units),
        "nodes": {"A": [0, 0], "B": [4, 0], "C": [4, 3], "D": [0, 3]},
        "members": [member(name, start, end, 60000, 2) for name, start, end in
                    [("AB", "A", "B"), ("BC", "B", "C"), ("CD", "C", "D"),
                     ("DA", "D", "A"), ("AC", "A", "C"), ("BD", "B", "D")]],
        "restraints": [restraint("A", 0), restraint("A", 1), restraint("B", 1)],
        "loads": [{"node": "C", "components": [0, -24]},
                  {"node": "D", "components": [0, -24]}],
    }
    spatial = {
        "id": "synthetic_spatial_three_bar_control", "dimension": 3, "units": dict(units),
        "nodes": {"O": [0, 0, 0], "A": [3, 4, 0], "B": [3, 0, 4], "C": [0, 3, 4]},
        "members": [member("O"+name, "O", name, 25000, 5) for name in ("A", "B", "C")],
        "restraints": [restraint(name, axis) for name in ("A", "B", "C") for axis in range(3)],
        "loads": [{"node": "O", "components": ["18/5", "47/20", 7]}],
    }
    settlement = {
        "id": "synthetic_prescribed_support_movement", "dimension": 2, "units": dict(units),
        "nodes": {"A": [0, 0], "B": [3, 0], "C": [6, 0]},
        "members": [member("AB", "A", "B", 200, 2), member("BC", "B", "C", 200, 2)],
        "restraints": [restraint("A", 0), restraint("A", 1), restraint("B", 1),
                       restraint("C", 0, "1/100"), restraint("C", 1)],
        "loads": [{"node": "A", "components": [1, 0]}],
    }
    return [triangle, rectangle, spatial, settlement]


def determinant(matrix):
    """Independent permutation expansion, suitable only for small controls."""
    n = len(matrix)
    require(n <= 6, "Independent verification fixture exceeds the small determinant limit")
    answer = F(0)
    for order in permutations(range(n)):
        inversions = sum(order[i] > order[j] for i in range(n) for j in range(i+1, n))
        term = F(-1 if inversions % 2 else 1)
        for i, j in enumerate(order):
            term *= matrix[i][j]
        answer += term
    return answer


def independent_system(case):
    """Assemble with unit directions and EA/L, unlike the kernel's raw R."""
    d = case["dimension"]
    names = list(case["nodes"])
    offset = {name: d*i for i, name in enumerate(names)}
    size = d*len(names)
    K = [[F(0)]*size for _ in range(size)]
    members = []
    for edge in case["members"]:
        start, end = edge["start"], edge["end"]
        vector = [F(b)-F(a) for a, b in zip(case["nodes"][start], case["nodes"][end])]
        squared = sum(x*x for x in vector)
        length = F(isqrt(squared.numerator), isqrt(squared.denominator))
        require(length*length == squared, "Verifier fixture has an irrational length")
        direction = [x/length for x in vector]
        E, area = F(edge["young_modulus"]), F(edge["area"])
        scalar = E*area/length
        for a in range(d):
            for b in range(d):
                component = scalar*direction[a]*direction[b]
                for node_i, sign_i in ((start, -1), (end, 1)):
                    for node_j, sign_j in ((start, -1), (end, 1)):
                        K[offset[node_i]+a][offset[node_j]+b] += sign_i*sign_j*component
        members.append((edge, length, direction, E, area))
    fixed = {offset[item["node"]]+item["axis"]: F(item["displacement"])
             for item in case["restraints"]}
    free = [i for i in range(size) if i not in fixed]
    force = [F(0)]*size
    for load in case["loads"]:
        for axis, value in enumerate(load["components"]):
            force[offset[load["node"]]+axis] += F(value)
    reduced = [[K[i][j] for j in free] for i in free]
    rhs = [force[i]-sum(K[i][j]*value for j, value in fixed.items()) for i in free]
    return d, names, offset, K, members, fixed, free, force, reduced, rhs


def validate_solution(case, receipt, symbolic=None):
    d, names, offset, K, members, fixed, free, force, reduced, rhs = independent_system(case)
    require(receipt["status"] == "unique_linear_solution", "Unique solution status missing")
    require([[F(x) for x in row] for row in receipt["stiffness_matrix"]] == K,
            "Independent unit-direction stiffness differs")
    denominator = determinant(reduced)
    require(denominator > 0, "Fixture lacks positive material stiffness")
    u = [F(v) for name in names for v in receipt["nodal_displacements"][name]]
    reactions = [F(v) for name in names for v in receipt["support_reactions"][name]]
    for i, value in fixed.items():
        require(u[i] == value, "Prescribed displacement changed")
    for column, i in enumerate(free):
        replaced = [list(row) for row in reduced]
        for row in range(len(free)):
            replaced[row][column] = rhs[row]
        require(u[i] == determinant(replaced)/denominator, "Independent Cramer solution differs")
        require(reactions[i] == 0, "Reaction reported at a free coordinate")
    if symbolic is not None:
        sp = symbolic
        answer = sp.Matrix(reduced).inv()*sp.Matrix(rhs)
        require([F(x) for x in answer] == [u[i] for i in free], "Symbolic solve differs")
        require(sp.Matrix(reduced).rank() == len(free), "Symbolic stiffness rank differs")
    internal = [F(0)]*len(u)
    energies = []
    returned = {item["id"]: item for item in receipt["members"]}
    require(set(returned) == {item[0]["id"] for item in members}, "Member output coverage changed")
    for edge, length, direction, E, area in members:
        start, end = offset[edge["start"]], offset[edge["end"]]
        extension = sum(direction[k]*(u[end+k]-u[start+k]) for k in range(d))
        axial = E*area*extension/length
        energy = axial*extension/2
        row = returned[edge["id"]]
        require(F(row["length"]) == length and F(row["extension"]) == extension,
                "Member geometric compatibility failed")
        require(F(row["axial_force"]) == axial and F(row["axial_stress"]) == axial/area,
                "Member force or stress violates the constitutive equation")
        require(F(row["axial_strain"]) == extension/length and F(row["strain_energy"]) == energy,
                "Member strain or energy changed")
        require(row["force_sign"] == ("tension" if axial > 0 else "compression" if axial < 0 else "zero"),
                "Tension/compression sign changed")
        for axis in range(d):
            internal[start+axis] -= axial*direction[axis]
            internal[end+axis] += axial*direction[axis]
        energies.append(energy)
    require(internal == [f+r for f, r in zip(force, reactions)], "Independent joint equilibrium failed")
    require(F(receipt["strain_energy"]) == sum(energies), "Total member energy changed")
    require(2*sum(energies) == sum((f+r)*v for f, r, v in zip(force, reactions, u)),
            "Load/displacement and member energies disagree")
    require(F(receipt["applied_force_work"]) == sum(f*v for f, v in zip(force, u))
            and F(receipt["support_reaction_work"]) == sum(r*v for r, v in zip(reactions, u)),
            "Reported work changed")
    require(all(F(v) == 0 for v in receipt["free_equilibrium_residuals"]), "Nonzero free residual")
    require(all(sum(force[offset[name]+axis]+reactions[offset[name]+axis] for name in names) == 0
                for axis in range(d)), "Total external force does not balance")
    moments = []
    pairs = [(0, 1)] if d == 2 else [(1, 2), (2, 0), (0, 1)]
    for a, b in pairs:
        value = sum(F(case["nodes"][name][a])*(force[offset[name]+b]+reactions[offset[name]+b])
                    - F(case["nodes"][name][b])*(force[offset[name]+a]+reactions[offset[name]+a])
                    for name in names)
        moments.append(value)
    require(all(v == 0 for v in moments), "Total reference-geometry moment does not balance")
    require([F(x) for x in receipt["reference_geometry_moment_equilibrium_residual"]] == moments
            and [F(x) for x in receipt["total_force_equilibrium_residual"]] == [F(0)]*d,
            "Reported global equilibrium residual changed")


def validate_mechanism(case, receipt, nullity, compatible, symbolic=None):
    d, names, offset, K, members, fixed, free, force, reduced, rhs = independent_system(case)
    require(receipt["status"] == "mechanism" and receipt["mechanism_dimension"] == nullity,
            "Mechanism classification changed")
    require(receipt["equilibrium_compatible"] is compatible, "Load compatibility classification changed")
    require("nodal_displacements" not in receipt, "Singular model presented an arbitrary displacement as unique")
    require(receipt["supported_rigidity_rank"] == len(free)-nullity, "Fixture rank changed")
    require(determinant(reduced) == 0, "Independent determinant does not show singularity")
    vectors, works = [], []
    for witness in receipt["mechanisms"]:
        vector = [F(v) for name in names for v in witness["nodal_motion"][name]]
        require(all(vector[i] == 0 for i in fixed), "Mechanism moves a restrained coordinate")
        for edge, length, direction, E, area in members:
            start, end = offset[edge["start"]], offset[edge["end"]]
            require(sum(direction[k]*(vector[end+k]-vector[start+k]) for k in range(d)) == 0,
                    "Mechanism stretches a bar to first order")
        work = sum(a*b for a, b in zip(vector, force))
        require(F(witness["load_virtual_work"]) == work, "Mechanism virtual-work certificate changed")
        vectors.append(vector)
        works.append(work)
    require(len(vectors) == nullity, "Incomplete mechanism basis")
    gram = [[sum(a*b for a, b in zip(v, w)) for w in vectors] for v in vectors]
    require(determinant(gram) > 0, "Mechanism witnesses are linearly dependent")
    require(all(work == 0 for work in works) is compatible, "Witnesses do not support load compatibility")
    if symbolic is not None:
        require(symbolic.Matrix(reduced).rank() == len(free)-nullity,
                "Independent symbolic mechanism rank differs")


def verify_all(symbolic=False):
    sp = importlib.import_module("sympy") if symbolic else None
    cases = example_cases()
    receipts = [solve_truss(case) for case in cases]
    for case, receipt in zip(cases, receipts):
        validate_solution(case, receipt, sp)
    triangle, rectangle, spatial, settlement = receipts
    require(triangle["nodal_displacements"] == {"A": ["0", "0"], "B": ["2/125", "0"],
                                                "C": ["1/125", "-63/2000"]}, "Analytic triangle deflection differs")
    require([F(m["axial_force"]) for m in triangle["members"]] == [4, -5, -5],
            "Analytic triangle forces differ")
    require(triangle["support_reactions"] == {"A": ["0", "3"], "B": ["0", "3"], "C": ["0", "0"]},
            "Analytic triangle reactions differ")
    require(F(triangle["strain_energy"]) == F(189, 2000), "Analytic triangle energy differs")
    require(rectangle["nodal_displacements"] == {"A": ["0", "0"], "B": ["1/7500", "0"],
                "C": ["1/7500", "-21/40000"], "D": ["0", "-21/40000"]},
            "Redundant rectangle compatibility differs")
    require([F(m["axial_force"]) for m in rectangle["members"]] == [4, -21, 4, -21, -5, -5]
            and rectangle["supported_self_stress_dimension"] == 1,
            "Redundant rectangle force sharing differs")
    require(F(rectangle["strain_energy"]) == F(63, 5000), "Redundant rectangle energy differs")
    require(spatial["nodal_displacements"]["O"] == ["1/10000", "-1/20000", "1/5000"],
            "Non-axis-aligned spatial displacement differs")
    require([F(m["axial_force"]) for m in spatial["members"]] == [-F(1, 2), -F(11, 2), -F(13, 4)]
            and F(spatial["strain_energy"]) == F(657, 800000), "Spatial force/energy differs")
    require(settlement["nodal_displacements"]["B"] == ["1/200", "0"]
            and settlement["support_reactions"]["A"] == ["-5/3", "0"]
            and settlement["support_reactions"]["C"] == ["2/3", "0"]
            and [F(m["axial_force"]) for m in settlement["members"]] == [F(2, 3)]*2
            and F(settlement["strain_energy"]) == F(1, 300),
            "Prescribed support movement or load-at-support reaction differs")

    mechanism_cases = []
    missing = deepcopy(cases[0]); missing["members"].pop()
    mechanism_cases.append((missing, 1, False))
    unloaded = deepcopy(missing); unloaded["loads"] = []
    mechanism_cases.append((unloaded, 1, True))
    free_body = deepcopy(cases[0]); free_body["restraints"] = []; free_body["loads"] = []
    mechanism_cases.append((free_body, 3, True))
    mechanism_receipts = []
    for case, nullity, compatible in mechanism_cases:
        receipt = solve_truss(case)
        validate_mechanism(case, receipt, nullity, compatible, sp)
        mechanism_receipts.append(receipt)

    # Response must be independent of how a member is oriented in the input.
    reversed_case = deepcopy(cases[1])
    for edge in reversed_case["members"]:
        edge["start"], edge["end"] = edge["end"], edge["start"]
    reversed_receipt = solve_truss(reversed_case)
    validate_solution(reversed_case, reversed_receipt, sp)
    require(reversed_receipt["nodal_displacements"] == rectangle["nodal_displacements"]
            and [m["axial_force"] for m in reversed_receipt["members"]]
                == [m["axial_force"] for m in rectangle["members"]], "Member orientation changed physical response")
    split = deepcopy(cases[0]); split["loads"] = [{"node": "C", "components": [0, -2]},
                                                 {"node": "C", "components": [0, -4]}]
    require(solve_truss(split)["nodal_displacements"] == triangle["nodal_displacements"],
            "Equivalent nodal-load decomposition changed response")

    invalid = []
    for edit in (
        lambda c: c["nodes"]["C"].__setitem__(0, 4.0),
        lambda c: c["nodes"].__setitem__("C", [0, 0]),
        lambda c: c["members"][0].__setitem__("young_modulus", -1),
        lambda c: c["members"][0].__setitem__("area", 0),
        lambda c: c["members"][0].__setitem__("end", "missing"),
        lambda c: c["restraints"].append(deepcopy(c["restraints"][0])),
        lambda c: c["loads"][0].__setitem__("components", [0]),
        lambda c: c.__setitem__("dimension", True),
        lambda c: c["units"].__setitem__("length", ""),
        lambda c: c["members"][0].__setitem__("area", True),
        lambda c: c["members"][0].__setitem__("young_modulus", "1/0"),
        lambda c: c["restraints"][0].__setitem__("axis", 2),
    ):
        case = deepcopy(cases[0]); edit(case); invalid.append(case)
    for case in invalid:
        try: solve_truss(case)
        except ModelError: pass
        else: raise ArithmeticError("Malformed or invalid physical-domain input accepted")
    irrational = deepcopy(cases[0]); irrational["nodes"]["C"] = [4, 4]
    try: solve_truss(irrational)
    except UnsupportedGeometryError: pass
    else: raise ArithmeticError("Irrational member length was silently accepted")

    tampered = []
    for edit in (
        lambda r: r["nodal_displacements"]["C"].__setitem__(0, "0"),
        lambda r: r["support_reactions"]["A"].__setitem__(1, "4"),
        lambda r: r["members"][0].__setitem__("axial_force", "-4"),
        lambda r: r.__setitem__("strain_energy", "1"),
    ):
        changed = deepcopy(triangle); edit(changed); tampered.append(changed)
    for changed in tampered:
        try: validate_solution(cases[0], changed, sp)
        except ArithmeticError: pass
        else: raise ArithmeticError("Tampered solution certificate accepted")
    changed = deepcopy(mechanism_receipts[0]); changed["mechanisms"][0]["nodal_motion"]["A"][0] = "1"
    try: validate_mechanism(mechanism_cases[0][0], changed, 1, False, sp)
    except ArithmeticError: pass
    else: raise ArithmeticError("Tampered mechanism witness accepted")

    return {
        "checks_passed": True,
        "verification_mode": "standard_library_plus_symbolic" if symbolic else "standard_library",
        "controls": {"analytic_unique_solutions": 4, "independent_cramer_solves": 5,
                     "mechanism_and_load_compatibility_cases": 3, "orientation_and_load_decomposition": 2,
                     "malformed_or_invalid_inputs_rejected": len(invalid),
                     "unsupported_irrational_lengths_rejected": 1, "tampered_certificates_rejected": 5},
        "optional_symbolic_checks": {"enabled": symbolic,
                    "unique_solves": 5 if symbolic else 0, "mechanism_ranks": 3 if symbolic else 0},
        "example_summaries": [{key: r[key] for key in ("case_id", "status", "dimension", "free_dof_count",
                            "supported_rigidity_rank", "supported_self_stress_dimension", "strain_energy")}
                              for r in receipts],
        "scope": "Synthetic exact axial-model verification only; no site, capacity or physical approval",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbolic", action="store_true", help="add optional independent SymPy checks")
    args = parser.parse_args()
    print(json.dumps(verify_all(symbolic=args.symbolic), indent=2))
