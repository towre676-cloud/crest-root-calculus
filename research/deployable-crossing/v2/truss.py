"""Small exact 2D/3D linear elastic pin-jointed truss kernel.

This is a research analysis model, not a bridge design, load rating, code
compliance check, or permission to build or operate equipment. All geometry,
material properties, cross-sectional areas, restraints and nodal loads must
be declared. Loads such as self-weight are NOT inferred. The model has axial
bars and translational DOFs; it has no beam bending, joint stiffness, member
buckling, geometric nonlinearity, capacity limits, dynamics or foundations.

solve_truss(case) accepts a JSON-compatible mapping with exactly these keys:
  id: nonempty string
  dimension: 2 or 3
  units: {"length": nonempty label, "force": nonempty label}
  nodes: {node_id: [coordinates in the declared length unit]}
  members: [{"id", "start", "end", "young_modulus", "area"}]
  restraints: [{"node", "axis", "displacement"}]
  loads: [{"node", "components": [force components]}]

Axis indices are zero-based. E has force/length^2 units and area length^2.
Unit labels are recorded, not converted or independently validated. Multiple
loads at one node are summed. A restraint specifies an exact displacement,
which can be nonzero. Every scalar is int, Fraction, or a rational/decimal
string. Floats and booleans are rejected. Member lengths must be rational:
unsupported irrational lengths raise UnsupportedGeometryError, without any
rounding. Optional symbolic verification lives in verify_truss.py.

Results contain exact rational strings, assembly data, free-equilibrium
residuals, member forces/stresses/extensions, and reactions. Singular models
return an exact supported infinitesimal-mechanism basis and its virtual work
against the declared loading; they do not return a spurious unique solution.
"""

import argparse
from fractions import Fraction
from math import isqrt
import json
from pathlib import Path


class ModelError(ValueError):
    """Malformed or unsupported declared model."""


class UnsupportedGeometryError(ModelError):
    """A member length is outside this rational arithmetic kernel."""


NOT_CHECKED = [
    "Material strength, yield, fatigue or permissible stress",
    "Member or global buckling and geometric nonlinear stability",
    "Connection, lock, bearing, weld or fastener behavior",
    "Beam bending, torsion, joint stiffness or deck plate behavior",
    "Self-weight or other omitted loads, load combinations and moving loads",
    "Dynamics, vibration, impact, wind response or deployment transients",
    "Ground contact, foundations, soil, anchorage or support capacity",
    "Clearance, accessibility, traffic control, code compliance or construction approval",
]


def exact(value):
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise ModelError("Scalars must be exact integers or rational/decimal strings; floats are unsupported")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ModelError("Invalid rational scalar") from error


def _keys(value, required, name):
    if not isinstance(value, dict) or set(value) != set(required):
        raise ModelError(name + " must contain exactly: " + ", ".join(required))


def _label(value, name):
    if not isinstance(value, str) or not value.strip():
        raise ModelError(name + " must be a nonempty string")
    return value


def _array(value, name, length=None):
    if not isinstance(value, (list, tuple)) or (length is not None and len(value) != length):
        raise ModelError(name + " has an invalid array shape")
    return value


def _dot(a, b):
    return sum((x*y for x, y in zip(a, b)), Fraction(0))


def _matrix_vector(matrix, vector):
    return [_dot(row, vector) for row in matrix]


def _rref(matrix, ncolumns):
    """Exact row reduction; augmented columns, if present, are never pivots."""
    rows = [list(row) for row in matrix]
    pivots = []
    position = 0
    for column in range(ncolumns):
        found = next((i for i in range(position, len(rows)) if rows[i][column]), None)
        if found is None:
            continue
        rows[position], rows[found] = rows[found], rows[position]
        divisor = rows[position][column]
        rows[position] = [v/divisor for v in rows[position]]
        for i in range(len(rows)):
            if i != position and rows[i][column]:
                multiplier = rows[i][column]
                rows[i] = [a-multiplier*b for a, b in zip(rows[i], rows[position])]
        pivots.append(column)
        position += 1
        if position == len(rows):
            break
    return rows, pivots


def _nullspace(matrix, ncolumns):
    reduced, pivots = _rref(matrix, ncolumns)
    basis = []
    for column in range(ncolumns):
        if column not in pivots:
            vector = [Fraction(0)] * ncolumns
            vector[column] = Fraction(1)
            for row, pivot in enumerate(pivots):
                vector[pivot] = -reduced[row][column]
            basis.append(vector)
    return len(pivots), basis


def _solve_nonsingular(matrix, rhs):
    size = len(rhs)
    if not size:
        return []
    reduced, pivots = _rref([row+[value] for row, value in zip(matrix, rhs)], size)
    if len(pivots) != size:
        raise ArithmeticError("Positive-definite reduced stiffness unexpectedly singular")
    answer = [Fraction(0)] * size
    for row, column in enumerate(pivots):
        answer[column] = reduced[row][-1]
    if _matrix_vector(matrix, answer) != rhs:
        raise ArithmeticError("Exact free-equilibrium solve failed")
    return answer


def _rational_length(squared, member_id):
    if squared <= 0:
        raise ModelError("Zero-length member: " + member_id)
    numerator, denominator = isqrt(squared.numerator), isqrt(squared.denominator)
    if numerator*numerator != squared.numerator or denominator*denominator != squared.denominator:
        raise UnsupportedGeometryError(
            "Member " + member_id + " has irrational length sqrt(" + str(squared)
            + "); rational-length geometry is required and no rounding was performed")
    return Fraction(numerator, denominator)


def _parse(case):
    _keys(case, ["id", "dimension", "units", "nodes", "members", "restraints", "loads"], "case")
    _label(case["id"], "id")
    dimension = case["dimension"]
    if type(dimension) is not int or dimension not in (2, 3):
        raise ModelError("dimension must be the integer 2 or 3")
    _keys(case["units"], ["length", "force"], "units")
    for label in ("length", "force"):
        _label(case["units"][label], "units."+label)
    if not isinstance(case["nodes"], dict) or len(case["nodes"]) < 2:
        raise ModelError("At least two distinct node identifiers are required")
    nodes = {}
    for name, coordinates in case["nodes"].items():
        _label(name, "node identifier")
        nodes[name] = tuple(exact(v) for v in _array(coordinates, "node coordinates", dimension))
    order = list(nodes)
    indices = {name: i for i, name in enumerate(order)}
    ndof = dimension * len(nodes)
    member_input = _array(case["members"], "members")
    if not member_input:
        raise ModelError("At least one member is required")
    members, names = [], set()
    for member in member_input:
        _keys(member, ["id", "start", "end", "young_modulus", "area"], "member")
        name = _label(member["id"], "member.id")
        if name in names:
            raise ModelError("Duplicate member identifier: " + name)
        names.add(name)
        start, end = member["start"], member["end"]
        if not isinstance(start, str) or not isinstance(end, str) or start not in nodes or end not in nodes:
            raise ModelError("Member endpoint does not identify an existing node")
        if start == end:
            raise ModelError("Member endpoints must differ")
        E, area = exact(member["young_modulus"]), exact(member["area"])
        if E <= 0 or area <= 0:
            raise ModelError("Member modulus and area must be strictly positive")
        delta = tuple(b-a for a, b in zip(nodes[start], nodes[end]))
        length = _rational_length(_dot(delta, delta), name)
        row = [Fraction(0)] * ndof
        for axis, value in enumerate(delta):
            row[dimension*indices[start]+axis] = -value
            row[dimension*indices[end]+axis] = value
        members.append({"id": name, "start": start, "end": end, "delta": delta,
                        "length": length, "E": E, "area": area, "EA": E*area, "row": row})
    prescribed = {}
    for restraint in _array(case["restraints"], "restraints"):
        _keys(restraint, ["node", "axis", "displacement"], "restraint")
        name, axis = restraint["node"], restraint["axis"]
        if not isinstance(name, str) or name not in nodes or type(axis) is not int or not 0 <= axis < dimension:
            raise ModelError("Restraint node or coordinate axis is invalid")
        dof = dimension*indices[name] + axis
        if dof in prescribed:
            raise ModelError("A displacement coordinate was restrained more than once")
        prescribed[dof] = exact(restraint["displacement"])
    force = [Fraction(0)] * ndof
    for load in _array(case["loads"], "loads"):
        _keys(load, ["node", "components"], "load")
        name = load["node"]
        if not isinstance(name, str) or name not in nodes:
            raise ModelError("Load node does not exist")
        for axis, value in enumerate(_array(load["components"], "load.components", dimension)):
            force[dimension*indices[name]+axis] += exact(value)
    return dimension, nodes, order, members, prescribed, force


def _by_node(values, order, dimension):
    return {name: [str(v) for v in values[i*dimension:(i+1)*dimension]]
            for i, name in enumerate(order)}


def _external_resultants(nodes, order, dimension, forces):
    force_sum = [Fraction(0)] * dimension
    moment_sum = [Fraction(0)] * (1 if dimension == 2 else 3)
    for index, name in enumerate(order):
        position = nodes[name]
        force = forces[index*dimension:(index+1)*dimension]
        for axis, value in enumerate(force):
            force_sum[axis] += value
        if dimension == 2:
            moment_sum[0] += position[0]*force[1] - position[1]*force[0]
        else:
            for axis, j, k in ((0, 1, 2), (1, 2, 0), (2, 0, 1)):
                moment_sum[axis] += position[j]*force[k] - position[k]*force[j]
    return force_sum, moment_sum


def solve_truss(case):
    """Analyze one declared rational-length, small-displacement axial truss."""
    dimension, nodes, order, members, prescribed, force = _parse(case)
    ndof = len(force)
    free = [i for i in range(ndof) if i not in prescribed]
    K = [[Fraction(0)] * ndof for _ in range(ndof)]
    for member in members:
        weight = member["EA"] / member["length"]**3
        row = member["row"]
        nonzero = [i for i, value in enumerate(row) if value]
        for i in nonzero:
            for j in nonzero:
                K[i][j] += weight*row[i]*row[j]
    reduced_rigidity = [[member["row"][i] for i in free] for member in members]
    rank, mechanisms = _nullspace(reduced_rigidity, len(free))
    Kff = [[K[i][j] for j in free] for i in free]
    rhs = [force[i] - sum((K[i][j]*value for j, value in prescribed.items()), Fraction(0))
           for i in free]
    common = {
        "case_id": case["id"], "model": "linear_elastic_pin_jointed_axial_truss",
        "dimension": dimension, "units": dict(case["units"]), "node_order": order,
        "degree_order": [{"node": name, "axis": axis} for name in order for axis in range(dimension)],
        "free_dof_indices": free, "prescribed_dof_indices": sorted(prescribed),
        "free_dof_count": len(free), "supported_rigidity_rank": rank,
        "member_count": len(members),
        "supported_self_stress_dimension": len(members)-rank,
        "rigidity_matrix": [[str(v) for v in member["row"]] for member in members],
        "stiffness_matrix": [[str(v) for v in row] for row in K],
        "applied_nodal_forces": _by_node(force, order, dimension),
        "checks_not_performed": list(NOT_CHECKED),
        "scope": "Exact solution of declared linear axial equations only; no physical approval",
    }
    if mechanisms:
        witnesses = []
        compatible = True
        for mechanism in mechanisms:
            full = [Fraction(0)] * ndof
            for i, value in zip(free, mechanism):
                full[i] = value
            virtual_work = _dot(mechanism, rhs)
            compatible = compatible and virtual_work == 0
            if any(_dot(member["row"], full) for member in members):
                raise ArithmeticError("Mechanism witness failed the bar-extension equations")
            witnesses.append({"nodal_motion": _by_node(full, order, dimension),
                              "load_virtual_work": str(virtual_work)})
        common.update(status="mechanism", mechanism_dimension=len(mechanisms),
                      equilibrium_compatible=compatible, mechanisms=witnesses,
                      interpretation=("Declared load has equilibria but no unique displacement solution"
                                      if compatible else "Declared load is incompatible with linear static equilibrium"))
        return common

    displacement = [Fraction(0)] * ndof
    for i, value in prescribed.items():
        displacement[i] = value
    for i, value in zip(free, _solve_nonsingular(Kff, rhs)):
        displacement[i] = value
    internal = _matrix_vector(K, displacement)
    residual = [a-b for a, b in zip(internal, force)]
    if any(residual[i] for i in free):
        raise ArithmeticError("Free displacement solution violates equilibrium")
    reactions = [residual[i] if i in prescribed else Fraction(0) for i in range(ndof)]
    member_results = []
    energy_sum = Fraction(0)
    for member in members:
        projection = _dot(member["row"], displacement)
        extension = projection/member["length"]
        strain = extension/member["length"]
        axial_force = member["EA"]*strain
        energy = member["EA"]*extension**2/(2*member["length"])
        energy_sum += energy
        member_results.append({
            "id": member["id"], "start": member["start"], "end": member["end"],
            "length": str(member["length"]), "axial_rigidity": str(member["EA"]),
            "extension": str(extension), "axial_strain": str(strain),
            "axial_force": str(axial_force), "axial_stress": str(axial_force/member["area"]),
            "strain_energy": str(energy),
            "force_sign": "tension" if axial_force > 0 else "compression" if axial_force < 0 else "zero",
        })
    energy_global = _dot(displacement, internal)/2
    if energy_sum != energy_global:
        raise ArithmeticError("Member energies do not sum to assembled elastic energy")
    resultant_force, resultant_moment = _external_resultants(
        nodes, order, dimension, [a+b for a, b in zip(force, reactions)])
    if any(resultant_force) or any(resultant_moment):
        raise ArithmeticError("External forces and reactions fail total reference-geometry equilibrium")
    common.update(
        status="unique_linear_solution", mechanism_dimension=0,
        nodal_displacements=_by_node(displacement, order, dimension),
        support_reactions=_by_node(reactions, order, dimension),
        free_equilibrium_residuals=[str(residual[i]) for i in free],
        members=member_results, strain_energy=str(energy_sum),
        applied_force_work=str(_dot(force, displacement)),
        support_reaction_work=str(_dot(reactions, displacement)),
        total_force_equilibrium_residual=[str(v) for v in resultant_force],
        reference_geometry_moment_equilibrium_residual=[str(v) for v in resultant_moment],
        interpretation="Unique displacement solution within the declared linear axial model",
    )
    return common


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("case", type=Path, help="one JSON truss model using the documented schema")
    args = parser.parse_args()
    print(json.dumps(solve_truss(json.loads(args.case.read_text(encoding="utf-8"))), indent=2))
