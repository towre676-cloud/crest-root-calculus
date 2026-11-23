"""Independent exact controls for the illustrative crossing research screens.

Run: python research/deployable-crossing/verify.py
The default uses only the standard library, hand-derived fixture ranks and
corner reactions, and exact constraint/equilibrium residuals. Add --symbolic
for extra SymPy rank, equilibrium and integral cross-checks. No physical
measurements, real machine specification, or bridge approval is provided by
these tests. All numerical fixtures are synthetic algebra cases.
"""

from copy import deepcopy
from fractions import Fraction as F
from itertools import product
import argparse
import importlib
import json

from model import (
    DEPLOYMENT_RECORD_ORDER, Interval, cantilever_uniform_load,
    deployment_sequence_screen, exact, horizontal_load_moment,
    join_geometry_screen, ramp_slope_screen, rigid_half_kinematics,
    support_reaction_screen,
)


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def bounds(receipt):
    return F(receipt["lower"]), F(receipt["upper"])


def corners(value):
    if isinstance(value, Interval):
        return (value.lower, value.upper)
    return (F(value), F(value))


def validate_support_receipt(receipt, *, load, moment, left, right,
                             known_reactions, symbolic_module=None):
    """Verify hand-computed corner reactions by both equilibrium equations.

    Reactions are affine in (W,M), so these corners attain their extrema.
    This is a justified exact-box test, not a sampling claim for nonlinear
    load paths or for parameters absent from the declared model.
    """
    input_corners = list(product(corners(load), corners(moment)))
    require(len(known_reactions) == len(input_corners), "Incomplete independent fixture")
    answers = [tuple(F(r) for r in pair) for pair in known_reactions]
    for (w, m), (rl, rr) in zip(input_corners, answers):
        require(rl + rr == w and F(left)*rl + F(right)*rr == m,
                "Hand-derived reaction fixture fails equilibrium")
    if symbolic_module is not None:
        sp = symbolic_module
        coefficients = sp.Matrix([[1, 1], [sp.Rational(left), sp.Rational(right)]])
        for (w, m), expected in zip(input_corners, answers):
            answer = coefficients.inv() * sp.Matrix([sp.Rational(w), sp.Rational(m)])
            require(tuple(F(r) for r in answer) == expected,
                    "Independent symbolic equilibrium solve differs")
    extrema = [(min(a[i] for a in answers), max(a[i] for a in answers))
               for i in range(2)]
    states = []
    for field, check, expected in zip(
            ("reaction_left", "reaction_right"),
            ("left_contact", "right_contact"), extrema):
        require(bounds(receipt[field]) == expected, "Reaction enclosure changed")
        require(bounds(receipt["checks"][check]["margin"]) == expected,
                "Contact margin differs from its reaction")
        lo, hi = expected
        state = ("inside" if lo > 0 else "outside" if hi < 0
                 else "marginal" if lo == hi == 0 else "indeterminate")
        require(receipt["checks"][check]["status"] == state,
                "Contact status differs from exact corner signs")
        states.append(state)
    expected = next((s for s in ("outside", "indeterminate", "marginal")
                     if s in states), "inside")
    require(receipt["status"] == expected, "Overall support status changed")


def validate_kinematic_receipt(receipt, expected_rank, symbolic_module=None):
    matrix = [[F(v) for v in row] for row in receipt["constraint_matrix"]]
    require(expected_rank == receipt["constraint_rank"],
            "Hand-derived rigid-body rank mismatch")
    if symbolic_module is not None:
        require(symbolic_module.Matrix(matrix).rank() == expected_rank,
                "Independent symbolic rigid-body rank mismatch")
    nullity = 6 - expected_rank
    require(receipt["first_order_degrees_of_freedom"] == nullity,
            "Incorrect first-order degree count")
    require(len(receipt["nullspace_basis"]) == nullity,
            "Missing nullspace certificate")
    vectors = [[F(x) for x in row]
               for row in receipt["nullspace_basis"]]
    for vector in vectors:
        require(len(vector) == 6 and all(sum(a*b for a, b in zip(row, vector)) == 0
                                        for row in matrix),
                "Reported motion violates a constraint")
    # These known fixtures have nullity at most two. Their Gram determinant
    # tests basis independence without reimplementing the model's row reduction.
    require(nullity <= 2, "Unexpected fixture nullity")
    if nullity == 1:
        require(sum(v*v for v in vectors[0]) > 0, "Zero nullspace certificate")
    if nullity == 2:
        aa = sum(v*v for v in vectors[0])
        bb = sum(v*v for v in vectors[1])
        ab = sum(a*b for a, b in zip(vectors[0], vectors[1]))
        require(aa*bb - ab*ab > 0, "Dependent nullspace certificate")
    expected = ("first_order_mechanism" if nullity
                else "no_first_order_rigid_body_mechanism")
    require(receipt["status"] == expected, "Kinematic interpretation changed")


def verify_all(symbolic=False):
    sp = importlib.import_module("sympy") if symbolic else None
    counts = {}
    rank_cases = [
        ("fixed", "fixed", "open", (3, 0), (-3, 0), 6),
        ("fixed", "fixed", "hinge", (3, 0), (-3, 0), 6),
        ("fixed", "pin", "open", (3, 0), (-3, 0), 5),
        ("pin", "pin", "open", (3, 0), (-3, 0), 4),
        ("pin", "pin", "hinge", (3, 0), (-3, 0), 5),
        ("pin", "pin", "moment_lock", (3, 0), (-3, 0), 6),
        ("pin", "vertical_roller", "hinge", (3, 0), (-3, 0), 5),
        ("pin", "vertical_roller", "moment_lock", (3, 0), (-3, 0), 6),
        ("vertical_roller", "vertical_roller", "moment_lock", (3, 0), (-3, 0), 5),
        ("pin", "fixed", "hinge", (3, 0), (-3, 0), 6),
        ("pin", "pin", "hinge", (3, 2), (-3, 2), 6),
        ("pin", "pin", "hinge", (3, F(1, 10**20)), (-3, 0), 6),
    ]
    kinematic_receipts = []
    for left_support, right_support, connection, left, right, rank in rank_cases:
        receipt = rigid_half_kinematics(
            left_tip_from_base=left, right_tip_from_base=right,
            left_support=left_support, right_support=right_support,
            connection=connection)
        validate_kinematic_receipt(receipt, rank, sp)
        kinematic_receipts.append(receipt)
    counts["independent_rank_and_nullspace_cases"] = len(rank_cases)

    # The collinear pinned case has a first-order motion, NOT a claim that
    # perfectly rigid fixed-length halves can finitely fold with fixed ends.
    collinear = kinematic_receipts[4]
    C = [[F(v) for v in row] for row in collinear["constraint_matrix"]]
    known_motion = [0, 0, 1, 0, 0, -1]
    require(all(sum(a*b for a, b in zip(row, known_motion)) == 0 for row in C),
            "Known equal-half counterrotation was lost")

    support_cases = [
        (10, 9, -1, 1, "inside", ((F(1, 2), F(1, 2)), (F(19, 2), F(19, 2)))),
        (10, 10, -1, 1, "marginal", ((F(0), F(0)), (F(10), F(10)))),
        (10, 11, -1, 1, "outside", ((-F(1, 2), -F(1, 2)), (F(21, 2), F(21, 2)))),
        (Interval(9, 11), Interval(8, 12), -1, 1, "indeterminate",
         ((-F(3, 2), F(3, 2)), (F(17, 2), F(23, 2)))),
        (100, 120, -2, 3, "inside", ((F(36), F(36)), (F(64), F(64)))),
        (Interval(90, 110), Interval(100, 130), -2, 3, "inside",
         ((F(28), F(46)), (F(56), F(70)))),
        (10, -10, -1, 1, "marginal", ((F(10), F(10)), (F(0), F(0)))),
        (10, -11, -1, 1, "outside", ((F(21, 2), F(21, 2)), (-F(1, 2), -F(1, 2)))),
        (10, 20, 1, 4, "inside", ((F(20, 3), F(20, 3)), (F(10, 3), F(10, 3)))),
        (10, -20, -4, -1, "inside", ((F(10, 3), F(10, 3)), (F(20, 3), F(20, 3)))),
    ]
    # Literal reactions for the two interval fixtures, ordered as
    # (Wlower,Mlower), (Wlower,Mupper), (Wupper,Mlower), (Wupper,Mupper).
    known_corner_tables = {
        3: [(F(1, 2), F(17, 2)), (-F(3, 2), F(21, 2)),
            (F(3, 2), F(19, 2)), (-F(1, 2), F(23, 2))],
        5: [(34, 56), (28, 62), (46, 64), (40, 70)],
    }
    support_receipts = []
    for index, (load, moment, left, right, expected, expected_reactions) in enumerate(support_cases):
        receipt = support_reaction_screen(total_downward_load=load,
            total_clockwise_moment=moment, x_left=left, x_right=right)
        known = known_corner_tables.get(index,
                    [(expected_reactions[0][0], expected_reactions[1][0])] * 4)
        validate_support_receipt(receipt, load=load, moment=moment, left=left, right=right,
                                 known_reactions=known, symbolic_module=sp)
        require(receipt["status"] == expected, "Known support classification changed")
        require((bounds(receipt["reaction_left"]), bounds(receipt["reaction_right"]))
                == expected_reactions, "Known support reaction changed")
        support_receipts.append(receipt)
    counts["independent_support_equilibrium_cases"] = len(support_cases)

    root = cantilever_uniform_load(line_load=3, reach=4, direction=1)
    require(bounds(root["root_vertical_force_magnitude"]) == (12, 12), "Wrong uniform shear")
    require(bounds(root["root_clockwise_moment"]) == (24, 24), "Wrong uniform root moment")
    doubled = cantilever_uniform_load(line_load=3, reach=8, direction=1)
    require(bounds(doubled["root_clockwise_moment"]) == (96, 96), "Reach-squared scaling failed")
    opposite = cantilever_uniform_load(line_load=3, reach=4, direction=-1)
    require(bounds(opposite["root_clockwise_moment"]) == (-24, -24), "Left-facing sign failed")
    root_box = cantilever_uniform_load(line_load=Interval(2, 3), reach=Interval(4, 5), direction=1)
    require(bounds(root_box["root_vertical_force_magnitude"]) == (8, 15), "Wrong shear enclosure")
    require(bounds(root_box["root_clockwise_moment"]) == (16, F(75, 2)), "Wrong root-moment enclosure")
    horizontal = horizontal_load_moment(force=Interval(-2, 4), height=Interval(3, 5))
    require(bounds(horizontal["clockwise_moment"]) == (-10, 20), "Signed force-height enclosure failed")
    retracted = cantilever_uniform_load(line_load=3, reach=0, direction=1)
    require(bounds(retracted["root_vertical_force_magnitude"]) == (0, 0)
            and bounds(retracted["root_clockwise_moment"]) == (0, 0),
            "Zero deployed reach has a nonzero own-load resultant")
    if sp is not None:
        # Independent symbolic integral of the distributed load about the root.
        w, length, x = sp.symbols("w length x", positive=True)
        require(sp.integrate(w*x, (x, 0, length)) == w*length**2/2,
                "Independent distributed-load integral failed")
    counts["load_resultant_and_scaling_cases"] = 6

    join_args = dict(roadside_separation=10, left_reach=F(11, 2), right_reach=5,
                     left_tip_height=4, right_tip_height=F(401, 100),
                     min_overlap=F(1, 4), max_overlap=F(3, 4),
                     max_vertical_mismatch=F(1, 20))
    join = join_geometry_screen(**join_args)
    require(join["status"] == "inside" and bounds(join["overlap"]) == (F(1, 2), F(1, 2))
            and bounds(join["vertical_mismatch"]) == (F(1, 100), F(1, 100)), "Known join fit changed")
    altered = {**join_args, "left_reach": 4}
    require(join_geometry_screen(**altered)["status"] == "outside", "A gap was accepted as overlap")
    altered = {**join_args, "left_reach": Interval(F(26, 5), F(27, 5))}
    require(join_geometry_screen(**altered)["status"] == "indeterminate", "Uncertain overlap was certified")
    altered = {**join_args, "left_reach": F(21, 4)}
    require(join_geometry_screen(**altered)["status"] == "marginal", "Exact lap limit lost")
    altered = {**join_args, "right_tip_height": F(41, 10)}
    require(join_geometry_screen(**altered)["status"] == "outside", "Excess vertical mismatch accepted")
    # |height difference| reaches its minimum inside a sign-changing box.
    altered = {**join_args, "left_tip_height": Interval(F(399, 100), F(401, 100)),
               "right_tip_height": 4}
    require(bounds(join_geometry_screen(**altered)["vertical_mismatch"])
            == (0, F(1, 100)), "Absolute-value interior minimum missed")
    counts["join_envelope_cases"] = 6

    ramp = ramp_slope_screen(rise=3, horizontal_run=40, declared_max_slope=F(1, 8))
    require(ramp["status"] == "inside" and bounds(ramp["slope"]) == (F(3, 40), F(3, 40))
            and bounds(ramp["required_horizontal_run"]) == (24, 24), "Known ramp arithmetic failed")
    require(ramp_slope_screen(rise=3, horizontal_run=24, declared_max_slope=F(1, 8))
            ["status"] == "marginal", "Slope equality was not marked marginal")
    require(ramp_slope_screen(rise=3, horizontal_run=20, declared_max_slope=F(1, 8))
            ["status"] == "outside", "Steeper declared slope accepted")
    require(ramp_slope_screen(rise=Interval(2, 4), horizontal_run=24,
                             declared_max_slope=F(1, 8))["status"] == "indeterminate",
            "Slope uncertainty straddle was certified")
    counts["ramp_arithmetic_cases"] = 4

    for length in range(1, len(DEPLOYMENT_RECORD_ORDER) + 1):
        require(deployment_sequence_screen(DEPLOYMENT_RECORD_ORDER[:length])["status"]
                == "consistent_record", "Valid record prefix rejected")
    for index in range(1, len(DEPLOYMENT_RECORD_ORDER)):
        skipped = [state for i, state in enumerate(DEPLOYMENT_RECORD_ORDER) if i != index]
        if index < len(DEPLOYMENT_RECORD_ORDER) - 1:
            require(deployment_sequence_screen(skipped)["status"] == "inconsistent_record",
                    "Skipped prerequisite declaration accepted")
    require(deployment_sequence_screen(["stowed", "access_requested"])["status"]
            == "inconsistent_record", "Access request skipped the declaration sequence")
    counts["record_sequence_cases"] = 13

    malformed = [
        lambda: exact(0.1), lambda: exact(True), lambda: exact("1/0"),
        lambda: Interval(2, 1), lambda: Interval(1, 2) / Interval(-1, 1),
        lambda: support_reaction_screen(total_downward_load=0, total_clockwise_moment=0, x_left=-1, x_right=1),
        lambda: support_reaction_screen(total_downward_load=10, total_clockwise_moment=0, x_left=1, x_right=1),
        lambda: cantilever_uniform_load(line_load=-1, reach=4, direction=1),
        lambda: cantilever_uniform_load(line_load=1, reach=4, direction=True),
        lambda: horizontal_load_moment(force=1, height=-1),
        lambda: join_geometry_screen(**{**join_args, "max_overlap": 0}),
        lambda: join_geometry_screen(**{**join_args, "max_vertical_mismatch": -1}),
        lambda: ramp_slope_screen(rise=1, horizontal_run=0, declared_max_slope=F(1, 8)),
        lambda: ramp_slope_screen(rise=1, horizontal_run=1, declared_max_slope=0),
        lambda: rigid_half_kinematics(left_tip_from_base=(0, 0), right_tip_from_base=(-3, 0),
                 left_support="pin", right_support="pin", connection="hinge"),
        lambda: rigid_half_kinematics(left_tip_from_base=(3, 0), right_tip_from_base=(-3, 0),
                 left_support="wheel_brakes_assumed_fixed", right_support="pin", connection="hinge"),
        lambda: deployment_sequence_screen([]), lambda: deployment_sequence_screen(["ready"]),
    ]
    for operation in malformed:
        try:
            operation()
        except (TypeError, ValueError):
            pass
        else:
            raise ArithmeticError("Malformed or unmodelled input was accepted")
    counts["malformed_inputs_rejected"] = len(malformed)

    tampered = deepcopy(support_receipts[0])
    tampered["reaction_left"]["lower"] = "1"
    changed_status = deepcopy(support_receipts[1])
    changed_status["status"] = "inside"
    changed_rank = deepcopy(collinear)
    changed_rank["constraint_rank"] = 6
    changed_motion = deepcopy(collinear)
    changed_motion["nullspace_basis"][0][0] = "1"
    attempts = [
        lambda: validate_support_receipt(tampered, load=10, moment=9, left=-1, right=1,
                    known_reactions=[(F(1, 2), F(19, 2))]*4, symbolic_module=sp),
        lambda: validate_support_receipt(changed_status, load=10, moment=10, left=-1, right=1,
                    known_reactions=[(0, 10)]*4, symbolic_module=sp),
        lambda: validate_kinematic_receipt(changed_rank, 5, sp),
        lambda: validate_kinematic_receipt(changed_motion, 5, sp),
    ]
    for operation in attempts:
        try:
            operation()
        except ArithmeticError:
            pass
        else:
            raise ArithmeticError("Tampered certificate accepted")
    counts["tampered_receipts_rejected"] = len(attempts)

    return {"checks_passed": True,
            "verification_mode": "standard_library_plus_symbolic" if symbolic else "standard_library",
            "checks": counts,
            "optional_symbolic_checks": {
                "enabled": symbolic, "rank_cases": len(rank_cases) if symbolic else 0,
                "equilibrium_corner_solves": 4*len(support_cases) if symbolic else 0,
                "distributed_load_integrals": 1 if symbolic else 0,
            },
            "examples": {
                "collinear_pins_center_hinge": collinear,
                "pin_and_vertical_roller_center_moment_lock": kinematic_receipts[7],
                "fixed_bases_center_hinge": kinematic_receipts[1],
                "support_inside": support_receipts[0],
                "support_marginal": support_receipts[1],
                "support_outside": support_receipts[2],
                "support_indeterminate": support_receipts[3],
                "cantilever": root, "join_geometry": join, "ramp": ramp,
            },
            "scope": "Synthetic exact model controls only; not a bridge design or physical approval"}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbolic", action="store_true",
                        help="add independent SymPy rank, equilibrium and integral checks")
    args = parser.parse_args()
    print(json.dumps(verify_all(symbolic=args.symbolic), indent=2))
