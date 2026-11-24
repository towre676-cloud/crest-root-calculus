"""Independent analytic review of the bounded engineering challenge.

The default uses only the standard library. --symbolic adds SymPy checks of
quadratic extrema and rotation matrices. The saved receipt binds the reviewed
source and report bytes; this is not a physical or code-compliance approval.
"""

import argparse
from copy import deepcopy
from fractions import Fraction as F
import hashlib
import importlib
from itertools import product
import json
from pathlib import Path

from engineering import support_polygon, minimum_quadratic, join_capture, access_envelope
from engineering_study import run as regenerate_report
from verify_engineering import verify as original_checks


HERE = Path(__file__).resolve().parent


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def _key(row):
    return (F(row["half_reach_m"]), F(row["support_half_width_across_m"]),
            F(row["support_half_length_along_m"]), F(row["ballast_weight_kN"]))


def check_candidate_report(report):
    """Independent closed-form endpoint/wind inequalities for the actual grid."""
    reaches = (F(19, 2), F(23, 2), F(27, 2))
    expected_grid = set(product(reaches, (F(3, 2), F(2), F(5, 2)),
                                (F(3, 2), F(5, 2)), map(F, range(0, 101, 10))))
    observed = [_key(row) for row in report["candidates"]]
    require(len(observed) == len(set(observed)) == 198 and set(observed) == expected_grid,
            "Candidate grid omitted, duplicated, or changed a design")
    feasible_count, edge_checks = 0, 0
    for row in report["candidates"]:
        r, a, b, ballast = _key(row)
        W = 130 + 4*r + ballast
        # The deck's affine CG increases with s. All weights remain constant.
        qx_low = -4*r - F(5, 2)*ballast - 60
        qx_high = 2*r*r - F(5, 2)*ballast + 60
        expected = -a*W < qx_low and qx_high < a*W and 36 < b*W
        require(row["strict_vertical_contact_feasibility"] is expected,
                "Closed-form whole-path/wind feasibility differs")
        require(F(row["support_rectangle_area_m2"]) == 4*a*b, "Support area changed")
        feasible_count += expected
        wind_rows = row["wind_corners"]
        require(len(wind_rows) == 4 and {(F(c["wind_x_kN"]), F(c["wind_y_kN"]))
                    for c in wind_rows} == set(product((-10, 10), (-6, 6))),
                "Wind rectangle corners are incomplete")
        for corner in wind_rows:
            wx, wy = F(corner["wind_x_kN"]), F(corner["wind_y_kN"])
            q0, q1, qy = -4*r-F(5, 2)*ballast+6*wx, 4*r+2*r*r, 6*wy
            # CCW hull edges: bottom, right, top, left. These follow directly
            # from x in [-a,a], y in [-b,b], independently of the hull code.
            polynomials = [(2*a*(qy+b*W), F(0), F(0)),
                           (2*b*(a*W-q0), -2*b*q1, F(0)),
                           (2*a*(b*W-qy), F(0), F(0)),
                           (2*b*(a*W+q0), 2*b*q1, F(0))]
            result = corner["result"]
            require(F(result["total_downward_load"]) == W
                    and list(map(F, result["first_moment_x_polynomial"])) == [q0, q1, F(0)]
                    and list(map(F, result["first_moment_y_polynomial"])) == [qy, F(0), F(0)],
                    "First-moment path coefficients differ from direct statics")
            require(len(result["edges"]) == 4, "Rectangular support edge count changed")
            all_minima = []
            for edge, coefficients in zip(result["edges"], polynomials):
                minimum = min(coefficients[0], sum(coefficients))
                argmin = F(1) if coefficients[1] < 0 else F(0)
                require(tuple(map(F, edge["polynomial"])) == coefficients
                        and F(edge["minimum"]) == minimum and F(edge["at"]) == argmin,
                        "Reported edge polynomial, minimum, or minimizing parameter changed")
                all_minima.append(minimum)
                edge_checks += 1
            smallest = min(all_minima)
            status = "outside" if smallest < 0 else "boundary" if smallest == 0 else "inside"
            require(result["status"] == status, "Corner path status differs from its minima")
    require(feasible_count == report["strict_feasible_count"] == 182,
            "Feasible-count receipt changed")
    # Hand-derived complete front, including the smallest ballast grid value
    # strictly above each analytic threshold. This checks completeness, not
    # merely that the displayed designs are individually nondominated.
    expected_front = {(F(19,2), F(3,2), F(3,2), F(0)),
                      (F(23,2), F(3,2), F(3,2), F(20)),
                      (F(23,2), F(2), F(3,2), F(0)),
                      (F(27,2), F(3,2), F(3,2), F(40)),
                      (F(27,2), F(2), F(3,2), F(20)),
                      (F(27,2), F(5,2), F(3,2), F(0))}
    front = report["pareto_front_under_declared_model"]
    require(len(front) == 6 and {_key(row) for row in front} == expected_front,
            "Finite-grid Pareto front is incomplete or incorrect")
    return {"candidates": 198, "wind_corner_paths": 792, "edge_minima": edge_checks,
            "feasible": feasible_count, "complete_grid_pareto_points": 6}


def hamilton(a, b):
    """Quaternion multiplication, independent of the rotation-matrix formula."""
    w, x, y, z = a
    v, p, q, r = b
    return (w*v-x*p-y*q-z*r, w*p+x*v+y*r-z*q,
            w*q-x*r+y*v+z*p, w*r+x*q-y*p+z*v)


def review(symbolic=False):
    files = ["engineering.py", "engineering_study.py", "verify_engineering.py",
             "engineering-report.json", "verification-engineering.json", "site_geometry.py",
             "verify_engineering_independent.py", "truss.py", "verify_truss.py",
             "verification-truss.json", "verification-truss-symbolic.json"]
    hashes = {name: hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files}
    report = json.loads((HERE/"engineering-report.json").read_text())
    require(report == regenerate_report(), "Saved engineering report differs from fresh execution")
    base = original_checks()
    require(base["status"] == "passed", "Original verification did not pass")
    counts = {"original_verifier_checks_replayed": base["check_count"],
              "full_report_replay": 1, "grid": check_candidate_report(report)}
    sp = importlib.import_module("sympy") if symbolic else None

    # Explicit 3D cross-product sign check: 7 downward at (2,-3), and
    # horizontal (4,-5) at height 6. External Mx=51, My=38, so the
    # vertical reaction first moments must be Qx=38, Qy=-51.
    support = support_polygon([[-10,-10],[10,-10],[10,10],[-10,10]], 7, 38, -51)
    require(support["status"] == "inside", "Known signed-moment resultant escaped support hull")
    witness = support["feasible_reaction_witness"]
    require(sum(row["reaction"] for row in witness) == 7, "Vertical force witness changed")
    require(sum(row["point"][1]*row["reaction"] for row in witness)+51 == 0
            and -sum(row["point"][0]*row["reaction"] for row in witness)+38 == 0,
            "Independent spatial torque balance exposes a sign error")
    for qx, qy, expected in ((6,6,"inside"),(18,18,"boundary"),(24,24,"outside")):
        result = support_polygon([[0,0],[6,0],[0,6],[0,0]], 6, qx, qy)
        require(result["status"] == expected, "Nonrectangular support classification differs")
    counts["support_moment_and_nonrectangular_cases"] = 4

    extrema = [([1,-4,4], (0,1), F(0), F(1,2)),
               ([1,4,-4], (0,1), F(1), F(0)),
               ([0,1,1], (0,1), F(0), F(0)),
               ([4], (0,1), F(4), F(0)),
               ([1,-4,1], (0,1), -F(2), F(1)),
               ([3,-3,2], (0,1), F(15,8), F(3,4)),
               ([2,-1,3], (-2,3), F(23,12), F(1,6)),
               ([0,0,-1], (-2,3), -F(9), F(3))]
    for coefficients, interval, value, at in extrema:
        result = minimum_quadratic(coefficients, interval)
        require(result["minimum"] == value and result["at"] == at, "Analytic quadratic extremum differs")
        if sp is not None:
            s = sp.Symbol("s", real=True)
            polynomial = sum(sp.Rational(c)*s**i for i, c in enumerate(coefficients))
            candidates = [sp.Rational(x) for x in interval]
            candidates += [root for root in sp.solve(sp.diff(polynomial,s),s)
                           if interval[0] <= root <= interval[1]]
            require(min(polynomial.subs(s,x) for x in candidates) == sp.Rational(value),
                    "Independent symbolic polynomial minimization differs")
    counter = report["endpoint_only_counterexample"]
    require(F(counter["endpoint_CoP_x"]) == F(3,4)
            and F(counter["interior_CoP_x_at_half"]) == F(5,4)
            and min(F(e["minimum"]) for e in counter["result"]["edges"]) == -F(1,2),
            "Whole-path counterexample lost its interior violation")
    counts["hand_derived_quadratic_extrema"] = len(extrema)
    counts["interior_path_violation_counterexample"] = 1

    quaternions = [(1,0,0,0),(0,0,0,1),(F(3,5),F(4,5),0,0),
                   (F(3,5),0,F(4,5),0),(F(3,5),0,0,F(4,5)),
                   (F(1,2),F(1,2),F(1,2),F(1,2)),
                   (F(39999,40001),0,F(400,40001),0)]
    points = [(1,2,3),(-3,F(1,5),2),(0,-1,-2),(0,1,2)]
    translation = (F(1,100),-F(1,200),F(3,400))
    for quaternion in quaternions:
        case = {"relative_rotation_quaternion": list(quaternion), "relative_translation": list(translation),
                "capture_half_widths": [10,10,10],
                "connectors": [{"id":str(i),"point":list(point)} for i,point in enumerate(points)]}
        result = join_capture(case)
        conjugate = (quaternion[0], *(-x for x in quaternion[1:]))
        for point, row in zip(points, result["connectors"]):
            rotated = hamilton(hamilton(quaternion,(0,*point)),conjugate)
            delta = [translation[i]+rotated[i+1]-point[i] for i in range(3)]
            require(row["displacement"] == delta and row["margins"] == [10-abs(v) for v in delta],
                    "Hamilton-product connector transform differs")
        opposite = join_capture({**case,"relative_rotation_quaternion":[-x for x in quaternion]})
        require(opposite == result, "Quaternion sign changed the physical rotation")
        R = result["rotation_matrix"]
        require(all(sum(R[k][i]*R[k][j] for k in range(3)) == (1 if i==j else 0)
                    for i in range(3) for j in range(3)), "Rotation is not orthogonal")
        if sp is not None:
            require(sp.Matrix(R).det() == 1, "Symbolic determinant does not describe a proper rotation")
    for example in report["join_examples"]:
        require(example["result"] == json.loads(json.dumps(join_capture(example["input"]), default=str)),
                "Stored connector example changed")
    counts["finite_quaternion_cases"] = len(quaternions)
    counts["hamilton_connector_transforms"] = len(quaternions)*len(points)
    counts["quaternion_sign_invariance_cases"] = len(quaternions)

    access = access_envelope("5.2832","0.1","0.6","0.1",0,"3/40")
    require(access["rise"] == F(3802,625) and access["sloped_horizontal_run"] == F(30416,375)
            and access["runs"] == 8 and access["intermediate_landings"] == 7,
            "Ramp rise, slope, or integer run count differs")
    require(access["straight_envelope_length"] == F(71119,750)
            and access["straight_envelope_width"] == F(9,5)
            and access["unstacked_switchback_envelope_length"] == F(1651,125)
            and access["unstacked_switchback_envelope_width"] == F(72,5),
            "Ramp bounding dimensions differ from direct geometry")
    for rise, expected_runs in ((F(1),1),(F(2),2),(F(2001,1000),3)):
        shape = access_envelope(rise,0,0,0,0,F(1,10),max_run_rise=1,
                                landing_length=F(3,2),outside_lane_width=2)
        require(shape["runs"] == expected_runs
                and shape["straight_envelope_length"] == 10*rise+(expected_runs+1)*F(3,2),
                "Ramp landing count failed an exact-boundary case")
    counts["ramp_dimension_and_run_boundary_cases"] = 4

    changed = deepcopy(report); changed["candidates"][0]["strict_vertical_contact_feasibility"] = False
    missing_front = deepcopy(report); missing_front["pareto_front_under_declared_model"].pop()
    changed_minimum = deepcopy(report)
    edge = changed_minimum["candidates"][0]["wind_corners"][0]["result"]["edges"][0]
    edge["minimum"] = str(F(edge["minimum"])+1)
    for altered in (changed, missing_front, changed_minimum):
        try: check_candidate_report(altered)
        except ArithmeticError: pass
        else: raise ArithmeticError("Altered engineering report accepted by independent audit")
    counts["altered_report_controls_rejected"] = 3
    require(hashes == {name:hashlib.sha256((HERE/name).read_bytes()).hexdigest() for name in files},
            "A source or reviewed report changed during independent verification")
    return {
        "status":"passed", "checks_passed":True,
        "verification_mode":"standard_library_plus_symbolic" if symbolic else "standard_library",
        "reviewed_sha256":hashes, "checks":counts,
        "optional_symbolic_checks":{"enabled":symbolic,"quadratic_extrema":len(extrema) if symbolic else 0,
                                    "proper_rotation_determinants":len(quaternions) if symbolic else 0},
        "substantive_bugs_found":[], "reviewed_source_changes_made":False,
        "review_notes":[
            "First moments Qx,Qy are not Cartesian torque components: vertical balance requires external Mx=-Qy, My=Qx. The implementation's weight and horizontal-force signs are consistent.",
            "Support-hull interior certifies existence of compressive vertical reactions, not actual pad forces, stiffness, capacities, friction or complete 3D equilibrium.",
            "Every declared quadratic path parameter is covered by endpoint/stationary-point minimization; the included counterexample correctly exposes endpoint-only checking.",
            "Four wind corners suffice because every edge margin is affine in the independent horizontal-force components at each parameter value. The force bounds are hypotheses, not a safe wind-speed determination.",
            "The six reported Pareto points are complete only for the 198 enumerated designs and the two stated objectives; no continuous or physical design optimum is established.",
            "Connector capture uses an exact active finite rigid transform in one declared reference frame. Identical nominal connector correspondence and the transform origin remain input assumptions; lock engagement, strength and elasticity are absent.",
            "Ramp results are horizontal bounding footprints with end landings. The declared run-rise/landing/width values and unstacked arrangement do not establish route, headroom, detailing or accessibility approval.",
            "Geometry predicates and site data receive a separate site audit; this review focuses on the engineering math and saved engineering report.",
        ],
        "scope":"Independent exact model/code review only; no physical validation or load rating",
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--symbolic",action="store_true",help="add optional independent SymPy cross-checks")
    args = parser.parse_args()
    receipt = review(symbolic=args.symbolic)
    (HERE/"verification-engineering-independent.json").write_text(json.dumps(receipt,indent=2)+"\n")
    print(json.dumps(receipt,indent=2))
