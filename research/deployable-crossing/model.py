"""Exact research screens for a declared deployable-crossing model.

This module is an illustrative feasibility kernel, not a structural design,
code-compliance check, deployment controller, or physical-use approval. Inputs
must describe one consistent planar model and use a consistent length/force
unit system. No force, clearance, allowable slope, or joining tolerance is
supplied by default. Exact rational arithmetic does not make assumptions true.

Public API:
  Interval(lower, upper): exact independent input enclosure; floats rejected.
  rigid_half_kinematics(...): ideal rigid-body constraint rank and nullspace.
  cantilever_uniform_load(...): horizontal half's own root shear and moment.
  horizontal_load_moment(...): signed force times elevation.
  support_reaction_screen(...): two ideal ground-contact reactions.
  join_geometry_screen(...): declared overlap and vertical-mismatch bounds.
  ramp_slope_screen(...): rise/run against a caller-supplied slope bound.
  deployment_sequence_screen(events): order of an abstract declared record.

All public results are JSON-serializable and exact numbers are rational strings.
Inequality status meanings: 'inside' means every enclosed value satisfies the
strict inequalities; 'outside' means a stated inequality is violated throughout
the enclosure; 'marginal' means a relevant margin is identically zero; otherwise
the answer is 'indeterminate'. No status establishes structural readiness.
"""

from dataclasses import dataclass
from fractions import Fraction
from collections.abc import Sequence


def exact(value):
    """Accept exact rationals only; decimal strings are exact, floats are not."""
    if isinstance(value, bool) or not isinstance(value, (int, str, Fraction)):
        raise TypeError("Use an int, Fraction, or exact rational/decimal string")
    try:
        return Fraction(value)
    except (ValueError, ZeroDivisionError) as error:
        raise ValueError("Invalid exact rational value") from error


@dataclass(frozen=True)
class Interval:
    """A closed rational enclosure, interpreted independently of other inputs."""

    lower: Fraction
    upper: Fraction

    def __post_init__(self):
        lower, upper = exact(self.lower), exact(self.upper)
        if lower > upper:
            raise ValueError("Interval lower endpoint exceeds upper endpoint")
        object.__setattr__(self, "lower", lower)
        object.__setattr__(self, "upper", upper)

    @classmethod
    def point(cls, value):
        value = exact(value)
        return cls(value, value)

    def as_dict(self):
        return {"lower": str(self.lower), "upper": str(self.upper)}

    def __add__(self, other):
        other = _interval(other)
        return Interval(self.lower + other.lower, self.upper + other.upper)

    __radd__ = __add__

    def __neg__(self):
        return Interval(-self.upper, -self.lower)

    def __sub__(self, other):
        return self + (-_interval(other))

    def __rsub__(self, other):
        return _interval(other) + (-self)

    def __mul__(self, other):
        other = _interval(other)
        corners = [a * b for a in (self.lower, self.upper)
                   for b in (other.lower, other.upper)]
        return Interval(min(corners), max(corners))

    __rmul__ = __mul__

    def __truediv__(self, other):
        other = _interval(other)
        if other.lower <= 0 <= other.upper:
            raise ValueError("Cannot divide by an interval containing zero")
        reciprocal = Interval(1 / other.upper, 1 / other.lower)
        return self * reciprocal

    def square(self):
        endpoints = [self.lower * self.lower, self.upper * self.upper]
        lower = Fraction(0) if self.lower <= 0 <= self.upper else min(endpoints)
        return Interval(lower, max(endpoints))

    def absolute(self):
        endpoints = [abs(self.lower), abs(self.upper)]
        lower = Fraction(0) if self.lower <= 0 <= self.upper else min(endpoints)
        return Interval(lower, max(endpoints))


def _interval(value):
    return value if isinstance(value, Interval) else Interval.point(value)


def _nonnegative(value, name, strictly=False):
    value = _interval(value)
    if value.lower < 0 or (strictly and value.lower == 0):
        raise ValueError(name + (" must be strictly positive" if strictly
                                 else " must be nonnegative"))
    return value


def _margin_status(value):
    if value.lower > 0:
        return "inside"
    if value.upper < 0:
        return "outside"
    if value.lower == value.upper == 0:
        return "marginal"
    return "indeterminate"


def _screen(margins):
    checks = {name: {"margin": value.as_dict(),
                     "status": _margin_status(value)}
              for name, value in margins.items()}
    states = {entry["status"] for entry in checks.values()}
    status = next((s for s in ("outside", "indeterminate", "marginal")
                   if s in states), "inside")
    return {"status": status, "checks": checks}


def _rref(matrix):
    matrix = [list(row) for row in matrix]
    pivot_columns = []
    row = 0
    for column in range(6):
        pivot = next((i for i in range(row, len(matrix))
                      if matrix[i][column] != 0), None)
        if pivot is None:
            continue
        matrix[row], matrix[pivot] = matrix[pivot], matrix[row]
        scale = matrix[row][column]
        matrix[row] = [entry / scale for entry in matrix[row]]
        for i in range(len(matrix)):
            if i != row:
                scale = matrix[i][column]
                matrix[i] = [a - scale * b
                             for a, b in zip(matrix[i], matrix[row])]
        pivot_columns.append(column)
        row += 1
        if row == len(matrix):
            break
    return matrix, pivot_columns


def _vector(value, name):
    if (not isinstance(value, Sequence) or isinstance(value, (str, bytes))
            or len(value) != 2):
        raise ValueError(name + " must have exactly two rational coordinates")
    result = tuple(exact(v) for v in value)
    if result == (0, 0):
        raise ValueError(name + " must be nonzero")
    return result


def rigid_half_kinematics(*, left_tip_from_base, right_tip_from_base,
                          left_support, right_support, connection):
    """Screen six planar rigid-half DOFs: (ux, uz, theta) for each half.

    Supports are 'pin', 'fixed', or 'vertical_roller'; connection is 'open', 'hinge', or
    'moment_lock'. A closed connection assumes both vectors refer to the same
    nominal joint; endpoint fit is a separate geometry question. Fixed means
    all three ideal rigid-body DOFs are restrained, not that outriggers or
    wheel brakes have established that restraint. A full rank does not check
    member flexibility, prestress, buckling, lock strength, or foundation.
    """
    left = _vector(left_tip_from_base, "left_tip_from_base")
    right = _vector(right_tip_from_base, "right_tip_from_base")
    support_dofs = {"pin": (0, 1), "fixed": (0, 1, 2), "vertical_roller": (1,)}
    if left_support not in support_dofs or right_support not in support_dofs:
        raise ValueError("Support must be 'pin', 'fixed', or 'vertical_roller'")
    if connection not in {"open", "hinge", "moment_lock"}:
        raise ValueError("Connection must be 'open', 'hinge', or 'moment_lock'")
    rows = []
    for offset, support in ((0, left_support), (3, right_support)):
        for local in support_dofs[support]:
            row = [Fraction(0)] * 6
            row[offset + local] = Fraction(1)
            rows.append(row)
    if connection != "open":
        rows.extend([
            [Fraction(1), Fraction(0), -left[1], Fraction(-1), Fraction(0), right[1]],
            [Fraction(0), Fraction(1), left[0], Fraction(0), Fraction(-1), -right[0]],
        ])
    if connection == "moment_lock":
        rows.append([Fraction(0), Fraction(0), Fraction(1), Fraction(0),
                     Fraction(0), Fraction(-1)])
    reduced, pivots = _rref(rows)
    basis = []
    for column in range(6):
        if column not in pivots:
            vector = [Fraction(0)] * 6
            vector[column] = Fraction(1)
            for row, pivot in enumerate(pivots):
                vector[pivot] = -reduced[row][column]
            basis.append([str(v) for v in vector])
    return {
        "status": ("first_order_mechanism" if basis
                   else "no_first_order_rigid_body_mechanism"),
        "constraint_rank": len(pivots),
        "first_order_degrees_of_freedom": 6 - len(pivots),
        "constraint_matrix": [[str(v) for v in row] for row in rows],
        "nullspace_basis": basis,
        "scope": "Ideal rigid-body kinematics only; not a finite-motion or stability proof",
    }


def cantilever_uniform_load(*, line_load, reach, direction):
    """Own-load root resultants of a horizontal cantilever before connection.

    line_load: force/length; reach: length; direction: +1 (right) or -1 (left).
    The signed clockwise load-moment demand is direction*w*reach**2/2; the
    ideal restraint's reaction couple has the opposite sign. This is NOT
    the aggregate ground-overturning moment: do not add it a second time when
    the same load is already represented by its global position in that sum.
    This calculation does not assert a root can supply the computed moment.
    """
    if type(direction) is not int or direction not in {-1, 1}:
        raise ValueError("direction must be the integer +1 or -1")
    load = _nonnegative(line_load, "line_load")
    reach = _nonnegative(reach, "reach")
    return {
        "root_vertical_force_magnitude": (load * reach).as_dict(),
        "root_clockwise_moment": (direction * load * reach.square() / 2).as_dict(),
        "scope": "Separate cantilever root demand; not a chassis moment total",
    }


def horizontal_load_moment(*, force, height):
    """Clockwise external moment F*h; F positive to the right, h above origin."""
    force = _interval(force)
    height = _nonnegative(height, "height")
    return {"clockwise_moment": (force * height).as_dict()}


def support_reaction_screen(*, total_downward_load, total_clockwise_moment,
                            x_left, x_right):
    """Two ideal vertical, nonadhesive point supports in a single plane.

    W includes ALL declared downward loads. M is their total signed clockwise
    external moment about x=0, including downward loads times their horizontal
    positions and horizontal loads times elevation. Both are caller-declared;
    no omissions or double-counting can be detected here. Positive M shifts
    reaction toward increasing x. The independent (W,M) enclosure may be
    conservative when the original loads are correlated.

    Exact equilibrium: RL+RR=W; xL*RL+xR*RR=M. Strict positive contact is the
    screen. Zero reaction is marginal, never an inside result. This does not
    check friction, soil bearing, support capacity, lateral restraint, dynamic
    loading, or whether the ideal point supports are physically present.
    """
    load = _nonnegative(total_downward_load, "total_downward_load", strictly=True)
    moment = _interval(total_clockwise_moment)
    left, right = exact(x_left), exact(x_right)
    if left >= right:
        raise ValueError("x_left must be strictly below x_right")
    span = right - left
    reaction_left = (right * load - moment) / span
    reaction_right = (moment - left * load) / span
    result = _screen({"left_contact": reaction_left, "right_contact": reaction_right})
    result.update(reaction_left=reaction_left.as_dict(),
                  reaction_right=reaction_right.as_dict(),
                  equivalent_condition="x_left*W < M < x_right*W",
                  scope="Declared aggregate-load, two-contact planar screen only")
    return result


def join_geometry_screen(*, roadside_separation, left_reach, right_reach,
                         left_tip_height, right_tip_height, min_overlap,
                         max_overlap, max_vertical_mismatch):
    """Declared horizontal lap-joint envelope, not a lock-engagement proof.

    All reaches, heights, overlaps and mismatches have length units. Positive
    overlap = left_reach+right_reach-separation; negative means a gap. Heights
    refer to one common datum. This screen does not model angular mismatch,
    transverse mismatch, connection internals, motion, or connection strength.
    Independent interval arithmetic encloses all stated input combinations.
    """
    separation = _nonnegative(roadside_separation, "roadside_separation", strictly=True)
    left = _nonnegative(left_reach, "left_reach")
    right = _nonnegative(right_reach, "right_reach")
    height_left, height_right = _interval(left_tip_height), _interval(right_tip_height)
    minimum, maximum, tolerance = map(exact, (min_overlap, max_overlap, max_vertical_mismatch))
    if minimum < 0 or maximum < minimum or tolerance < 0:
        raise ValueError("Require 0 <= min_overlap <= max_overlap and nonnegative mismatch tolerance")
    overlap = left + right - separation
    mismatch = (height_left - height_right).absolute()
    result = _screen({"minimum_overlap": overlap - minimum,
                      "maximum_overlap": maximum - overlap,
                      "vertical_mismatch": tolerance - mismatch})
    result.update(overlap=overlap.as_dict(), vertical_mismatch=mismatch.as_dict(),
                  scope="Declared lap and vertical fit only; no lock or strength inference")
    return result


def ramp_slope_screen(*, rise, horizontal_run, declared_max_slope):
    """Exact rise/horizontal-run arithmetic, with no accessibility-code default.

    Both rise and run have length units; slope is dimensionless. This checks
    neither landings nor width, cross-slope, guarding, route feasibility or any
    jurisdiction's compliance. Required run is rise/declared_max_slope.
    """
    rise = _nonnegative(rise, "rise")
    run = _nonnegative(horizontal_run, "horizontal_run", strictly=True)
    maximum = exact(declared_max_slope)
    if maximum <= 0:
        raise ValueError("declared_max_slope must be strictly positive")
    slope = rise / run
    result = _screen({"declared_slope_margin": maximum - slope})
    result.update(slope=slope.as_dict(),
                  required_horizontal_run=(rise / maximum).as_dict(),
                  scope="Declared rise/run arithmetic only; no compliance inference")
    return result


DEPLOYMENT_RECORD_ORDER = (
    "stowed", "supports_declared", "extension_in_progress", "alignment_declared",
    "lock_engaged_declared", "lock_verified_declared", "access_requested",
)


def deployment_sequence_screen(events):
    """Check a prefix of an abstract forward deployment declaration record.

    These strings assert no facts. A consistent record does not verify any
    physical support, clearance, connection, inspection, or permission to use
    a bridge. 'access_requested' is a request, never an authorization. Actual
    controls, emergency/retraction sequences and fault handling are outside
    this example. Skipping a declaration is detected, not repaired.
    """
    if (not isinstance(events, Sequence) or isinstance(events, (str, bytes))
            or not events):
        raise ValueError("events must be a nonempty sequence of known state strings")
    for state in events:
        if not isinstance(state, str) or state not in DEPLOYMENT_RECORD_ORDER:
            raise ValueError("Unknown deployment-record state")
    for index, state in enumerate(events):
        if index >= len(DEPLOYMENT_RECORD_ORDER) or state != DEPLOYMENT_RECORD_ORDER[index]:
            return {"status": "inconsistent_record", "first_invalid_index": index,
                    "scope": "Record-order check only; no physical readiness inference"}
    return {"status": "consistent_record", "recorded_states": len(events),
            "scope": "Record-order check only; no physical readiness inference"}
