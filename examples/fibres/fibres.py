"""Pull a rational point of v² = u³ − 29u − 59 back to the fixed quartic C⋆."""

import argparse
from dataclasses import dataclass
from fractions import Fraction

import sympy as sp

t = sp.Symbol("t")


def quartic(x, y, z):
    """The homogeneous quartic from scan/cstar.json, grouped by powers of Z."""
    return (
        x**4 - 6*x**3*y - 27*x**2*y**2 + 54*x*y**3 + 162*y**4
        + z*(2*x**3 - 18*x**2*y + 108*y**3)
        + z**2*(-3*x**2 + 837*y**2)
        + z**3*(-31*x + 186*y)
        + 31*z**4
    )


def map_numerators(x, y, z):
    """Return numerators for u and v over 3Y² and 27Y³, respectively."""
    a = x - 6*y - z
    return -9*y**2 + z*a, 27*y**3 - 9*y*z*a - 9*y*a**2 - 3*z*a**2 - a**3


def map_at(point):
    """Evaluate the extended map on an ordinary rational source point."""
    x, y, z = map(sp.Rational, point)
    if x == y == z == 0:
        raise ValueError("Choose projective coordinates with at least one nonzero entry.")
    if quartic(x, y, z) != 0:
        raise ValueError("That point does not satisfy the equation of C⋆.")
    if y == 0:
        return (sp.Integer(7), sp.Integer(-9)) if x == z else None
    a, b = map_numerators(x, y, z)
    return a / (3*y**2), b / (27*y**3)


@dataclass(frozen=True)
class Orbit:
    """Each root of polynomial gives one point. Multiplicity applies to each root."""
    polynomial: sp.Poly
    coordinates: tuple
    multiplicity: int


def split_points(equation, coords):
    factors = sp.Poly(equation, t, domain=sp.QQ).factor_list()[1]
    return [Orbit(factor, coords, multiplicity) for factor, multiplicity in factors]


def fibre(target):
    if target is None:
        # The fourth point on Y=0 maps to (7,-9).
        return split_points(t**3 + 3*t**2 - 31, (t, 0, 1))

    u, v = target
    if v**2 != u**3 - 29*u - 59:
        raise ValueError("The target must satisfy v² = u³ − 29u − 59.")
    r = u + 3
    if r == 0 and v == 1:
        # The coordinate h used for other fibres coalesces at h=0 here.
        rational = split_points(t + 3, (t, 1, 0))
        conjugates = split_points(3*t**2 + t - 5, (6*t + 1, t, 1))
        return rational + conjugates

    equation = t**3 + 3*t**2 + r*t + v + r - 1
    if r == 0:
        return split_points(equation, (6 + 3*t, 1, 0))
    # Dividing by r keeps this formula regular even at t=0.
    return split_points(equation, (1 + (3*t**2 + 6*t)/r, t/r, 1))


def check_orbit(orbit, target):
    f = orbit.polynomial

    def vanishes(expression):
        return sp.Poly(expression, t, domain=sp.QQ).rem(f).is_zero

    x, y, z = orbit.coordinates
    if all(vanishes(c) for c in (x, y, z)) or not vanishes(quartic(x, y, z)):
        raise ArithmeticError("Source equation or projective check failed.")
    if vanishes(y):
        observed = (7, -9) if vanishes(x - z) else None
        if observed != target:
            raise ArithmeticError("Extended-map check at Y=0 failed.")
    else:
        if target is None:
            raise ArithmeticError("A finite image was found above infinity.")
        u, v = target
        a, b = map_numerators(x, y, z)
        if not vanishes(a - 3*y**2*u) or not vanishes(b - 27*y**3*v):
            raise ArithmeticError("Exact image check failed.")


def rational(text):
    try:
        return sp.Rational(Fraction(text))
    except (ValueError, ZeroDivisionError) as exc:
        raise argparse.ArgumentTypeError("Use an integer, a decimal, or a fraction such as 3/2.") from exc


def projective(coords):
    return "(" + " : ".join(str(c) for c in coords) + ")"


def show_fibre(target):
    orbits = fibre(target)
    total = sum(g.polynomial.degree()*g.multiplicity for g in orbits)
    if total != 3:
        raise ArithmeticError("The fibre degree check failed.")
    for g in orbits:
        check_orbit(g, target)

    print("Target:", "infinity" if target is None else str(target))
    rational_count = 0
    distinct = 0
    for g in orbits:
        degree = g.polynomial.degree()
        distinct += degree
        count = f"; multiplicity {g.multiplicity}" if g.multiplicity > 1 else ""
        if degree == 1:
            rational_count += 1
            value = -g.polynomial.nth(0)/g.polynomial.nth(1)
            point = [sp.sympify(c).subs(t, value) for c in g.coordinates]
            print(f"Rational point: {projective(point)}{count}")
        else:
            print(f"{degree} conjugate points: {projective(g.coordinates)}{count}")
            print(f"  One point for each root of {g.polynomial.as_expr()} = 0.")
            if degree == 2:
                print("  t =", " or ".join(str(r) for r in sp.solve(g.polynomial.as_expr(), t)))
    print(f"Distinct points: {distinct}; rational points: {rational_count}; degree counting multiplicity: {total}.")
    print("Exact source equation and image checks passed for every point.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    commands = parser.add_subparsers(dest="command")
    sub = commands.add_parser("target", help="Find all source points over a rational target (u,v).")
    sub.add_argument("u", type=rational)
    sub.add_argument("v", type=rational)
    commands.add_parser("infinity", help="Find the sources of the elliptic point at infinity.")
    sub = commands.add_parser("forward", help="Map rational projective coordinates X Y Z to the target.")
    for key in ("X", "Y", "Z"):
        sub.add_argument(key, type=rational)
    args = parser.parse_args()
    try:
        if args.command == "forward":
            result = map_at((args.X, args.Y, args.Z))
            print("Target:", "infinity" if result is None else str(result))
        elif args.command == "infinity":
            show_fibre(None)
        elif args.command == "target":
            show_fibre((args.u, args.v))
        else:
            show_fibre((sp.Integer(-3), sp.Integer(1)))
    except ValueError as exc:
        parser.error(str(exc))


if __name__ == "__main__":
    main()
