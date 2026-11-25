"""Certify optimal elliptic degree spectra from two explicit rational maps."""

import argparse
from collections import Counter
from fractions import Fraction
import json
from math import gcd, isqrt, lcm
from pathlib import Path
import sys

import sympy as sp

from absence import good_model
from kernel_certificate import certify_kernel, point_key


HERE = Path(__file__).resolve().parent
PERIOD_STUDY = HERE.parent / "period-reconstruction"
sys.path.insert(0, str(PERIOD_STUDY))
from certificate import certify_record, x  # noqa: E402


def load_maps(f, records):
    maps, directions = [], set()
    for record in records:
        degree, field, slope, _ = certify_record(f, record)
        if field != "Q" or slope in directions:
            continue
        directions.add(slope)
        differential = record["differential"].removeprefix("(").removesuffix(") x dx/y")
        a, b = map(Fraction, differential.split(") dx/y + ("))
        mp = dict(a=a, b=b, degree=degree)
        for key in ("N", "D"):
            poly = sp.Poly(record[key], x, domain=sp.QQ)
            mp[key] = [Fraction(str(value)) for value in reversed(poly.all_coeffs())]
        for key in ("g2", "g3"):
            mp[key] = Fraction(record[key])
        maps.append(mp)
    if len(maps) != 2:
        raise ValueError("Supply exactly two independent rational maps.")
    return maps


def j_invariant(record):
    g2, g3 = record["g2"], record["g3"]
    return 1728*g2**3/(g2**3-27*g3**2)


def trace_over_prime(record, p):
    G2, G3 = (value.numerator*pow(value.denominator, -1, p) % p
              for value in (record["g2"], record["g3"]))
    squares = Counter(y*y % p for y in range(p))
    return p-sum(squares[(4*t**3-G2*t-G3) % p] for t in range(p))


def trace_of_power(trace, p, exponent):
    a, b = 2, trace
    for _ in range(1, exponent):
        a, b = b, trace*b-p*a
    return a if not exponent else b


def nonisogeny_check(maps, p):
    if p in (2, 3) or not sp.isprime(p):
        raise ValueError("Use a good odd prime for the quotient curves.")
    for record in maps:
        g2, g3 = record["g2"], record["g3"]
        if any(value.denominator % p == 0 for value in (g2, g3)):
            raise ValueError("The quotient coefficients have poles at the witness prime.")
        if (g2**3-27*g3**2).numerator % p == 0:
            raise ValueError("The witness is a bad prime for an elliptic target.")
    traces = [trace_over_prime(record, p) for record in maps]
    if trace_of_power(traces[0], p, 120) == trace_of_power(traces[1], p, 120):
        raise ValueError("The witness does not separate the geometric factors.")
    return {"p": p, "traces": traces, "exponent": 120}


def spectrum_from_kernel(n, degrees, torsion_data, maps, bound):
    field, (E1, _), (points, images), kernel = torsion_data
    g21, g31 = maps[0]["g2"], maps[0]["g3"]
    g22, g32 = maps[1]["g2"], maps[1]["g3"]
    scale = (g32*g21)/(g31*g22)
    if scale**2*g21 != g22 or scale**3*g31 != g32:
        raise ArithmeticError("The geometric identification fails the curve identity.")
    u = field.q(scale)
    if not u.is_square():
        raise ValueError("The elliptic identification needs a larger residue field.")
    v = u*u.sqrt()
    target_points = {point_key(Q) for Q in images}

    def rho(P):
        if P is None:
            return None
        X, Y = P
        return u*X, v*Y

    if any(point_key(rho(P)) not in target_points for P in points):
        raise ArithmeticError("The identification fails on the enumerated torsion.")
    kernel_set = {(point_key(P), point_key(Q)) for P, Q in kernel}
    weights = {}
    for a in range(n):
        for b in range(n):
            if gcd(gcd(a, b), n) == 1:
                weights[a, b] = sum(
                    (point_key(E1.mul(a, P)), point_key(rho(E1.mul(b, P)))) in kernel_set
                    for P in points
                )
    # Over the algebraic closure the identification is only fixed up to sign.
    normal = tuple(weights[a, b] for a, b in weights)
    flipped = tuple(weights[a, -b % n] for a, b in weights)
    weights = dict(zip(weights, min(normal, flipped)))
    d1, d2 = degrees
    counts = Counter()
    for a in range(isqrt(d2*bound)+1):
        for b in range(-isqrt(d1*bound), isqrt(d1*bound)+1):
            if (a == 0 and b != 1) or gcd(a, b) != 1:
                continue
            d, rem = divmod(d1*a*a+d2*b*b, weights[a % n, b % n])
            if rem:
                raise ArithmeticError("The finite kernel order does not divide the polarization degree.")
            if d <= bound:
                counts[d] += 1
    factor = sp.Rational(4*sum(weights.values()), n*n)/sp.sqrt(d1*d2)
    return {"residue_K_v": {f"{a},{b}": m for (a,b), m in weights.items()},
            "spectrum": {str(d): counts[d] for d in sorted(counts)},
            "leading_times_pi": str(factor)}


def certify(f, records, p, r, witness=None, bound=150):
    if bound < 0 or r < 1:
        raise ValueError("Use a nonnegative bound and a positive extension degree.")
    maps = load_maps(f, records)
    degrees = [row["degree"] for row in maps]
    n = lcm(*degrees)
    if n not in (2, 4):
        raise ValueError("This implementation handles degree-two and degree-four input maps.")
    j = [j_invariant(record) for record in maps]
    result = {"p": p, "r": r, "j": list(map(str, j)), "degrees_of_maps": degrees}
    if j[0] == j[1]:
        if j[0].denominator == 1:
            raise ValueError("Need a non-CM certificate for the equal-j case.")
    else:
        if witness is None:
            raise ValueError("Need a good prime separating the elliptic targets.")
        result["non_isogeny"] = nonisogeny_check(maps, witness)
    data = certify_kernel(f, maps, p, r, n)
    kernel = data[-1]
    result["kernel_order"] = len(kernel)
    if len(kernel) != degrees[0]*degrees[1]:
        raise ValueError("The full product-degree kernel bound was not attained.")
    if j[0] == j[1]:
        result.update(spectrum_from_kernel(n, degrees, data, maps, bound))
    else:
        axis_sizes = [sum(Q is None for P, Q in kernel), sum(P is None for P, Q in kernel)]
        result["axis_kernel_orders"] = axis_sizes
        if any(degree % size for degree, size in zip(degrees, axis_sizes)):
            raise ArithmeticError("The axis kernel order does not divide the map degree.")
        result["optimal_degrees"] = [degree//size for degree, size in zip(degrees, axis_sizes)]
    return result


def data():
    return json.loads((PERIOD_STUDY / "replay.json").read_text(encoding="utf-8"))["curves"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    fixtures = ("C", "Cprime", "T1", "T2")
    source = parser.add_mutually_exclusive_group(required=True)
    source.add_argument("curve", nargs="?", choices=fixtures)
    source.add_argument("--input", type=Path, help="JSON object with f and maps")
    source.add_argument("--coefficients", type=json.loads, help="Propose maps from the coefficients with periods")
    parser.add_argument("--period-bound", type=int, default=6)
    parser.add_argument("--prime", type=int, required=True)
    parser.add_argument("--extension-degree", type=int, default=1)
    parser.add_argument("--non-isogeny-prime", type=int)
    parser.add_argument("--bound", type=int, default=150)
    args = parser.parse_args()
    if args.bound < 0 or args.extension_degree < 1:
        parser.error("Use a nonnegative bound and a positive extension degree.")
    if args.coefficients is not None:
        from maps import first_quotient
        record = first_quotient(args.coefficients, args.period_bound)
        record["f"] = args.coefficients
    elif args.input is not None:
        record = json.loads(args.input.read_text(encoding="utf-8"))
    else:
        record = list(data().values())[fixtures.index(args.curve)]
    try:
        result = certify(record["f"], record.get("maps", []), args.prime, args.extension_degree,
                         args.non_isogeny_prime, args.bound)
    except ValueError as error:
        parser.exit(2, f"Unresolved: {error}\n")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
