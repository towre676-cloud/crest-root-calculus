"""Find rational polynomial quotients of hyperelliptic equations from coefficients."""

import argparse
import json
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.ring_series import rs_mul, rs_nth_root
from sympy.polys.rings import ring


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "bounded-correspondences"))
from finite import coefficients, rational_points


X = sp.Symbol("x")
SIEVE_PRIMES = (1009, 1013, 1019)


def polynomial(values):
    if (not isinstance(values, list) or not values
            or any(type(value) not in (str, int) for value in values)):
        raise ValueError("Use a nonempty list of integer or rational-string coefficients.")
    return sp.Poly.from_list([sp.Rational(value) for value in reversed(values)], X, domain=sp.QQ)


def normalization(source):
    """Use the odd-multiplicity factors and keep the original leading unit."""
    if not source:
        raise ValueError("The zero radicand does not define a double cover.")
    unit, factors = source.sqf_list()
    square = sp.Poly(1, X, domain=sp.QQ)
    reduced = sp.Poly(unit, X, domain=sp.QQ)
    for factor, multiplicity in factors:
        square *= factor ** (multiplicity // 2)
        if multiplicity % 2:
            reduced *= factor
    return square, reduced


def genus(radicand):
    # Constants have no positive-genus component, even if the algebra is split.
    return max(0, (int(radicand.degree()) - 1) // 2)


def constant_degrees(source):
    n = int(source.degree())
    return [(d, n // d) for d in range(2, n // 3 + 1) if n % d == 0]


def left_divide(source, inner):
    """Write in base inner; every digit must be constant in a composition."""
    zero = source.domain.zero
    coefficients_ascending = []
    remainder = source
    while remainder:
        remainder, digit = remainder.div(inner)
        coefficients_ascending.append(digit.rep.nth(0) if digit else zero)
        if digit.degree() > 0:
            return None
    return sp.Poly.from_list(list(reversed(coefficients_ascending)), X, domain=source.domain)


def exact_inner(source, d, m):
    """Polynomial part of the normalized root at infinity, with constant zero."""
    series_ring, t = ring("t", source.domain)
    series = sum((source.rep.nth(d * m - i) * t**i for i in range(d)), series_ring.zero)
    series = series * source.domain.exquo(source.domain.one, source.rep.LC())
    root = rs_nth_root(series, m, t, d)
    descending = [root.get((i,), source.domain.zero) for i in range(d)]
    return sp.Poly.from_list(descending + [source.domain.zero], X, domain=source.domain)


def modular_inner(source, d, m, prime):
    """The normalized coefficients are integral at a good sieve prime."""
    domain = sp.GF(prime)
    leading = source.nth(d * m)
    heads = [domain.convert(int(c.p)) / domain.convert(int(c.q))
             for c in (source.nth(d * m - i) / leading for i in range(d))]
    inner = [domain.one]
    for i in range(1, d):
        numerator = sum(((i - (m + 1) * j) * heads[i - j] * inner[j]
                         for j in range(i)), domain.zero)
        inner.append(numerator / domain.convert(i * m))
    return sp.Poly.from_list(inner + [domain.zero], X, domain=domain)


def constant_maps(source, *, sieve=True):
    rows = []
    for d, m in constant_degrees(source):
        for prime in SIEVE_PRIMES if sieve else ():
            if prime <= max(d, m) or any(c.q % prime == 0 for c in source.all_coeffs()):
                continue
            if source.LC().p % prime == 0:
                continue
            modular = sp.Poly.from_list([int(c.p) * pow(int(c.q), -1, prime) % prime
                                         for c in source.all_coeffs()], X, modulus=prime)
            inner = modular_inner(source, d, m, prime)
            if left_divide(modular, inner) is None:
                rows.append({"degree": d, "target_degree": m, "excluded_modulo": prime})
                break
            # One modular pass is sufficient for the sieve; all survivors are exact.
            break
        if rows and rows[-1]["degree"] == d:
            continue
        inner = exact_inner(source, d, m)
        outer = left_divide(source, inner)
        row = {"degree": d, "target_degree": m, "inner": coefficients(inner)}
        if outer is not None:
            row["outer"] = coefficients(outer)
        rows.append(row)
    return rows


def nonconstant_degrees(source):
    n = int(source.degree())
    return [(d, m, (m * d - n) // 2) for m in range(3, n, 2)
            for d in range(3, (n - 2) // (m - 2) + 1, 2)
            if n < m * d <= n + 2 * (d - 1)]


def vertical_from_complement(source, d, m, e, domain):
    """Solve the top e equations of U' = d V W for monic V."""
    series_ring, t = ring("t", domain)
    n = source.degree()
    f = sum((source.rep.nth(n - j) * t**j for j in range(e + 1)), series_ring.zero)
    f *= domain.exquo(domain.one, source.rep.LC())
    w_coefficients = domain.gens if domain.is_PolynomialRing else ()
    w = series_ring.one + sum((c * t**j for j, c in enumerate(w_coefficients, 1)), series_ring.zero)
    v = series_ring.one
    for j in range(1, e + 1):
        inflated = rs_mul(f, rs_mul(v, v, t, j + 1), t, j + 1)
        root = rs_nth_root(inflated, m, t, j + 1).get((j,), domain.zero)
        product = rs_mul(v, w, t, j + 1).get((j,), domain.zero)
        coefficient = domain.exquo(m * ((d - j) * root - d * product),
                                    domain.convert(m * d - 2 * (d - j)))
        v += coefficient * t**j
    return (sp.Poly.from_list([v.get((j,), domain.zero) for j in range(e + 1)], X, domain=domain),
            sp.Poly.from_list([domain.one, *w_coefficients],
                              X, domain=domain))


def variable_incidence(source, d, m, e, *, derivative=True):
    """Eliminate U and G, then use the shorter of V and its derivative cofactor."""
    w = d - 1 - e
    complement = derivative and w < e
    prefix, count = ("w", w) if complement else ("v", e)
    variables = sp.symbols(f"{prefix}0:{count}")
    domain = sp.QQ.poly_ring(*variables) if variables else sp.QQ
    source = source.set_domain(domain)
    if complement:
        vertical, cofactor = vertical_from_complement(source, d, m, e, domain)
    else:
        vertical = sp.Poly.from_list([domain.one, *domain.gens], X, domain=domain)
    inflated = source * vertical**2
    inner = exact_inner(inflated, d, m)
    if derivative:
        differential = inner.diff() - d * vertical * cofactor if complement else inner.diff().rem(vertical)
        derivative_terms = differential.coeffs() if differential else []
    else:
        derivative_terms = []

    if not variables:
        valid = not derivative_terms and left_divide(inflated, inner) is not None
        return variables, [] if valid else [sp.Integer(1)], inner, vertical
    basis = sp.groebner(derivative_terms, *variables, domain=sp.QQ) if derivative_terms else None
    if basis is not None:
        if list(basis) == [1]:
            return variables, [sp.Integer(1)], inner, vertical
        inner = sp.Poly.from_list([basis.reduce(c)[1] for c in inner.all_coeffs()], X, domain=domain)

    # Division is over a polynomial coefficient ring. Since the inner
    # polynomial is monic, this commutes with every specialization.
    equations = list(basis) if basis is not None else []
    remainder = inflated
    while remainder:
        remainder, digit = remainder.div(inner)
        for (power,), expression in digit.terms():
            if power:
                residual = basis.reduce(expression)[1] if basis is not None else expression
                if residual:
                    equations.append(residual)
    return variables, equations, inner, vertical


def variable_maps(source, d, m, e, *, derivative=True):
    variables, equations, inner, vertical = variable_incidence(source, d, m, e, derivative=derivative)
    certificate = rational_points(equations, variables) if variables else None
    points = certificate["rational_points"] if variables else ([] if equations else [[]])
    maps = []
    for point in points:
        assignment = dict(zip(variables, map(sp.Rational, point), strict=True))
        u = sp.Poly(inner.as_expr().subs(assignment), X, domain=sp.QQ)
        v = sp.Poly(vertical.as_expr().subs(assignment), X, domain=sp.QQ)
        inflated = source * v**2
        g = left_divide(inflated, u)
        if (g is None or g.degree() != m or g.compose(u) != inflated
                or g.gcd(g.diff()).degree() != 0 or not u.diff().rem(v).is_zero):
            raise ArithmeticError("A recovered map failed its original identity.")
        maps.append({"inner": coefficients(u), "outer": coefficients(g), "vertical": coefficients(v)})
    return {"degree": d, "target_degree": m, "vertical_degree": e,
            "maps": sorted(maps, key=lambda row: (row["inner"], row["outer"], row["vertical"])),
            "geometric_points": certificate["geometric_points"] if variables else len(points),
            "certificate": certificate}


def discover(values, *, nonconstant=False, sieve=True):
    square, source = normalization(polynomial(values))
    if nonconstant and source.degree() % 2 != 1:
        raise ValueError("The nonconstant-multiplier search requires a normalized odd-degree source.")
    result = {"square": coefficients(square), "source": coefficients(source),
              "source_genus": genus(source), "constant": constant_maps(source, sieve=sieve)}
    if nonconstant:
        result["nonconstant"] = [variable_maps(source, d, m, e)
                                 for d, m, e in nonconstant_degrees(source)]
    return result


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON list of coefficients, constant term first")
    parser.add_argument("--nonconstant", action="store_true", help="also search with nonconstant V")
    args = parser.parse_args()
    values = json.loads(args.input.read_text(encoding="utf-8"))
    print(json.dumps(discover(values, nonconstant=args.nonconstant), indent=2))


if __name__ == "__main__":
    main()
