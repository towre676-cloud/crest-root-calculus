"""Exact unknown-target searches for cyclic covers over the rationals."""

import argparse
import json
from math import gcd
from pathlib import Path
import sys

import sympy as sp
from sympy.polys.ring_series import rs_mul, rs_nth_root, rs_pow
from sympy.polys.rings import ring


HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "coefficient-quotients"))
from discover import X, coefficients, exact_inner, left_divide, polynomial
from finite import rational_points, verify as verify_finite


def power_reduce(source, exponent, power=1):
    """Factor source**power = removed**exponent * branch without full factorization."""
    if not source:
        raise ValueError("The branch polynomial must be nonzero.")
    unit, factors = source.sqf_list()
    removed = sp.Poly(1, X, domain=sp.QQ)
    branch = sp.Poly(unit**power, X, domain=sp.QQ)
    for factor, multiplicity in factors:
        quotient, remainder = divmod(multiplicity * power, exponent)
        removed *= factor**quotient
        branch *= factor**remainder
    return removed, branch


def twice_genus(exponent, factors):
    return 1 - exponent + sum(f.degree()*(exponent-gcd(exponent, i)) for f, i in factors)


def source_model(values, exponent):
    if type(exponent) is not int or exponent < 2:
        raise ValueError("The source exponent must be an integer at least two.")
    removed, branch = power_reduce(polynomial(values), exponent)
    if branch.degree() < 2 or gcd(exponent, branch.degree()) != 1:
        raise ValueError("The normalized branch degree must be at least two and coprime to the source exponent.")
    _, factors = branch.sqf_list()
    # The vertical eigenspaces have distinct residues of their pole orders.
    first_noninvariant = min(sum(f.degree()*(i*j % exponent) for f, i in factors)
                             for j in range(1, exponent))
    return removed, branch, twice_genus(exponent, factors) // 2, first_noninvariant


def eligible(branch, k, minimum=2, maximum=None):
    """List (deg U, deg G, deg V) for one normalized character branch."""
    n = branch.degree()
    genus_twice = twice_genus(k, branch.sqf_list()[1])
    if minimum == 1 and genus_twice >= 2 and genus_twice == (k-1)*(n-1):
        yield 1, n, 0
    for m in range(2, n):
        slope = (k-1)*m-k
        if gcd(k, m) != 1 or slope <= 0:
            continue
        bound = (genus_twice-1)//slope
        if maximum is not None:
            bound = min(bound, maximum)
        for d in range(2, bound+1):
            e, remainder = divmod(m*d-n, k)
            if not remainder and 0 <= (k-1)*e <= d-1:
                yield d, m, e


def complementary_vertical(source, k, d, m, e, domain, names):
    series, t = ring("t", domain)
    n = source.degree()
    f = sum((source.rep.nth(n-j)*t**j for j in range(e+1)), series.zero)
    f *= domain.exquo(domain.one, source.rep.LC())
    w = series.one + sum((c*t**j for j, c in enumerate(names, 1)), series.zero)
    v = series.one
    for j in range(1, e+1):
        inflated = rs_mul(f, rs_pow(v, k, t, j+1), t, j+1)
        root = rs_nth_root(inflated, m, t, j+1).get((j,), domain.zero)
        cross = rs_mul(rs_pow(v, k-1, t, j+1), w, t, j+1).get((j,), domain.zero)
        denom = m*d*(k-1) - k*(d-j)
        coeff = domain.exquo(m*((d-j)*root-d*cross), domain.convert(denom))
        v += coeff*t**j
    return (sp.Poly.from_list([v.get((j,), domain.zero) for j in range(e+1)], X, domain=domain),
            sp.Poly.from_list([domain.one, *names], X, domain=domain))


def basic_equations(equations):
    equations = [term for term in equations if term != 0]
    if any(term.is_number for term in equations):
        return [sp.Integer(1)]
    return list(dict.fromkeys(equations))


def discriminant_cut(equations, variables, outer):
    """Retain exactly the squarefree-target open locus."""
    if equations == [1]:
        return variables, equations
    disc = outer.discriminant()
    if disc.is_number:
        return variables, equations if disc else [sp.Integer(1)]
    inverse = sp.Symbol("disc_inverse")
    return (*variables, inverse), [*equations, inverse*disc - 1]


def incidence(source, k, d, m, e, *, derivative=True):
    """Build the finite incidence, optionally omitting all derivative reductions."""
    w = d-1-(k-1)*e
    complement = derivative and w < e
    prefix, count = ("w", w) if complement else ("v", e)
    parameters = sp.symbols(f"{prefix}0:{count}")
    domain = sp.QQ.poly_ring(*parameters) if parameters else sp.QQ
    names = domain.gens if parameters else ()
    lifted = source.set_domain(domain)
    if complement:
        vertical, residual_factor = complementary_vertical(lifted, k, d, m, e, domain, names)
    else:
        vertical = sp.Poly.from_list([domain.one, *names], X, domain=domain)
    inflated = lifted*vertical**k
    inner = exact_inner(inflated, d, m)
    equations = []
    if derivative:
        derivative_power = vertical**(k-1)
        difference = (inner.diff() - d*derivative_power*residual_factor
                      if complement else inner.diff().rem(derivative_power))
        equations = basic_equations(difference.coeffs())
    if equations == [1]:
        return parameters, equations, inner, vertical, count
    if not parameters:
        outer = left_divide(inflated, inner)
        valid = outer is not None and outer.discriminant() != 0
        return parameters, [] if valid else [sp.Integer(1)], inner, vertical, count

    basis = sp.groebner(equations, *parameters, domain=sp.QQ) if equations else None
    if basis is not None:
        if list(basis) == [1]:
            return parameters, [sp.Integer(1)], inner, vertical, count
        equations = list(basis)
        inner = sp.Poly.from_list([basis.reduce(c)[1] for c in inner.all_coeffs()], X, domain=domain)
    outer_digits = []
    remainder = inflated
    while remainder:
        remainder, digit = remainder.div(inner)
        outer_digits.append(digit.nth(0))
        equations.extend(basis.reduce(c)[1] if basis is not None else c for (j,), c in digit.terms() if j)
    equations = basic_equations(equations)
    variables = parameters
    # With the derivative and a squarefree source, the degree bound implies
    # a squarefree target. Mixed source or the ablation needs this open cut.
    if not derivative or source.gcd(source.diff()).degree() > 0:
        outer = sp.Poly.from_list(list(reversed(outer_digits)), X, domain=domain)
        variables, equations = discriminant_cut(equations, parameters, outer)
    return variables, equations, inner, vertical, count


def exported_map(source, k, m, inner, vertical, assignment):
    u = sp.Poly(inner.as_expr().subs(assignment), X, domain=sp.QQ)
    v = sp.Poly(vertical.as_expr().subs(assignment), X, domain=sp.QQ)
    inflated = source*v**k
    g = left_divide(inflated, u)
    if (g is None or g.degree() != m or g.gcd(g.diff()).degree() != 0 or g.compose(u) != inflated
            or not u.diff().rem(v**(k-1)).is_zero):
        raise ArithmeticError("An exported map fails its original identity or ramification check.")
    return {"inner": coefficients(u), "outer": coefficients(g), "vertical": coefficients(v)}


def solve_stratum(source, k, d, m, e, *, derivative=True, basis_solver=rational_points):
    variables, equations, inner, vertical, count = incidence(source, k, d, m, e, derivative=derivative)
    certificate = basis_solver(equations, variables) if variables else None
    return read_stratum(source, k, m, variables, equations, inner, vertical, count, certificate)


def read_stratum(source, k, m, variables, equations, inner, vertical, count, certificate):
    points = certificate["rational_points"] if variables else ([] if equations else [[]])
    maps = [exported_map(source, k, m, inner, vertical,
                         dict(zip(variables, map(sp.Rational, point), strict=True)))
            for point in points]
    return {"factor_parameters": count, "elimination_parameters": len(variables),
            "geometric_points": certificate["geometric_points"] if variables else len(points),
            "maps": sorted(maps, key=lambda row: (row["inner"], row["outer"], row["vertical"])),
            "certificate": certificate}


def search_plan(values, exponent, target_exponent=None, character=None, maximum_degree=None):
    if target_exponent is not None and (type(target_exponent) is not int or target_exponent < 2):
        raise ValueError("The target exponent must be an integer at least two.")
    if character is not None and (target_exponent is None or type(character) is not int
                                  or not 1 <= character < target_exponent or gcd(character, target_exponent) != 1):
        raise ValueError("The character must be a unit modulo the specified target exponent.")
    if maximum_degree is not None and (type(maximum_degree) is not int or maximum_degree < 2):
        raise ValueError("The degree cap must be an integer at least two.")
    removed, source, genus, minimum_pole = source_model(values, exponent)
    divisors = ([target_exponent] if target_exponent is not None else
                [k for k in range(2, exponent+1) if exponent % k == 0])
    plan = []
    for k in divisors:
        if exponent % k:
            continue
        index = exponent//k
        for j in ([character] if character is not None else [j for j in range(1, k) if gcd(j, k) == 1]):
            power_factor, branch = power_reduce(source, k, j)
            upper = maximum_degree//index if maximum_degree is not None else None
            if upper is not None and upper < 1:
                continue
            strata = list(eligible(branch, k, 1 if index > 1 else 2, upper))
            chart = {"target_exponent": k, "character": j, "working_branch": coefficients(branch),
                     "vertical_input_power": index*j,
                     "vertical_denominator": coefficients(removed**(index*j)*power_factor),
                     "vertical_input_degree_index": index}
            plan.append((chart, branch, strata))
    report = {"source_exponent": exponent, "source_genus": genus, "source": coefficients(source),
              "removed_power_factor": coefficients(removed), "minimum_noninvariant_pole": minimum_pole,
              "target_exponent_filter": target_exponent, "character_filter": character,
              "maximum_total_degree": maximum_degree,
              "automatic_equivariance_through": {str(k): (minimum_pole-1)//k for k in divisors}}
    return report, plan


def compute(values, exponent, target_exponent=None, character=None, maximum_degree=None):
    report, plan = search_plan(values, exponent, target_exponent, character, maximum_degree)
    report["charts"] = []
    for chart, branch, strata in plan:
        k, index = chart["target_exponent"], chart["vertical_input_degree_index"]
        report["charts"].append({**chart, "strata": [
            {"degrees": [d, m, e], "total_degree": index*d, **solve_stratum(branch, k, d, m, e)}
            for d, m, e in strata]})
    return report


def verify_report(values, exponent, report):
    expected, plan = search_plan(values, exponent, report["target_exponent_filter"],
                                 report["character_filter"], report["maximum_total_degree"])
    if {key: value for key, value in report.items() if key != "charts"} != expected:
        raise ArithmeticError("The certificate has incorrect source data or degree bounds.")
    for (chart, branch, strata), row in zip(plan, report["charts"], strict=True):
        if {key: value for key, value in row.items() if key != "strata"} != chart:
            raise ArithmeticError("A character chart is missing or altered.")
        k = chart["target_exponent"]
        for (d, m, e), entry in zip(strata, row["strata"], strict=True):
            if entry["degrees"] != [d, m, e] or entry["total_degree"] != chart["vertical_input_degree_index"]*d:
                raise ArithmeticError("A certified degree stratum is missing or altered.")
            variables, equations, inner, vertical, count = incidence(branch, k, d, m, e)
            certificate = entry["certificate"]
            if variables:
                if (certificate["variables"] != list(map(str, variables))
                        or certificate["equations"] != list(map(str, equations))):
                    raise ArithmeticError("The certificate uses a different incidence.")
                verify_finite(certificate)
            elif certificate is not None:
                raise ArithmeticError("A zero-variable stratum has an unexpected certificate.")
            expected_row = read_stratum(branch, k, m, variables, equations, inner, vertical, count, certificate)
            if {key: entry[key] for key in expected_row} != expected_row:
                raise ArithmeticError("An exported map list differs from the certified list.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON with exponent and ascending coefficients")
    parser.add_argument("--target-exponent", type=int)
    parser.add_argument("--character", type=int)
    parser.add_argument("--max-degree", type=int)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    if args.verify:
        verify_report(data["coefficients"], data["exponent"], json.loads(args.verify.read_text(encoding="utf-8")))
        print("Verified all cyclic quotient strata in the certificate.")
    else:
        result = compute(data["coefficients"], data["exponent"], args.target_exponent, args.character, args.max_degree)
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
