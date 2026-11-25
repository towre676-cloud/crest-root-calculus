"""Exact cyclic quotient search with finite horizontal poles on the branch locus."""

import argparse
from dataclasses import dataclass
from itertools import combinations, product
import json
from math import gcd
from pathlib import Path

import sympy as sp
from sympy.polys.ring_series import rs_mul, rs_nth_root, rs_pow
from sympy.polys.rings import ring

from cyclic import (X, basic_equations, coefficients, polynomial, rational_points,
                    search_plan, twice_genus, verify_finite)


@dataclass(frozen=True)
class Poles:
    horizontal: sp.Poly
    vertical: sp.Poly
    removed: sp.Poly
    derivative_factor: sp.Poly


def pole_sets(components, k, m, degree_bound, fixed):
    """Every finite pole divisor supported on the normalized branch locus."""
    factors = [(p.monic(), r, (-r*pow(m, -1, k)) % k)
               for component, r in components for p, _ in component.factor_list()[1]]
    one = sp.Poly(1, X, domain=sp.QQ)
    if fixed is not None:
        remaining = fixed
        orders = []
        for p, _, residue in factors:
            order = 0
            while remaining.rem(p).is_zero:
                remaining = remaining.exquo(p)
                order += 1
            if order and order % k != residue:
                return
            orders.append(order)
        if remaining.degree() == 0:
            yield build_poles(factors, orders, k, m, one)
        return
    for size in range(1, len(factors)+1):
        for support in combinations(factors, size):
            least_q = sum(p.degree()*residue for p, _, residue in support)
            if least_q >= degree_bound:
                continue
            orders = [range(residue, 1+(degree_bound-1-least_q)//p.degree()+residue, k)
                      for p, _, residue in support]
            for choices in product(*orders):
                if sum(p.degree()*b for (p, _, _), b in zip(support, choices)) < degree_bound:
                    yield build_poles(support, choices, k, m, one)


def build_poles(support, orders, k, m, one):
    horizontal = vertical = removed = derivative_factor = one
    for (prime, residue, _), order in zip(support, orders, strict=True):
        if order:
            horizontal *= prime**order
            vertical *= prime**((m*order+residue)//k)
            removed *= prime**residue
            derivative_factor *= prime**(order-1)
    return Poles(horizontal, vertical, removed, derivative_factor)


def eligible(branch, k, maximum, fixed=None):
    """Yield a pole divisor and every admissible degree triple over it."""
    _, components = branch.sqf_list()
    genus = twice_genus(k, components)
    support_degree = sum(p.degree() for p, _ in components)
    n = branch.degree()
    for m in range(2, 2+genus//(k-1)):
        if gcd(k, m) != 1:
            continue
        slope = (k-1)*m-k
        # Since each chosen branch place contributes at least one to deg removed,
        # (k-1)*(n-deg removed)+k*(chosen support degree-1) is at most this.
        universal = ((k-1)*n+support_degree-k)//slope
        if maximum is not None:
            universal = min(universal, maximum)
        genus_term = (k-1)*(m-1)-2
        if genus_term > 0:
            universal = min(universal, (genus-2)//genus_term)
        if universal < 2:
            continue
        for poles in pole_sets(components, k, m, universal, fixed):
            q = poles.horizontal.degree()
            support = q-poles.derivative_factor.degree()
            weight = n-poles.removed.degree()
            bound = min(universal, ((k-1)*weight+k*(support-1))//slope)
            for d in range(q+1, bound+1):
                e, remainder = divmod(m*d-weight, k)
                if (not remainder and 0 <= (k-1)*e <= d+support-1
                        and weight >= m*(d % k)):
                    yield poles, d, m, e


def reverse(poly, length, series, t):
    degree = poly.degree()
    return sum((poly.rep.nth(degree-j)*t**j for j in range(min(length, degree+1))), series.zero)


def wronskian(inner, denominator):
    return inner.diff()*denominator-inner*denominator.diff()


def recover_numerator(wanted, denominator, d, domain):
    """Invert the linear Wronskian map on the normalized numerator chart."""
    q = denominator.degree()
    inner = sp.Poly.from_dict({(d,): domain.one}, X, domain=domain)
    remaining = wanted-wronskian(inner, denominator)
    for power in range(d-1, -1, -1):
        if power == q:
            continue
        coeff = domain.exquo(remaining.rep.nth(power+q-1), domain.convert(power-q))
        if coeff:
            term = sp.Poly.from_dict({(power,): coeff}, X, domain=domain)
            inner += term
            remaining -= wronskian(term, denominator)
    return inner, remaining


def complementary(f, denominator, derivative_factor, k, d, m, e, names, domain):
    q = denominator.degree()
    gap = d-q
    w = d+q-1-derivative_factor.degree()-(k-1)*e
    known = min(e, gap-1)
    series, t = ring("t", domain)
    fv = reverse(f, known+1, series, t)
    fv *= domain.exquo(domain.one, f.rep.LC())
    den = reverse(denominator, known+1, series, t)
    fixed = reverse(derivative_factor, known+1, series, t)
    residual = series.one+sum((c*t**i for i, c in enumerate(names[:w], 1)), series.zero)
    vertical = series.one
    for j in range(1, known+1):
        root = rs_nth_root(rs_mul(fv, rs_pow(vertical, k, t, j+1), t, j+1), m, t, j+1)
        lhs = gap*root*den-t*root.diff(t)*den+t*root*den.diff(t)
        rhs = gap*rs_mul(rs_pow(vertical, k-1, t, j+1), fixed*residual, t, j+1)
        divisor = m*gap*(k-1)-k*(gap-j)
        coeff = domain.exquo(m*(lhs-rhs).get((j,), domain.zero), domain.convert(divisor))
        vertical += coeff*t**j
    vertical += sum((c*t**j for j, c in enumerate(names[w:], known+1)), series.zero)
    v = sp.Poly.from_list([vertical.get((j,), domain.zero) for j in range(e+1)], X, domain=domain)
    factor = sp.Poly.from_list([domain.one, *names[:w]], X, domain=domain)
    inner, remainder = recover_numerator(gap*derivative_factor*v**(k-1)*factor, denominator, d, domain)
    return inner, v, remainder


def direct(f, denominator, derivative_factor, k, d, m, e, names, domain):
    vertical = sp.Poly.from_list([domain.one, *names[:e]], X, domain=domain)
    inflated = f*vertical**k
    gap = d-denominator.degree()
    series, t = ring("t", domain)
    reversed_f = reverse(inflated, gap, series, t)
    reversed_f *= domain.exquo(domain.one, inflated.rep.LC())
    root = rs_nth_root(reversed_f, m, t, gap)
    inner = sp.Poly.from_list([*[root.get((j,), domain.zero) for j in range(gap)],
                              domain.zero, *names[e:]], X, domain=domain)
    remainder = wronskian(inner, denominator).rem(derivative_factor*vertical**(k-1))
    return inner, vertical, remainder


def homogenized_outer(inflated, inner, denominator, m):
    """Extract the descending coefficients of G_h(A,Q) and its residual."""
    gap = inner.degree()-denominator.degree()
    top = m*denominator.degree()
    values = []
    for i in range(m, -1, -1):
        coeff = inflated.rep.nth(top+i*gap)
        values.append(coeff)
        if coeff:
            inflated -= (inner**i*denominator**(m-i)).mul_ground(coeff)
    return sp.Poly.from_list(values, X, domain=inner.domain), inflated


def incidence(branch, k, poles, d, m, e):
    q = poles.horizontal.degree()
    derivative_factor = poles.derivative_factor
    w = d+q-1-derivative_factor.degree()-(k-1)*e
    tail = max(0, e-(d-q)+1)
    use_complement = w+tail < e+q
    count = min(w+tail, e+q)
    params = sp.symbols(f"b0:{count}")
    domain = sp.QQ.poly_ring(*params) if params else sp.QQ
    names = domain.gens if params else ()
    f = branch.exquo(poles.removed).set_domain(domain)
    den = poles.horizontal.set_domain(domain)
    build = complementary if use_complement else direct
    inner, vertical, residual = build(f, den, derivative_factor.set_domain(domain), k, d, m, e, names, domain)
    equations = basic_equations(residual.coeffs())
    if equations != [1]:
        outer, residual = homogenized_outer(f*vertical**k, inner, den, m)
        equations = basic_equations([*equations, *residual.coeffs()])
        if equations != [1]:
            unit = inner.resultant(den)*outer.discriminant()
            if not unit:
                equations = [sp.Integer(1)]
            elif not unit.is_number:
                inverse = sp.Symbol("open_inverse")
                equations.append(inverse*unit-1)
                params = (*params, inverse)
    return params, equations, inner, vertical, count, use_complement


def export(branch, k, poles, m, inner, vertical, assignment):
    den = poles.horizontal
    a = sp.Poly(inner.as_expr().subs(assignment), X, domain=sp.QQ)
    v = sp.Poly(vertical.as_expr().subs(assignment), X, domain=sp.QQ)
    f = branch.exquo(poles.removed)
    outer, remainder = homogenized_outer(f*v**k, a, den, m)
    if (remainder or poles.vertical**k != den**m*poles.removed
            or a.gcd(den).degree() != 0 or not outer.discriminant()
            or a.degree() <= den.degree() or a.nth(den.degree()) or v.LC() != 1):
        raise ArithmeticError("A map violates the exact identity or its normalization.")
    return {"numerator": coefficients(a), "outer": coefficients(outer),
            "vertical": coefficients(v), "vertical_pole_denominator": coefficients(poles.vertical)}


def certified_stratum(branch, k, poles, d, m, e, *, check=False, supplied=None):
    variables, equations, a, v, count, complement = incidence(branch, k, poles, d, m, e)
    certificate = supplied if check else rational_points(equations, variables) if variables else None
    if check and variables:
        if (certificate is None or certificate["variables"] != list(map(str, variables))
                or certificate["equations"] != list(map(str, equations))):
            raise ArithmeticError("The finite certificate uses different equations.")
        verify_finite(certificate)
    elif not variables and certificate is not None:
        raise ArithmeticError("Unexpected certificate for a zero-variable identity.")
    points = certificate["rational_points"] if variables else ([] if equations else [[]])
    maps = [export(branch, k, poles, m, a, v, dict(zip(variables, map(sp.Rational, row), strict=True)))
            for row in points]
    return {"denominator": coefficients(poles.horizontal), "degrees": [d, m, e], "factor_parameters": count,
            "uses_complement": complement, "elimination_parameters": len(variables),
            "geometric_points": certificate["geometric_points"] if variables else len(maps),
            "maps": maps, "certificate": certificate}


def plan(values, exponent, target, character, degree, finite_poles):
    header, charts = search_plan(values, exponent, target, character, degree)
    # This solver requires equivariance and has a different infinity condition.
    del header["automatic_equivariance_through"]
    fixed = polynomial(finite_poles) if finite_poles is not None else None
    if fixed is not None:
        if not fixed or fixed.degree() < 1:
            raise ValueError("A denominator filter must have positive degree.")
        fixed = fixed.monic()
    header["denominator_filter"] = coefficients(fixed) if fixed is not None else None
    planned = []
    for chart, branch, _ in charts:
        index = chart["vertical_input_degree_index"]
        bound = degree//index if degree is not None else None
        strata = list(eligible(branch, chart["target_exponent"], bound, fixed))
        planned.append((chart, branch, strata))
    return header, planned


def compute(values, exponent, target=None, character=None, degree=None, finite_poles=None):
    header, charts = plan(values, exponent, target, character, degree, finite_poles)
    return {**header, "charts": [{**chart, "strata": [
        {"total_degree": chart["vertical_input_degree_index"]*d,
         **certified_stratum(branch, chart["target_exponent"], poles, d, m, e)}
        for poles, d, m, e in strata]} for chart, branch, strata in charts]}


def verify(values, exponent, report):
    header, charts = plan(values, exponent, report["target_exponent_filter"],
                          report["character_filter"], report["maximum_total_degree"], report["denominator_filter"])
    if {key: value for key, value in report.items() if key != "charts"} != header:
        raise ArithmeticError("Source or search scope differs from the supplied certificate.")
    for (chart, branch, strata), row in zip(charts, report["charts"], strict=True):
        if {key: value for key, value in row.items() if key != "strata"} != chart:
            raise ArithmeticError("A character chart is missing or altered.")
        for stratum, entry in zip(strata, row["strata"], strict=True):
            actual = certified_stratum(branch, chart["target_exponent"], *stratum,
                                       check=True, supplied=entry["certificate"])
            actual["total_degree"] = chart["vertical_input_degree_index"]*stratum[1]
            if actual != entry:
                raise ArithmeticError("The certified stratum or map list has changed.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path, help="JSON with exponent and ascending coefficients")
    parser.add_argument("--target-exponent", type=int)
    parser.add_argument("--character", type=int)
    parser.add_argument("--max-degree", type=int)
    parser.add_argument("--denominator", nargs="+", help="monic denominator filter, ascending coefficients")
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    values, exponent = data["coefficients"], data["exponent"]
    if args.verify:
        verify(values, exponent, json.loads(args.verify.read_text(encoding="utf-8")))
        print("Verified all branch-supported pole strata in the certificate.")
    else:
        result = compute(values, exponent, args.target_exponent, args.character, args.max_degree, args.denominator)
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
