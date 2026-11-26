"""Bounded unknown-target quotient search with arbitrary horizontal poles."""

import argparse
from dataclasses import dataclass
import json
from itertools import chain, combinations
from math import gcd
from pathlib import Path

import sympy as sp

from cyclic import (X, basic_equations, coefficients, polynomial, rational_points,
                    search_plan, twice_genus, verify_finite)
from rational import (Poles, build_poles, complementary, direct, homogenized_outer,
                      pole_sets)


@dataclass(frozen=True)
class FibreStratum:
    poles: Poles
    unbranched_degree: int
    numerator_degree: int
    degree: int
    target_degree: int
    multiplier_degree: int
    infinity: str
    profile: tuple[tuple[int, int], ...]
    first: int

    def denominator_degree(self, k):
        return self.poles.horizontal.degree() + k*self.unbranched_degree


def target_profiles(k, source_genus_twice, include_mixed):
    """Enumerate finite branch multiplicities with bounded target genus."""
    steps = [(r, k-gcd(k, r)) for r in range(1, k if include_mixed else 2)]
    budget = source_genus_twice + k-1

    def populate(i, left, profile, m):
        if i == len(steps):
            genus_twice = 1-k+budget-left
            if genus_twice >= 2 and gcd(m, k) == 1:
                yield profile, m, genus_twice
            return
        order, weight = steps[i]
        for n in range(left//weight, -1, -1):
            next_profile = (*profile, (order, n)) if n else profile
            yield from populate(i+1, left-weight*n, next_profile, m+order*n)

    yield from populate(0, budget, (), 0)


def strata(branch, k, index, cap, infinity, higher_genus_only, include_mixed):
    """Enumerate degree data and the rational branch-supported part of Q."""
    _, components = branch.sqf_list()
    source_genus_twice = twice_genus(k, components)
    one = polynomial([1])
    for profile, m, target_genus_twice in target_profiles(k, source_genus_twice, include_mixed):
        if higher_genus_only and target_genus_twice == 2:
            continue
        bound = cap//index
        if target_genus_twice > 2:
            bound = min(bound, (source_genus_twice-2)//(target_genus_twice-2))
        minimum = 2 if index == 1 else 1
        if bound < minimum:
            continue
        P = (k-1)//profile[-1][0]
        branch_poles = chain((build_poles((), (), k, m, one),),
                             pole_sets(components, k, m, bound+1, None))
        for poles in branch_poles:
            q0 = poles.horizontal.degree()
            n = branch.degree()-poles.removed.degree()
            for h in range((bound-q0)//k+1):
                q = q0+k*h
                fixed_derivative_degree = poles.derivative_factor.degree()+(k-1)*h
                if infinity in ("both", "pole"):
                    for d in range(max(minimum, q+1), bound+1):
                        e, rem = divmod(m*d-n, k)
                        if (not rem and 0 <= P*e <= d+q-1-fixed_derivative_degree
                                and n >= sum(count*(nu*d % k) for nu, count in profile)):
                            yield FibreStratum(poles, h, d, d, m, e, "pole", profile, m)
                if infinity in ("both", "root") and minimum <= q <= bound:
                    residue = sum(count*(nu*q % k) for nu, count in profile)
                    for first, _ in profile:
                        if gcd(first, k) != 1:
                            continue
                        for a in range(q):
                            e, rem = divmod((m-first)*q+first*a-n, k)
                            if (not rem and 0 <= P*e <= a+q-1-fixed_derivative_degree
                                    and n >= residue-(first*q % k)+(first*a % k)):
                                yield FibreStratum(poles, h, a, q, m, e, "root", profile, first)


def monic(names, domain):
    return sp.Poly.from_list([domain.one, *names], X, domain=domain)


def finite_root_outer(inflated, a, q, m, first):
    """Extract coefficients of G, where G(0)=0 and deg(a)<deg(q)."""
    gap = q.degree()-a.degree()
    values = [inflated.domain.zero]*first
    for i in range(first, m+1):
        coeff = inflated.rep.nth(m*q.degree()-i*gap)
        values.append(coeff)
        if coeff:
            inflated -= (a**i*q**(m-i)).mul_ground(coeff)
    return sp.Poly.from_list(values[::-1], X, domain=inflated.domain), inflated


def target_polynomial(inflated, a, q, s):
    if s.infinity == "pole":
        return homogenized_outer(inflated, a, q, s.target_degree)
    return finite_root_outer(inflated, a, q, s.target_degree, s.first)


def profile_polynomials(s, names, domain):
    """Encode each group of roots of a given multiplicity by a monic factor."""
    groups = []
    for order, degree in s.profile:
        if order == 1:
            continue
        at_zero = s.infinity == "root" and order == s.first
        size = degree-int(at_zero)
        factor = monic(names[:size], domain)
        if at_zero:
            factor *= sp.Poly(X, X, domain=domain)
        groups.append((order, factor))
        names = names[size:]
    return groups


def impose_open(equations, parameters, branch, s, inner, hidden, outer, groups, domain):
    if equations == [1]:
        return equations, parameters
    basis = sp.groebner(equations, *parameters, domain=sp.QQ, order="grevlex") if equations else None
    equations = list(basis) if basis is not None else equations
    if equations == [1]:
        return equations, parameters

    def reduced(poly):
        if basis is None:
            return poly
        return sp.Poly.from_list([basis.reduce(c)[1] for c in poly.all_coeffs()], X, domain=domain)

    a0, h0, g0 = map(reduced, (inner, hidden, outer))
    units = [h0.resultant(branch.set_domain(domain)), a0.resultant(h0),
             a0.resultant(s.poles.horizontal.set_domain(domain)), g0.nth(s.target_degree),
             g0.nth(s.first)]
    if groups:
        factors = [reduced(poly) for poly in groups]
        units.extend(poly.discriminant() for poly in factors)
        units.extend(a.resultant(b) for a, b in combinations(factors, 2))
    else:
        units.append(g0.discriminant())
    if any(unit == 0 for unit in units):
        return [sp.Integer(1)], parameters
    unit = sp.sympify(sp.prod(unit for unit in units if not unit.is_number))
    if basis is not None:
        unit = basis.reduce(unit)[1]
    if not unit:
        return [sp.Integer(1)], parameters
    if not unit.is_number:
        inverse = sp.Symbol("open_inverse")
        return [*equations, inverse*unit-1], (*parameters, inverse)
    return equations, parameters


def incidence(branch, k, s):
    h, a, m, e = (s.unbranched_degree, s.numerator_degree,
                   s.target_degree, s.multiplier_degree)
    P = (k-1)//s.profile[-1][0]
    q_degree = s.denominator_degree(k)
    gap = a-q_degree
    w = a+q_degree-1-s.poles.derivative_factor.degree()-(k-1)*h-P*e
    known = min(e, abs(gap)-1)
    direct_free = e+(q_degree if gap > 0 else max(0, a+1+gap))
    use_complement = w+e-known < direct_free
    free = w+e-known if use_complement else direct_free
    group_free = sum(count for order, count in s.profile if order > 1)-int(gap < 0 and s.first > 1)
    parameters = sp.symbols(f"b0:{h+free+group_free}")
    domain = sp.QQ.poly_ring(*parameters) if parameters else sp.QQ
    names = domain.gens if parameters else ()
    hidden = monic(names[:h], domain)
    f = branch.exquo(s.poles.removed).set_domain(domain)
    q = s.poles.horizontal.set_domain(domain)*hidden**k
    derivative_factor = s.poles.derivative_factor.set_domain(domain)*hidden**(k-1)
    builder = complementary if use_complement else direct
    inner, v, residual = builder(f, q, derivative_factor, k, a, m, e, names[h:h+free], domain,
                                 first=s.first, derivative_power=P)
    constraints = residual.coeffs()
    if gap < 0 and 0 < s.first*a < k:
        # A multiplier zero in the fibre over zero would force
        # the local order of G(A/Q) to be at least k.
        constraints.extend(f.rem(inner**s.first).coeffs())
    equations = basic_equations(constraints)
    if equations != [1]:
        outer, residual = target_polynomial(f*v**k, inner, q, s)
        equations = basic_equations([*equations, *residual.coeffs()])
        groups = profile_polynomials(s, names[h+free:], domain)
        factors = [poly for _, poly in groups]
        if groups and equations != [1]:
            repeated = sp.prod(poly**order for order, poly in groups)
            simple, residual = outer.div(repeated)
            equations = basic_equations([*equations, *residual.coeffs()])
            if s.profile[0][0] == 1:
                factors.append(simple)
        equations, parameters = impose_open(equations, parameters, branch, s, inner, hidden, outer,
                                             factors, domain)
    return parameters, equations, (inner, hidden, v), h+free, group_free, use_complement


def exported(branch, k, s, polynomials, assignment):
    a, hidden, v = [sp.Poly(p.as_expr().subs(assignment), X, domain=sp.QQ)
                    for p in polynomials]
    q = s.poles.horizontal*hidden**k
    d = s.poles.vertical*hidden**s.target_degree
    f = branch.exquo(s.poles.removed)
    outer, remainder = target_polynomial(f*v**k, a, q, s)
    actual_profile = tuple((order, poly.degree()) for poly, order in outer.sqf_list()[1])
    if (remainder or a.gcd(q).degree() != 0 or hidden.gcd(branch).degree() != 0
            or outer.degree() != s.target_degree or actual_profile != s.profile
            or not outer.nth(s.first) or d**k != q**s.target_degree*s.poles.removed):
        raise ArithmeticError("An exported map violates the homogeneous identity.")
    return {"numerator": coefficients(a), "denominator": coefficients(q),
            "unbranched_pole_factor": coefficients(hidden),
            "outer": coefficients(outer), "vertical": coefficients(v),
            "vertical_pole_denominator": coefficients(d)}


def certified_stratum(branch, k, index, s, *, supplied=None, check=False):
    variables, equations, polynomials, parameter_count, group_count, complement = incidence(branch, k, s)
    certificate = supplied if check else rational_points(equations, variables) if variables else None
    if check:
        if variables:
            if (certificate is None or certificate["variables"] != list(map(str, variables))
                    or certificate["equations"] != list(map(str, equations))):
                raise ArithmeticError("The finite certificate uses different equations.")
            verify_finite(certificate)
        elif certificate is not None:
            raise ArithmeticError("Unexpected certificate in a zero-variable stratum.")
    points = certificate["rational_points"] if variables else ([] if equations else [[]])
    maps = [exported(branch, k, s, polynomials, dict(zip(variables, map(sp.Rational, row), strict=True)))
            for row in points]
    return {
        "degree": index*s.degree, "numerator_degree": s.numerator_degree,
        "target_degree": s.target_degree, "multiplier_degree": s.multiplier_degree,
        "target_profile": [list(entry) for entry in s.profile],
        "target_root_multiplicity": s.first if s.infinity == "root" else None,
        "infinity_image": s.infinity, "branch_pole_denominator": coefficients(s.poles.horizontal),
        "unbranched_degree": s.unbranched_degree, "factor_parameters": parameter_count,
        "target_branch_parameters": group_count, "elimination_parameters": len(variables),
        "uses_complement": complement,
        "geometric_points": certificate["geometric_points"] if variables else len(maps),
        "maps": maps, "certificate": certificate,
    }


def plan(values, exponent, degree, target, character, infinity, higher_genus_only=False, include_mixed=False):
    if degree is None and not higher_genus_only:
        raise ValueError("Supply a degree bound for the unrestricted elliptic fibres.")
    if infinity not in ("both", "pole", "root"):
        raise ValueError("The image of source infinity must be a pole or a finite branch point.")
    header, charts = search_plan(values, exponent, target, character, degree)
    del header["automatic_equivariance_through"]
    header["infinity_image_filter"] = infinity
    header["higher_genus_only"] = higher_genus_only
    header["mixed_target_multiplicities"] = include_mixed
    cap = max(0, header["source_genus"]-1) if degree is None else degree
    return header, [(chart, branch, list(strata(branch, chart["target_exponent"],
                                              chart["vertical_input_degree_index"], cap,
                                              infinity, higher_genus_only, include_mixed)))
                    for chart, branch, _ in charts]


def compute(values, exponent, degree=None, target=None, character=None, infinity="both", *,
            higher_genus_only=False, include_mixed=False):
    header, charts = plan(values, exponent, degree, target, character, infinity,
                          higher_genus_only, include_mixed)
    return {**header, "charts": [
        {**chart, "strata": [certified_stratum(branch, chart["target_exponent"],
                                             chart["vertical_input_degree_index"], s) for s in entries]}
        for chart, branch, entries in charts
    ]}


def verify(values, exponent, report):
    header, charts = plan(values, exponent, report["maximum_total_degree"],
                          report["target_exponent_filter"], report["character_filter"],
                          report["infinity_image_filter"], report["higher_genus_only"],
                          report["mixed_target_multiplicities"])
    if {key: value for key, value in report.items() if key != "charts"} != header:
        raise ArithmeticError("Source or search scope differs from the supplied certificate.")
    for (chart, branch, entries), row in zip(charts, report["charts"], strict=True):
        if {key: value for key, value in row.items() if key != "strata"} != chart:
            raise ArithmeticError("A character chart differs from the supplied certificate.")
        for s, supplied in zip(entries, row["strata"], strict=True):
            actual = certified_stratum(branch, chart["target_exponent"], chart["vertical_input_degree_index"],
                                       s, supplied=supplied["certificate"], check=True)
            if actual != supplied:
                raise ArithmeticError("A degree stratum or a map differs from the supplied certificate.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--max-degree", type=int)
    parser.add_argument("--target-exponent", type=int)
    parser.add_argument("--higher-genus-only", action="store_true",
                        help="Exclude elliptic targets and derive a degree bound from Riemann–Hurwitz.")
    parser.add_argument("--character", type=int)
    parser.add_argument("--mixed-targets", action="store_true",
                        help="Include all normalized target branch multiplicity profiles.")
    parser.add_argument("--infinity-image", choices=("both", "pole", "root"), default="both")
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text(encoding="utf-8"))
    values, exponent = data["coefficients"], data["exponent"]
    if args.verify:
        verify(values, exponent, json.loads(args.verify.read_text(encoding="utf-8")))
        print("Verified all bounded fibre strata in the certificate.")
    elif args.max_degree is None and not args.higher_genus_only:
        parser.error("supply --max-degree or use --higher-genus-only")
    else:
        result = compute(values, exponent, args.max_degree, args.target_exponent,
                         args.character, args.infinity_image, higher_genus_only=args.higher_genus_only,
                         include_mixed=args.mixed_targets)
        print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
