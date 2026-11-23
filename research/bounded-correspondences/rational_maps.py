"""All maps between hyperelliptic curves with target genus >= 2.

Two normalized coefficient charts in every Riemann--Hurwitz permitted degree
cover the complete finite map scheme. Exact elimination can be expensive.
"""
import argparse
import json
from pathlib import Path

import sympy as sp

from finite import coefficients, polynomial, rational_points, verify
from maps import X


def curve(values):
    polynomial_ = polynomial(values, X)
    if polynomial_.degree() < 5:
        raise ValueError("The unrestricted map theorem requires hyperelliptic genus at least two.")
    if sp.gcd(polynomial_, polynomial_.diff()).degree() > 0:
        raise ValueError("The branch polynomial must be squarefree.")
    return polynomial_


def chart(source, target, degree, leading_chart):
    h = (target.degree() - 1) // 2
    g = (source.degree() - 1) // 2
    if h < 2 or degree < 1 or degree > (g - 1) // (h - 1):
        raise ValueError("The genus or degree does not satisfy the complete-map bounds.")
    c_degree = ((2 * h + 2) * degree - source.degree()) // 2
    if c_degree < 0:
        return [sp.Integer(1)], [sp.Symbol("inverse")], sp.Integer(0), sp.Integer(1), sp.Integer(0)
    a = list(sp.symbols("a0:" + str(degree)))
    if leading_chart == 1:
        A = X**degree + sum(value * X**i for i, value in enumerate(a))
        b = list(sp.symbols("b0:" + str(degree + 1)))
        B = sum(value * X**i for i, value in enumerate(b))
    elif leading_chart == 2:
        A = sum(value * X**i for i, value in enumerate(a))
        b = list(sp.symbols("b0:" + str(degree)))
        B = X**degree + sum(value * X**i for i, value in enumerate(b))
    else:
        raise ValueError("There are exactly two leading-coefficient charts.")
    if g - 1 == degree * (h - 1):
        c = [sp.Symbol("differential_scale")]
        C = c[0] * (sp.diff(A, X) * B - A * sp.diff(B, X))
    else:
        c = list(sp.symbols("c0:" + str(c_degree + 1)))
        C = sum(value * X**i for i, value in enumerate(c))
    homogeneous = sum(target.nth(i) * A**i * B**(target.degree() - i)
                      for i in range(target.degree() + 1))
    difference = sp.Poly(B**(2 * h + 2 - target.degree()) * homogeneous - source.as_expr() * C**2, X)
    resultant = sp.resultant(A, B, X)
    inverse = sp.Symbol("inverse")
    equations = list(difference.all_coeffs()) + [inverse * resultant - 1]
    return equations, [*a, *b, *c, inverse], A, B, C


def recover(source_values, target_values, prime=None, only_degree=None, only_chart=None):
    source, target = curve(source_values), curve(target_values)
    g = (source.degree() - 1) // 2
    h = (target.degree() - 1) // 2
    if h < 2:
        raise ValueError("The unrestricted map theorem requires target genus at least two.")
    maximum_degree = (g - 1) // (h - 1)
    requested_degrees = range(1, maximum_degree + 1) if only_degree is None else [only_degree]
    requested_charts = [1, 2] if only_chart is None else [only_chart]
    result = {"source": coefficients(source), "target": coefficients(target),
              "source_genus": g, "target_genus": h, "degree_bound": maximum_degree,
              "scope": "all nonconstant rational morphisms" if only_degree is None and only_chart is None
                       else "specified degree or coefficient chart", "charts": []}
    for degree in requested_degrees:
        for leading_chart in requested_charts:
            equations, variables, A, B, C = chart(source, target, degree, leading_chart)
            certificate = rational_points(equations, variables, prime)
            maps = []
            for point in certificate["rational_points"]:
                assignment = dict(zip(variables, map(sp.Rational, point)))
                polys = [sp.Poly(value.subs(assignment), X) for value in [A, B, C]]
                a, b, c = polys
                homogeneous = sum(target.nth(i) * a.as_expr()**i * b.as_expr()**(target.degree() - i)
                                  for i in range(target.degree() + 1))
                if sp.expand(b.as_expr()**(2 * h + 2 - target.degree()) * homogeneous - source.as_expr() * c.as_expr()**2) != 0:
                    raise ArithmeticError("A recovered rational map failed its original equation.")
                maps.append({"A": coefficients(a), "B": coefficients(b), "C": coefficients(c),
                             "degree": degree, "vertical_denominator_exponent": h + 1})
            result["charts"].append({"degree": degree, "leading_chart": leading_chart,
                                     "geometric_maps": certificate["geometric_points"],
                                     "rational_maps": maps, "certificate": certificate})
    return result


def verify_maps(result):
    source, target = curve(result["source"]), curve(result["target"])
    g, h = (source.degree() - 1) // 2, (target.degree() - 1) // 2
    if h < 2:
        raise ArithmeticError("The certificate requires target genus at least two.")
    maximum_degree = (g - 1) // (h - 1)
    if (g, h, maximum_degree) != (result["source_genus"], result["target_genus"], result["degree_bound"]):
        raise ArithmeticError("Incorrect genus or degree bound.")
    if result["scope"] not in {"all nonconstant rational morphisms", "specified degree or coefficient chart"}:
        raise ArithmeticError("Unknown map completeness scope.")
    expected = {(degree, leading_chart) for degree in range(1, maximum_degree + 1)
                for leading_chart in [1, 2]}
    actual = [(row["degree"], row["leading_chart"]) for row in result["charts"]]
    if len(actual) != len(set(actual)):
        raise ArithmeticError("A map coefficient chart is repeated.")
    if result["scope"] == "all nonconstant rational morphisms" and set(actual) != expected:
        raise ArithmeticError("The record omits an admissible map chart.")
    for row in result["charts"]:
        equations, variables, A, B, C = chart(source, target, row["degree"], row["leading_chart"])
        certificate = row["certificate"]
        if certificate["equations"] != list(map(str, equations)) or certificate["variables"] != list(map(str, variables)):
            raise ArithmeticError("The certificate encodes a different map chart.")
        verify(certificate)
        expected_maps = []
        for point in certificate["rational_points"]:
            assignment = dict(zip(variables, map(sp.Rational, point)))
            values = [coefficients(sp.Poly(value.subs(assignment), X)) for value in [A, B, C]]
            expected_maps.append({"A": values[0], "B": values[1], "C": values[2],
                                  "degree": row["degree"], "vertical_denominator_exponent": h + 1})
        if expected_maps != row["rational_maps"] or row["geometric_maps"] != certificate["geometric_points"]:
            raise ArithmeticError("The exported map list differs from the certified points.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--degree", type=int)
    parser.add_argument("--chart", type=int, choices=[1, 2])
    parser.add_argument("--prime", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    parser.add_argument("--screen", action="store_true", help="Try exact Frobenius divisibility before elimination.")
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    if args.screen and not args.verify:
        from frobenius_filter import screen
        local = screen(data["source"], data["target"])
        if local["status"] == "rational_maps_excluded":
            text = json.dumps(local, indent=2) + "\n"
            if args.output:
                args.output.write_text(text)
            else:
                print(text, end="")
            raise SystemExit(0)
    if args.verify:
        saved = json.loads(args.verify.read_text())
        if coefficients(curve(data["source"])) != saved["source"] or coefficients(curve(data["target"])) != saved["target"]:
            raise ValueError("The certificate uses different curve equations.")
        verify_maps(saved)
        print("PASS")
    else:
        result = recover(data["source"], data["target"], args.prime, args.degree, args.chart)
        text = json.dumps(result, indent=2) + "\n"
        if args.output:
            args.output.write_text(text)
        else:
            print(text, end="")
