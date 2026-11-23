"""All hyperelliptic maps with a single point over the chosen infinity."""
import argparse
import json
from pathlib import Path

import sympy as sp

from finite import coefficients, polynomial, rational_points, verify


X = sp.Symbol("x")


def curve(values):
    f = polynomial(values, X)
    if f.degree() < 3 or f.degree() % 2 != 1:
        raise ValueError("The hyperelliptic equation must have odd degree at least three.")
    if sp.gcd(f, f.diff()).degree() > 0:
        raise ValueError("The branch polynomial must be squarefree.")
    return f


def degrees(source, target):
    n, m = source.degree(), target.degree()
    bound = (n - 2) // (m - 2)
    return [(d, (m * d - n) // 2) for d in range(1, bound + 1)
            if (m * d - n) >= 0 and (m * d - n) % 2 == 0]


def incidence(source, target, degree):
    if degree not in dict(degrees(source, target)):
        raise ValueError("The requested degree is excluded by the infinity and genus bounds.")
    v_degree = dict(degrees(source, target))[degree]
    u = sp.symbols("u0:" + str(degree + 1))
    v = sp.symbols("v0:" + str(v_degree + 1))
    inverse = sp.Symbol("inverse")
    U = sum(a * X**i for i, a in enumerate(u))
    V = sum(a * X**i for i, a in enumerate(v))
    difference = sp.Poly(target.as_expr().subs(X, U) - source.as_expr() * V**2, X)
    equations = list(difference.all_coeffs()) + [inverse * u[-1] - 1]
    return equations, [*u, *v, inverse], U, V


def recover(source_values, target_values, prime=None, only_degree=None):
    source, target = curve(source_values), curve(target_values)
    allowed = degrees(source, target)
    if only_degree is not None:
        allowed = [(d, e) for d, e in allowed if d == only_degree]
        if not allowed:
            raise ValueError("The requested degree is excluded by the proved bounds.")
    result = {"source": coefficients(source), "target": coefficients(target),
              "source_genus": (source.degree() - 1) // 2,
              "target_genus": (target.degree() - 1) // 2,
              "degree_bound": (source.degree() - 2) // (target.degree() - 2),
              "scope": "all rational maps with singleton infinity fibre"
                       if only_degree is None else "specified degree with singleton infinity fibre",
              "degrees": []}
    for degree, v_degree in allowed:
        if v_degree == 0:
            from composition import recover_constant
            decomposition = recover_constant(source_values, target_values, prime)
            result["degrees"].append({"degree": degree, "v_degree": 0,
                                      "method": "constant-composition",
                                      "geometric_maps": decomposition["geometric_maps"],
                                      "rational_maps": decomposition["rational_maps"],
                                      "decomposition": decomposition})
            continue
        method = "finite-incidence"
        if source.degree() == 5 and target.degree() == 3:
            from cubic import cubic_incidence
            equations, variables, U, V = cubic_incidence(source, target)
            method = "cubic-classification"
        else:
            equations, variables, U, V = incidence(source, target, degree)
        certificate = rational_points(equations, variables, prime)
        maps = []
        for point in certificate["rational_points"]:
            substitution = dict(zip(variables, map(sp.Rational, point)))
            u = sp.Poly(U.subs(substitution), X)
            v = sp.Poly(V.subs(substitution), X)
            if sp.expand(target.as_expr().subs(X, u.as_expr()) - source.as_expr() * v.as_expr()**2) != 0:
                raise ArithmeticError("The recovered map failed the original curve equation.")
            maps.append({"u": coefficients(u), "v_over_y": coefficients(v), "degree": degree})
        result["degrees"].append({"degree": degree, "v_degree": v_degree, "method": method,
                                  "geometric_maps": certificate["geometric_points"],
                                  "rational_maps": maps, "certificate": certificate})
    return result


def verify_maps(result):
    source, target = curve(result["source"]), curve(result["target"])
    expected = dict(degrees(source, target))
    if result["scope"] not in {"all rational maps with singleton infinity fibre", "specified degree with singleton infinity fibre"}:
        raise ArithmeticError("Unknown map completeness scope.")
    if (result["source_genus"], result["target_genus"], result["degree_bound"]) != ((source.degree() - 1) // 2, (target.degree() - 1) // 2, (source.degree() - 2) // (target.degree() - 2)):
        raise ArithmeticError("Incorrect curve genera or degree bound.")
    if len(result["degrees"]) != len({row["degree"] for row in result["degrees"]}):
        raise ArithmeticError("The record repeats a map degree.")
    if result["scope"] == "all rational maps with singleton infinity fibre" and set(expected) != {r["degree"] for r in result["degrees"]}:
        raise ArithmeticError("The record omits an admissible degree.")
    for row in result["degrees"]:
        if row["degree"] not in expected or row["v_degree"] != expected[row["degree"]]:
            raise ArithmeticError("The record gives an inadmissible coefficient degree.")
        if row["method"] == "constant-composition":
            from composition import verify_constant
            decomposition = row["decomposition"]
            if row["v_degree"] != 0 or decomposition["source"] != result["source"] or decomposition["target"] != result["target"]:
                raise ArithmeticError("The decomposition uses different curves or degrees.")
            verify_constant(decomposition)
            if row["geometric_maps"] != decomposition["geometric_maps"] or row["rational_maps"] != decomposition["rational_maps"]:
                raise ArithmeticError("The exported map list differs from the decomposition.")
            continue
        if row["method"] == "cubic-classification":
            from cubic import cubic_incidence
            equations, variables, U, V = cubic_incidence(source, target)
        elif row["method"] == "finite-incidence":
            equations, variables, U, V = incidence(source, target, row["degree"])
        else:
            raise ArithmeticError("Unknown reconstruction method.")
        certificate = row["certificate"]
        if certificate["equations"] != list(map(str, equations)) or certificate["variables"] != list(map(str, variables)):
            raise ArithmeticError("The certificate is for different map equations.")
        verify(certificate)
        if row["geometric_maps"] != certificate["geometric_points"]:
            raise ArithmeticError("Wrong geometric map count.")
        if len(row["rational_maps"]) != len(certificate["rational_points"]):
            raise ArithmeticError("Wrong number of exported maps.")
        for point, mapping in zip(certificate["rational_points"], row["rational_maps"]):
            if mapping["degree"] != row["degree"]:
                raise ArithmeticError("Wrong degree on an exported map.")
            substitution = dict(zip(variables, map(sp.Rational, point)))
            if coefficients(sp.Poly(U.subs(substitution), X)) != mapping["u"]:
                raise ArithmeticError("Incorrect exported horizontal map.")
            if coefficients(sp.Poly(V.subs(substitution), X)) != mapping["v_over_y"]:
                raise ArithmeticError("Incorrect exported vertical map.")


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    parser.add_argument("--prime", type=int)
    parser.add_argument("--degree", type=int)
    parser.add_argument("--output", type=Path)
    parser.add_argument("--verify", type=Path)
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    if args.verify:
        stored = json.loads(args.verify.read_text())
        if coefficients(curve(data["source"])) != stored["source"] or coefficients(curve(data["target"])) != stored["target"]:
            raise ValueError("The certificate source or target differs from the input.")
        verify_maps(stored)
        print("PASS")
    else:
        result = recover(data["source"], data["target"], args.prime, args.degree)
        text = json.dumps(result, indent=2) + "\n"
        if args.output:
            args.output.write_text(text)
        else:
            print(text, end="")
