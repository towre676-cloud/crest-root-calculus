"""Reconstruct all four primitive covers at a=3 and check both allowable ranks."""
import json
from pathlib import Path
import sys

import sympy as sp

from reconstruct import Q, S, T, reconstruct


def normalized_line(line):
    scale = next(c for c in line if c != 0)
    return tuple(sp.simplify(sp.sympify(c)/scale) for c in line)


def involution_axes(source):
    sys.path.insert(0, str(source.parents[1] / "involution-atlas"))
    from classify import classify
    parameter = sp.Symbol("T")
    axes = set()
    for orbit in classify(json.loads(source.read_text()))["orbits"]:
        polynomial = sp.Poly(orbit["field_polynomial"], parameter)
        for value in polynomial.all_roots():
            line = [sp.sympify(c).subs(parameter,value) for c in orbit["axis"]]
            axes.add(normalized_line(line))
    return axes


def main():
    source = Path(__file__).resolve().parent.parent / "primitive-degree-four/ciani.json"
    raw = json.loads(source.read_text())["terms"]
    axes = involution_axes(source)
    primitive = set()
    rows = {}
    for name, origin, twisted, u, v, a, b, expected in [
        ("sum", 1, False, 0, -5, "-112/3", "-1280/27", 2),
        ("difference", -1, False, 0, -5, "-112/3", "-1280/27", 2),
        ("quadratic_plus", 1, True, 10, -10, "-7/12", "5/54", 2),
        ("quadratic_minus", 1, True, -10, -10, "-7/12", "5/54", 2),
        ("isogeny_after_double", 1, False, 0, -10, "-52/3", "736/27", 3),
        ("double", 1, False, 0, -10, "-7/3", "-20/27", 6),
        ("changed_target", 1, False, 0, -5, "-112/3", "-1253/27", 0)
    ]:
        chart = sum(c*(origin+S)**i*T**k * (Q(-2)**(-(k//2)) if twisted else 1)
                    for i,j,k,c in raw).expand()
        result = reconstruct(chart, u, v, a, b)
        if result["nullity"] != expected:
            raise ArithmeticError(f"{name}: expected nullity {expected}, got {result}")
        if result["degree_four"]:
            slope = chart.diff(S).subs({T:0,S:0})
            t_scale = sp.sqrt(-2) if twisted else 1
            line = normalized_line([v,slope-v*origin,u*t_scale])
            result["canonical_line"] = list(map(str,line))
            result["primitive"] = line not in axes
            if result["primitive"]:
                primitive.add(line)
            if result["primitive"] != (name != "isogeny_after_double"):
                raise ArithmeticError(f"{name}: wrong involution-axis classification")
        rows[name] = {"point":[origin,1,0], "quadratic_chart":twisted, "parameters":[u,v,a,b], **result}
    if len(axes) != 5 or len(primitive) != 4:
        raise ArithmeticError("Failed to separate the four classes from the five involution axes")
    print(json.dumps(rows, indent=2))


if __name__ == "__main__":
    main()
