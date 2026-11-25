"""Check rational double quotients against an independent branch-point enumeration."""

from copy import deepcopy
from itertools import permutations

import sympy as sp

from double_covers import X, coefficients, function_record, recover, verify


def check(condition, reason):
    if not condition:
        raise ArithmeticError(reason)


def projective_point(x):
    return (1, 0) if x is None else (x, 1)


def keys(records):
    return {str((entry["involution"]["x"], entry["involution"]["y_multiplier"]))
            for entry in records}


def independent_decks(source):
    """A split branch set fixes a projective map by any three images."""
    roots = list(source.ground_roots())
    check(len(roots) == source.degree(), "The independent test needs a completely split source.")
    if source.degree() % 2:
        roots.append(None)
    branches = set(roots)
    domain = [projective_point(x) for x in roots[:3]]
    h = len(roots) // 2
    answer = set()
    for positions in permutations(roots, 3):
        matrix = sp.Matrix([[x0*y1, x1*y1, -x0*y0, -x1*y0]
                            for (x0, x1), (y0, y1) in zip(domain, map(projective_point, positions))])
        kernel = matrix.nullspace()
        check(len(kernel) == 1, "Three distinct source and target points must fix a unique map.")
        a, b, c, d = kernel[0]
        if a + d != 0 or a*d == b*c:
            continue
        images = set()
        for x0, x1 in map(projective_point, roots):
            denominator = c*x0+d*x1
            images.add((a*x0+b*x1)/denominator if denominator else None)
        if images != branches:
            continue
        m = (a*X+b)/(c*X+d)
        ratio = sp.cancel(source.as_expr().subs(X, m)*(c*X+d)**(2*h)/source.as_expr())
        check(bool(ratio.is_Rational), "The branch symmetry failed to give a rational constant.")
        root = sp.sqrt(ratio)
        if not root.is_Rational:
            continue
        scale = root/(c*X+d)**h
        if sp.cancel(scale*scale.subs(X, m) - 1) == 0:
            for sign in (1, -1):
                answer.add(str((function_record(m), function_record(sign*scale))))
    return answer


def main():
    sextic = (X**2-1)*(X**2-4)*(X**2-9)
    odd = sp.cancel(X**6*sextic.subs(X, 1+1/X))
    octic = sextic*(X**2-16)
    twist = sp.Rational(-7, 5)
    cases = [
        ("split-even", sextic, [1, 1], True),
        ("split-odd-moved-infinity", odd, [1, 1], True),
        ("split-genus-three", octic, [1, 2], True),
        ("odd-involutions-without-rational-lifts", X*(X**2-1)*(X**2-4), [], True),
        ("reciprocal", (X**2+1)*(X**2+2)*(2*X**2+1), [1]*4, False),
        ("nonsquare-lift", (X**2+2)*(X**4+X**2+4), [1, 1], False),
        ("nonrational-fixed-points", (X**4+X**2+4)*(X**4+3*X**2+4), [1, 1, 1, 2, 2, 2], False),
        ("generic-odd", X**5-X+1, [], False),
    ]
    for name, expression, genera, split in cases:
        source = sp.Poly(expression, X, domain=sp.QQ)
        values = coefficients(source)
        result = recover(values)
        verify(values, result)
        maps = result["quotients"]
        check(sorted(row["target_genus"] for row in maps) == genera, f"{name} quotient genera differ.")
        check(len(keys(maps)) == len(maps), f"{name} contains duplicate deck transformations.")
        if split:
            check(independent_decks(source) == keys(maps), f"{name} disagrees with branch-point enumeration.")
        if name == "nonsquare-lift":
            check(["0", "2"] in result["nonaffine_certificate"]["rational_points"],
                  "The projective symmetry should exist even though its involution is not rational.")
    original = sp.Poly(twist*sextic*(X**2+X+1)**2, X, domain=sp.QQ)
    values = coefficients(original)
    result = recover(values)
    verify(values, result)
    maps = result["quotients"]
    check(result["square"] == ["1", "1", "1"] and result["source"][-1] == "-7/5" and len(maps) == 2,
          "The inflated input lost a square factor, its rational unit or one of its double quotients.")
    target_infinity = next(row for row in maps if len(row["target"]) == 4)
    check(target_infinity["map"]["u"] == function_record(-X**2),
          "The two-infinity test lost its two poles over the elliptic model's unique infinity.")
    tampered = deepcopy(result)
    tampered["quotients"].pop()
    try:
        verify(values, tampered)
    except ArithmeticError:
        pass
    else:
        raise ArithmeticError("A removed quotient passed the coverage check.")
    print("Verified nine double-quotient inputs, four independent split branch searches, and a tampered certificate.")


if __name__ == "__main__":
    main()
