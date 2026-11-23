"""Replay the degree-four discovery, verify the output formulas and test controls."""
import json
from pathlib import Path

import sympy as sp

from search import (X, Y, Z, U, V, W, r, s, t, recover, quartic_from_terms,
                    squareclass_check)


HERE = Path(__file__).resolve().parent


def require(test, why):
    if not test:
        raise ArithmeticError(why)


def expression(value, variable):
    # Parse only the verifier's internally generated arithmetic, using explicit
    # symbols. User input is parsed solely by the monomial reader in search.py.
    return sp.sympify(value, locals=variable)


def verify_chart(data, report):
    """Pull the two output points back to the input plane equation."""
    require(report["status"] == "certified", "No degree-four certificate produced")
    frame = sp.Matrix(report["input_coordinates_X_Y_Z_from_U_V_W"])
    require(frame.det() != 0, "Coordinate frame is singular")
    F = quartic_from_terms(data)
    xyz = frame*sp.Matrix([U, V, W])
    chart = sp.expand(F.subs(dict(zip((X, Y, Z), xyz)), simultaneous=True).subs(W, 1))
    R = expression(report["base_r_in_s_t"], {"s": s, "t": t}).subs({s: U**2, t: V**2})
    cubic = expression(report["target"]["right"], {"r": r})
    h, k, L = map(sp.Rational, report["isomorphism_h_k_L"])
    x1, y1 = R, R*V
    x2, y2 = h*R+k, L*R*U
    for x, y in ((x1, y1), (x2, y2)):
        require(sp.rem(sp.Poly(y**2-cubic.subs(r, x), V, U), sp.Poly(chart, V, U)).is_zero,
                "A recovered point fails the elliptic equation modulo the source")
    # Verify the third point of the secant on the cubic over its function
    # field. The roots of the line/cubic intersection are x1, x2 and x3.
    x0, xj, slope, intercept, x = sp.symbols("x0 xj slope intercept x")
    cubic_poly = sp.Poly(cubic.subs(r, x), x)
    c3, c2 = cubic_poly.nth(3), cubic_poly.nth(2)
    x3 = (slope**2-c2)/c3-x0-xj
    line_cut = cubic.subs(r, x) - (slope*x+intercept)**2
    remainder = sp.rem(line_cut, (x-x0)*(x-xj), x)
    require(sp.cancel(line_cut - remainder - c3*(x-x0)*(x-xj)*(x-x3)) == 0,
            "Secant addition polynomial identity failed")
    require(sp.rem(sp.Poly(x2-x1, V, U), sp.Poly(chart, V, U)).as_expr() != 0,
            "Secant denominator vanishes identically")

    # Rebuild all squareclasses independently via rational-function tests.
    cubic2 = sp.expand(cubic.subs(r, h*r+k)/L**2)
    f1 = sp.Poly(cubic, r)
    require(f1.degree() == 3 and sp.discriminant(f1.as_expr(), r) != 0,
            "Singular elliptic target")
    places = sp.polys.polytools.ground_roots(f1)
    require(len(places) == 3, "Expected three rational finite branch points")
    for root in places:
        ratio = (r-root)/(h*r+k-root)
        for twist in (1, cubic, cubic2, cubic*cubic2):
            n, d = sp.fraction(sp.cancel(ratio * twist))
            factors = sp.factor_list(n, r)[1] + sp.factor_list(d, r)[1]
            require(any(power % 2 for _, power in factors),
                    "Two-torsion character became a square on the source")


def fixed_formula():
    a = sp.symbols("a")
    source = (a*U**2+(a-1)*V**2+(2*a**2-2*a-1)*W**2)*(V**2-U**2+3*W**2) - (2*a-1)**2*W**4
    field = sp.QQ.frac_field(a, U)
    affine = sp.Poly(source.subs(W, 1), V, domain=field)
    R = (a*U**2+(a-1)*V**2+2*a**2-2*a-1)/(2*a-1)
    target = r*(r-1)*(r-a)
    for x, y in ((R, R*V), (R+1, R*U)):
        require(sp.rem(sp.Poly(y**2 - target.subs(r, x), V, domain=field), affine).is_zero,
                "Family map identity failed")
    normalized = quartic_from_terms(json.loads((HERE/"ciani.json").read_text()))
    expected = normalized.subs({X: U, Y: V, Z: W}, simultaneous=True)
    require(sp.expand(source.subs(a, 3)-expected) == 0, "Family does not specialize to the fixture")
    return source, a


def verify_ramification(family, a):
    """Verify the canonical line section and both finite branch values."""
    section = W**2*(3*(2*a-1)*U**2+2*(a+1)*(a-2)*W**2)
    require(sp.expand(family.subs(V, -U)-section) == 0, "Ramification line mismatch")
    u_squared = -2*(a+1)*(a-2)/(3*(2*a-1))
    base = (a*U**2+(a-1)*V**2+2*a**2-2*a-1)/(2*a-1)
    r_value = sp.factor(base.subs(V, -U).subs(U**2, u_squared))
    require(sp.cancel(r_value-(2*a-1)/3) == 0, "Wrong base coordinate at ramification")
    target_x = sp.factor(4*r_value**2*u_squared+a-2*r_value)
    target_y_over_u = sp.factor(2*r_value*(r_value-target_x)+r_value)
    require(sp.cancel(target_x+(a-2)*(4*a+1)**2/27) == 0, "Wrong finite branch x")
    require(sp.cancel(target_y_over_u-(2*a-1)**2*(4*a-5)*(4*a+1)/81) == 0,
            "Wrong finite branch y")
    require(sp.cancel(u_squared*target_y_over_u**2-target_x*(target_x-1)*(target_x-a)) == 0,
            "Branch values must lie on the elliptic target")
    return {"geometric_monodromy": "S4", "ramification_indices": [3, 2, 2],
            "a3_triple_branch_value": "O",
            "a3_simple_branch_x": str(target_x.subs(a, 3)),
            "a3_simple_branch_y_per_U": str(target_y_over_u.subs(a, 3)),
            "a3_simple_branch_U_squared": str(u_squared.subs(a, 3))}


def main():
    family, parameter = fixed_formula()
    results = {"family_ramification": verify_ramification(family, parameter)}
    for name, height in (("ciani", 1), ("sheared", 2)):
        data = json.loads((HERE/f"{name}.json").read_text(encoding="utf-8"))
        report = recover(data, height)
        verify_chart(data, report)
        results[name] = report
    # Height one already reconstructs the missing axis from the other two
    # centers on some shears. Use a unimodular change with all centers high.
    mixed = quartic_from_terms(json.loads((HERE/"ciani.json").read_text()))
    matrix = sp.Matrix([[5, -2, 0], [4, 1, -2], [-3, 0, 1]])
    require(abs(matrix.det()) == 1, "Test change must be unimodular")
    plane = sp.Matrix([X, Y, Z])
    mixed = sp.Poly(mixed.subs(dict(zip(plane, matrix*plane)), simultaneous=True), X, Y, Z)
    transformed = [[i, j, k, int(c)] for (i, j, k), c in mixed.terms()]
    require(recover(transformed, 1)["status"] == "unresolved", "Height-one control failed")
    complete = recover(transformed, all_rational_centers=True)
    verify_chart(transformed, complete)
    results["all_rational_centers_on_mixed_input"] = {
        "matrix": [[int(v) for v in row] for row in matrix.tolist()], "rational_centers": complete["rational_centers"],
        "height_one": "unresolved", "map_degree": complete["degree"]}
    for value in (-2, 4, 5):
        polynomial = sp.Poly(family.subs(parameter, value), U, V, W)
        data = [[i, j, k, int(c)] for (i, j, k), c in polynomial.terms()]
        verify_chart(data, recover(data, height=1))
    require(squareclass_check(list(map(sp.Integer, (0, 1, 3))), list(map(sp.Integer, (0, 4, 12))),
                              sp.Rational(1, 4), 0) is None,
            "The aligned common two-torsion control should be rejected")
    negative = json.loads((HERE.parent/"row-23640"/"quartic.json").read_text(encoding="utf-8"))
    report = recover(negative, height=2)
    require(report["status"] == "unresolved", "Unexpected primitive map on the excluded control")
    results["controls"] = {"other_family_parameters": [-2, 4, 5], "aligned_kernel_rejected": True,
                           "row_23640": report}
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
