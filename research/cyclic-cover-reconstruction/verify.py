"""Independent geometric controls and certificate replay for cyclic quotients."""

from copy import deepcopy
from math import gcd
import json
from pathlib import Path

import sympy as sp

from cyclic import X, compute, coefficients, polynomial, solve_stratum, source_model, verify_report
from fixtures import mixed, normalization_example, nonsaturated, picard, poly, saturated, known_map
from raw import triangular
from study import audit, output


HERE = Path(__file__).resolve().parent


def require(condition, reason):
    if not condition:
        raise ArithmeticError(reason)


def roundtrip(F, exponent, target=None, character=None, degree=None):
    values = coefficients(F)
    report = compute(values, exponent, target, character, degree)
    report = json.loads(json.dumps(report))
    verify_report(values, exponent, report)
    T = sp.Symbol("T")
    for chart in report["charts"]:
        k = chart["target_exponent"]
        power = chart["vertical_input_power"]
        denominator = polynomial(chart["vertical_denominator"])
        require(F**power == polynomial(chart["working_branch"])**(exponent//k)*denominator**exponent,
                "Character conversion failed the original equation.")
        for row in chart["strata"]:
            for item in row["maps"]:
                U, V, G = (polynomial(item[key]) for key in ("inner", "vertical", "outer"))
                y_power = sp.rem((T**power)**k, T**exponent-F.as_expr(), T)
                require(sp.expand(y_power*V.as_expr()**k/G.LC() -
                                  (denominator**k*G.compose(U)).as_expr()/G.LC()) == 0,
                        "A map fails in the original cyclic function field.")
    return report


def maps(report):
    return [(chart["target_exponent"], chart["character"], row["total_degree"], value)
            for chart in report["charts"] for row in chart["strata"] for value in row["maps"]]


def coefficient_checks():
    F = picard()
    expected = {2: known_map(poly(X**2-X), poly(1), poly(X**2-5*X/12+sp.Rational(5,18))),
                5: known_map(poly(X**5-5*X**4/2+5*X**3/3), poly(X*(X-1)), poly(X*(X-sp.Rational(1,6))))}
    positive = roundtrip(F, 3)
    require(maps(positive) == [(3, 1, d, item) for d, item in expected.items()],
            "Picard quartic lost or gained a map.")
    require(positive["automatic_equivariance_through"] == {"3": 1}, "Wrong Picard pole threshold.")
    negative = roundtrip(F+poly(X), 3)
    require(all(row["geometric_points"] == 0 for chart in negative["charts"] for row in chart["strata"]),
            "The perturbed Picard quartic acquired a geometric map.")
    for source in (F, F+poly(X)):
        for index, degrees in enumerate(((2, 2, 0), (5, 2, 2))):
            row = (positive if source == F else negative)["charts"][0]["strata"][index]
            require(triangular(source, 3, *degrees) == output(row), "Raw Picard elimination disagreed.")

    U, V, G = nonsaturated(3)
    source = normalization_example(3, U, G, V)
    report = roundtrip(source, 3, 3, 1)
    rows = report["charts"][0]["strata"]
    require([r["total_degree"] for r in rows] == [5,8,11]
            and [r["geometric_points"] for r in rows] == [0,1,0]
            and rows[1]["maps"] == [known_map(U,V,G)]
            and rows[1]["factor_parameters"] == 1, "The cubic complement search changed.")

    for k in (3, 4, 5, 7, 9, 11):
        U, V, G = saturated(k)
        source = normalization_example(k, U, G, V)
        row = solve_stratum(source, k, U.degree(), G.degree(), V.degree())
        require(row["maps"] == [known_map(U,V,G)] and row["factor_parameters"] == 0,
                "A larger saturated stratum disagreed.")
    for k in (3, 4, 5, 7):
        U, V, G = mixed(k)
        source = normalization_example(k, U, G, V)
        row = solve_stratum(source, k, U.degree(), G.degree(), V.degree())
        require(row["maps"] == [known_map(U,V,G)], "A mixed branch map was lost.")
        if k == 4:
            require(triangular(source, k, U.degree(), G.degree(), V.degree()) == output(row),
                    "Raw mixed-branch elimination disagreed.")
    square_target = poly((X**4+X+1)**2)
    row = solve_stratum(square_target,3,4,2,0)
    require(row["geometric_points"] == 0,
            "A power-composition whose target is not squarefree survived.")
    return positive


def extra_character_checks():
    F = picard()
    changed = F**2
    require([(k,j,d) for k,j,d,_ in maps(roundtrip(changed,3))] == [(3,2,2),(3,2,5)],
            "Nonidentity character branch disagreed.")
    # Degree 9 would split at infinity for exponent 6. Use degree five
    # with squarefree factors to get both quadratic and cubic subfields.
    F = poly((X+1)*(X**4+X+1))
    k6 = F*poly(X-3)**6
    report = roundtrip(k6,6,degree=3)
    found = [(k,j,d) for k,j,d,_ in maps(report)]
    require((2,1,3) in found and (3,1,2) in found and report["removed_power_factor"] == ["-3","1"],
            "A composite-exponent subfield or original-coordinate normalizer was lost.")
    partial = poly((X-1)**2*(X+2)**3)
    _, source, g, pole = source_model(coefficients(partial),6)
    require(source == partial and (g,pole) == (1,2), "Mixed inertia or pole-order calculation disagreed.")
    high = poly((X**2+X+1)**3*(X**5-X+1))
    require(gcd(high.degree(),2)==1, "Invalid valuation control.")
    _, source, _, min_pole = source_model(coefficients(high),2)
    require(source.degree() == 7 and min_pole == 7, "Power normalization altered the pole threshold.")
    u, g = poly(X**2+X), poly(X**4+X+1)
    bound = roundtrip(g.compose(u),3,3,degree=2)
    require(bound["source_genus"] == 7 and bound["automatic_equivariance_through"] == {"3":2}
            and maps(bound) == [(3,1,2,known_map(u,poly(1),g))],
            "The cover in the automatic-equivariance range changed.")
    require(compute(coefficients(F),6,4)["charts"] == [], "A nondividing cyclic exponent survived the equivariant filter.")


def malformed(report):
    bad = deepcopy(report)
    bad["charts"][0]["strata"].pop()
    try:
        verify_report(coefficients(picard()),3,bad)
    except (ValueError, ArithmeticError):
        pass
    else:
        raise ArithmeticError("A removed degree stratum passed verification.")
    substitutions = (
        (["0"], 3),
        (["1","0","0","1"], 3),
        (["1","1","1"], 1),
    )
    for coefficients_, exponent in substitutions:
        try:
            compute(coefficients_, exponent)
        except ValueError:
            pass
        else:
            raise ArithmeticError("Malformed or out-of-scope input passed validation.")


if __name__ == "__main__":
    report = coefficient_checks()
    extra_character_checks()
    malformed(report)
    study = audit()
    print(f"Verified cyclic quotients, five independent raw identities, and {study['observations']} benchmark observations.")
