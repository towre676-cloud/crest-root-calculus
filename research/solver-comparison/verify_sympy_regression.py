"""Replay the exact SymPy 1.14 decomposition failure and its elementary witnesses.

The stock library is never changed. --observe reports the currently installed
version; default replay requires the archived source fingerprints and results.
"""
import argparse
import hashlib
import inspect
import json
from pathlib import Path

import sympy as sp
import sympy.polys.densetools as dense


HERE = Path(__file__).resolve().parent
X = sp.Symbol("x")


def ascending(poly):
    poly = sp.Poly(poly, X, domain=sp.QQ)
    return [str(poly.nth(i)) for i in range(poly.degree() + 1)]


def coefficient_recurrence(poly, inner_degree):
    """Independent coefficient form of the recurrence proved in the note."""
    n = poly.degree()
    outer_degree = n // inner_degree
    inner = {inner_degree: sp.S.One}
    for i in range(1, inner_degree):
        value = sum((i - (outer_degree + 1)*j) * poly.nth(n-i+j)
                    * inner[inner_degree-j] for j in range(i))
        inner[inner_degree-i] = value / (i*outer_degree*poly.LC())
    return sp.Poly(sum(value*X**degree for degree,value in inner.items()), X, domain=sp.QQ)


def observation(poly, inner_degree, exact_inner, exact_outer):
    poly = sp.Poly(poly, X, domain=sp.QQ)
    if sp.gcd(poly, poly.diff()).degree() != 0:
        raise ArithmeticError("The regression polynomial is not squarefree.")
    if sp.expand(exact_outer.subs(X, exact_inner) - poly.as_expr()) != 0:
        raise ArithmeticError("The exact composition witness failed.")
    corrected = coefficient_recurrence(poly.monic(), inner_degree)
    normalized_inner = sp.Poly(exact_inner, X).monic()
    normalized_inner -= sp.Poly(normalized_inner.nth(0), X)
    if corrected != normalized_inner.set_domain(sp.QQ):
        raise ArithmeticError("The independent coefficient recurrence failed.")
    values = poly.monic().rep.to_list()
    right = dense._dup_right_decompose(values, inner_degree, sp.QQ)
    left = dense._dup_left_decompose(values, right, sp.QQ)
    public = sp.decompose(poly.as_expr(), domain=sp.QQ)
    return {
        "source": ascending(poly),
        "squarefree": True,
        "exact_inner": ascending(exact_inner),
        "exact_outer": ascending(exact_outer),
        "stock_public_decomposition": [ascending(item) for item in public],
        "stock_fixed_degree_inner": ascending(sp.Poly.from_list(right, X)),
        "stock_fixed_degree_outer": None if left is None else ascending(sp.Poly.from_list(left, X)),
        "corrected_monic_zero_constant_inner": ascending(corrected),
        "public_missed_nontrivial_decomposition": len(public) == 1,
    }


def evidence():
    small_inner = X**3 + 2*X**2
    small = observation(sp.expand(small_inner**2+1), 3, small_inner, X**2+1)
    cases = json.loads((HERE/"map-panel.json").read_text())["cases"]
    case = next(row for row in cases if row["id"] == "g2-d3-twist")
    source = sum(sp.Rational(value)*X**i for i,value in enumerate(case["source"]))
    target = sum(sp.Rational(value)*X**i for i,value in enumerate(case["target"]))
    U = X**3 - 2*X**2 + 2
    scale_squared = sp.Rational(9,8)
    if sp.expand(target.subs(X,U) - scale_squared*source) != 0:
        raise ArithmeticError("The benchmark cover identity failed.")
    map_case = observation(source, 3, U, sp.Rational(8,9)*target)
    # Exhaust the residual affine target change independently of the map solver.
    a = sp.Symbol("a")
    H = X**3 - 2*X**2
    K = X**5 + 6*X**4 + 8*X**3 - 16*X**2 - 52*X - 44
    if sp.expand(K.subs(X,H) - sp.Rational(9,8)*source) != 0:
        raise ArithmeticError("The monic decomposition failed.")
    shift = (6*a+4)/5
    residual = sp.Poly(target.subs(X,a*X+shift)-a**5*K, X)
    equations = [sp.Poly(c,a,domain=sp.QQ) for c in residual.all_coeffs() if c != 0]
    eliminant = equations[0]
    for equation in equations[1:]:
        eliminant = sp.gcd(eliminant,equation)
    while eliminant.nth(0) == 0:
        eliminant = eliminant.exquo(sp.Poly(a,a))
    if eliminant.monic().as_expr() != a-1:
        raise ArithmeticError("The affine target automorphism calculation changed.")
    if sp.sqrt(scale_squared).is_Rational:
        raise ArithmeticError("The recorded twist became rational.")
    map_case.update({"id": case["id"], "target": case["target"],
                     "vertical_scale_squared": str(scale_squared),
                     "affine_scale_eliminant": str(eliminant.monic().as_expr()),
                     "geometric_maps": 2, "rational_maps": 0})
    module_path = Path(inspect.getsourcefile(dense))
    return {"schema":1, "sympy_version":sp.__version__,
            "module":"sympy/polys/densetools.py",
            "installed_module_sha256":hashlib.sha256(module_path.read_bytes()).hexdigest(),
            "installed_right_helper_sha256":hashlib.sha256(inspect.getsource(dense._dup_right_decompose).encode()).hexdigest(),
            "minimal_counterexample":small, "benchmark_counterexample":map_case}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--observe", action="store_true",
                        help="Print evidence for the currently installed library without requiring version 1.14.")
    args = parser.parse_args()
    result = evidence()
    if args.observe:
        print(json.dumps(result,indent=2))
        return
    archived = json.loads((HERE/"sympy-regression.json").read_text())
    if result != archived:
        raise ArithmeticError("The installed implementation or its behavior differs from the archived SymPy 1.14 evidence. Use --observe to inspect the new version.")
    print("PASS: archived stock SymPy decomposition behavior, exact identities and corrected recurrence.")


if __name__ == "__main__":
    main()
