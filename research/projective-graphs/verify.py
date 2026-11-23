"""Check the projective quotient graphs, projection degrees and boundary fibres."""
import json
from pathlib import Path

import sympy as sp

from graph import graph_basis, hilbert_polynomial, read_sparse


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def boundary_fibre(graph, source, target, point, k):
    t, u, v = source
    fixes = dict(zip(target, point))
    equations = [g.subs(fixes) for g in graph]
    finite = sp.groebner([g.subs(u, 1) for g in equations], t, v, k, domain=sp.QQ)
    # On the complementary line u=0, either t=1 or v=1.
    for coordinate, other in ((t, v), (v, t)):
        boundary = [g.subs({u: 0, coordinate: 1}) for g in equations]
        require(list(sp.groebner(boundary, other, k, domain=sp.QQ)) == [1],
                "Found an unexpected point outside the u=1 chart.")
    return list(finite)


def main():
    t, u, v, x, y, z, k = sp.symbols("t u v X Y Z k")
    source, target = (t, u, v), (x, y, z)
    a, b = u**2+7*u*v+2*v**2, u*v*(u**2+3*u*v+v**2)
    r = t**2 + 23*k*u*v
    equation = 3*r**2 + a*r + b
    w, g = 6*r+a, u**2+u*v+8*v**2
    h = 2*(w+g)
    forms = ((h-15*v**2)*v, h*(2*u+v), v**3)
    image_equation = x*(x-z)*(x+15*z)-y**2*z

    graph, specializes = graph_basis(source, target, [equation], forms, (k,))
    require(specializes and len(graph) == 7, "The universal graph basis failed its specialization check.")
    require(sp.gcd_list(forms) == 1, "Map coordinates have an unexpected common ambient factor.")
    eliminated = [f for f in graph if f.free_symbols.isdisjoint(source+(k,))]
    require(len(eliminated) == 1 and sp.expand(eliminated[0]-image_equation) == 0,
            "The target is not the expected Weierstrass cubic.")

    # Check both projections scheme-theoretically at the two points missing
    # from the dense v != 0 formula.
    at_infinity = boundary_fibre(graph, source, target, (0, 1, 0), k)
    at_torsion = boundary_fibre(graph, source, target, (-15, 0, 1), k)
    require(at_infinity == [t**2, v], "Unexpected fibre over the elliptic origin.")
    require(at_torsion == [t**2+sp.Rational(1, 3), v], "Unexpected fibre over (-15,0).")

    over_zero = boundary_fibre(graph, source, target, (0, 0, 1), k)
    over_one = boundary_fibre(graph, source, target, (1, 0, 1), k)
    require(over_zero == [t**2-46*k-sp.Rational(2, 3), v+2],
            "Unexpected fibre over (0,0).")
    require(over_one == [t**2-46*k-1, v+2], "Unexpected fibre over (1,0).")
    # In characteristic zero the diagonal meets the inverse graph on E
    # precisely in E[2], reduced of degree 4. Each pulled-back fibre has
    # degree two, verified by its monic quadratic and linear equation.
    fibres = [at_infinity, at_torsion, over_zero, over_one]
    intersection_length = sum(sp.degree(f[0], t) for f in fibres)
    require(intersection_length == 8, "Graph/negative-graph intersection mismatch.")

    # At either base point the three raw forms vanish. At target (0:0:1)
    # the naive minors vanish too, while a saturated graph equation excludes it.
    false_point = {u: 1, v: 0, x: 0, y: 0, z: 1, k: 0}
    modulus = t**2 + sp.Rational(1, 3)
    require(all(sp.rem(f.subs(false_point), modulus, t) == 0 for f in forms),
            "The base-point witness is invalid.")
    require(any(sp.rem(f.subs(false_point), modulus, t) != 0 for f in graph),
            "The saturated graph retained a false target point.")

    fixed = [f.subs(k, 0) for f in graph]
    hp = hilbert_polynomial(fixed, source, target)
    m, n = sp.symbols("m n")
    require(hp == 4*m+6*n-2, "Incorrect graph Hilbert polynomial.")
    full_fixed, _ = graph_basis(source, target, [equation.subs(k, 0)], [f.subs(k, 0) for f in forms])
    require(fixed == full_fixed, "The direct row does not match the universal specialization.")

    # Verify the runnable sparse input against the mathematical formula.
    data = json.loads(Path(__file__).with_name("row-23640.json").read_text(encoding="utf-8"))
    provided = data["equations"] + data["map"]
    expressions = [equation, *forms]
    require(data["source_coordinates"] == len(source) and len(provided) == len(expressions)
            and all(sp.expand(read_sparse(p, source)-e.subs(k, 0)) == 0
                    for p, e in zip(provided, expressions)), "Sparse fixture differs from row 23640.")

    proof = {
        "graph_basis_over_Q_k_polynomial_ring": [str(f) for f in graph],
        "monic_basis_specializes_for_every_k": specializes,
        "bigraded_hilbert_polynomial": str(hp),
        "projection_degrees_source_target": [1, 2],
        "fibre_at_origin_in_chart_u_1": [str(f) for f in at_infinity],
        "fibre_at_minus15_0_in_chart_u_1": [str(f) for f in at_torsion],
        "fibre_at_0_0_in_chart_u_1": [str(f) for f in over_zero],
        "fibre_at_1_0_in_chart_u_1": [str(f) for f in over_one],
        "graph_intersection_with_negative_graph": int(intersection_length),
        "spurious_naive_graph_point_removed": True,
    }
    print(json.dumps(proof, indent=2))


if __name__ == "__main__":
    main()
