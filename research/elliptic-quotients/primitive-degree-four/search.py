"""Recover primitive quartic-degree elliptic maps from smooth plane quartics.

The search enumerates rational involution centers up to a supplied height
or recovers all rational centers by exact projective elimination.
Its construction uses pairs of commuting involutions, rational conics with
two rational points at infinity, and split cubic quotients related by a
rational affine change of x and a rational scale on y.
"""
import argparse
from itertools import combinations, permutations
import json
from pathlib import Path
import sys

import sympy as sp

sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "row-23640"))
from discover import discover, parse_quartic

X, Y, Z = sp.symbols("X Y Z")
U, V, W = sp.symbols("U V W")
s, t, r = sp.symbols("s t r")


def quartic_from_terms(data):
    terms = parse_quartic(data)
    return sum(c*X**i*Y**j*Z**k for (i, j, k), c in terms.items())


def smooth(F):
    variables = (X, Y, Z)
    gradient = [sp.diff(F, x) for x in variables]
    for x in variables:
        affine = [g.subs(x, 1) for g in gradient]
        basis = sp.groebner(affine, *(y for y in variables if x != y), domain=sp.QQ)
        if list(basis) != [1]:
            return False
    return True


def frames(symmetries):
    """Get three common eigenlines, then try the three choices of last axis."""
    emitted = set()
    for a, b in combinations(symmetries, 2):
        A, B = sp.Matrix(a["center"]), sp.Matrix(b["center"])
        la, lb = sp.Matrix(a["axis"]), sp.Matrix(b["axis"])
        if la.dot(B) or lb.dot(A):
            continue
        C = la.cross(lb)
        if C.is_zero_matrix:
            continue
        lines = [tuple(int(x/sp.gcd_list(list(L))) for x in L) for L in (A, B, C)]
        lines = [L if next(x for x in L if x) > 0 else tuple(-x for x in L) for L in lines]
        for last in lines:
            others = sorted(L for L in lines if L != last)
            if len(others) != 2:
                continue
            frame = (*others, last)
            if frame not in emitted:
                emitted.add(frame)
                yield sp.Matrix.hstack(*(sp.Matrix(L) for L in frame))


def square_conic(F, frame):
    xyz = frame*sp.Matrix([U, V, W])
    transformed = sp.Poly(F.subs(dict(zip((X, Y, Z), xyz)), simultaneous=True), U, V, W)
    if any(any(e % 2 for e in powers) for powers, _ in transformed.terms()):
        return None
    return sp.expand(sum(c*s**(i//2)*t**(j//2) for (i, j, _), c in transformed.terms()))


def conic_coordinates(conic):
    """Write s=U²/W², t=V²/W² as Laurent polynomials in a base parameter."""
    p = sp.Poly(conic, s, t)
    top = sum(c*s**i*t**j for (i, j), c in p.terms() if i+j == 2)
    scalar, factors = sp.factor_list(top, s, t)
    if len(factors) != 2 or any(e != 1 or sp.Poly(L, s, t).total_degree() != 1
                                for L, e in factors):
        return
    for (l1, _), (l2, _) in permutations(factors):
        if sp.det(sp.Matrix([[L.coeff(s), L.coeff(t)] for L in (l1, l2)])) == 0:
            continue
        # Choose constants so (l1+c1)(l2+c2) gives both linear terms.
        c1, c2 = sp.symbols("c1 c2")
        remainder = sp.Poly(conic/scalar - (l1+c1)*(l2+c2), s, t)
        shifts = sp.solve([remainder.coeff_monomial(s), remainder.coeff_monomial(t)],
                          (c1, c2), dict=True)[0]
        n = -remainder.coeff_monomial(1).subs(shifts)
        if n == 0:
            continue
        A, B = (l1+c1).subs(shifts), (l2+c2).subs(shifts)
        solved = sp.solve([A-r, B-n/r], (s, t), dict=True)[0]
        if sp.cancel(conic.subs(solved, simultaneous=True)) != 0:
            raise ArithmeticError("Conic parametrization failed")
        yield A, solved[s], solved[t]


def finite_roots(cubic):
    poly = sp.Poly(cubic, r, domain=sp.QQ)
    if poly.degree() != 3:
        return []
    roots = poly.ground_roots()
    return sorted(roots) if len(roots) == 3 else []


def identifications(first, second, roots1, roots2):
    """All rational isomorphisms of these cubic models preserving infinity."""
    for image in permutations(roots1):
        h = (image[1]-image[0]) / (roots2[1]-roots2[0])
        k = image[0] - h*roots2[0]
        if h*roots2[2]+k != image[2]:
            continue
        ratio = sp.cancel(first.subs(r, h*r+k)/second)
        if not ratio.is_Rational or ratio <= 0:
            continue
        L = sp.sqrt(ratio)
        if L.is_Rational:
            yield h, k, L


def squareclass_check(roots1, roots2, h, k):
    # Each bit records an odd valuation, including the pole at infinity.
    B1 = set(roots1) | {"infinity"}
    B2 = set(roots2) | {"infinity"}
    squareclasses = (set(), B1, B2, B1 ^ B2)
    evidence = []
    for e in roots1:
        preimage = (e-k)/h
        # Both pullback representatives are linear functions with a simple
        # pole at infinity. The two poles cancel in their quotient.
        odd = {e} ^ {preimage}
        if odd in squareclasses:
            return None
        evidence.append({"root": str(e), "odd_places": sorted(map(str, odd))})
    return evidence


def recover(data, height=2, *, all_rational_centers=False):
    F = quartic_from_terms(data)
    if all_rational_centers:
        sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "involution-atlas"))
        from classify import classify
        atlas = classify(data)
        found = [{key: tuple(map(sp.Rational, orbit[key])) for key in ("center", "axis")}
                 for orbit in atlas["orbits"] if orbit["orbit_size"] == 1]
        scope = {"center_search": "all rational centers", "rational_centers": len(found)}
    else:
        if not smooth(F):
            raise ValueError("The input plane quartic is singular")
        found = discover(data, height)
        scope = {"center_search": "coordinate height bound", "height": height}
    for frame in frames(found):
        conic = square_conic(F, frame)
        if conic is None:
            continue
        for base, u_squared, v_squared in conic_coordinates(conic):
            c1 = sp.cancel(r**2*v_squared)
            c2 = sp.cancel(r**2*u_squared)
            roots1, roots2 = finite_roots(c1), finite_roots(c2)
            if len(roots1) != 3 or len(roots2) != 3 or len(set(roots1) & set(roots2)) != 1:
                continue
            for h, k, L in identifications(c1, c2, roots1, roots2):
                proof = squareclass_check(roots1, roots2, h, k)
                if proof is None:
                    continue
                return {"status": "certified", "degree": 4, "primitive_over": "Qbar",
                        "defined_over": "Q", **scope,
                        "input_coordinates_X_Y_Z_from_U_V_W": [[int(x) for x in row] for row in frame.tolist()],
                        "base_r_in_s_t": str(base), "s": "(U/W)^2", "t": "(V/W)^2",
                        "U_squared_over_W_squared": str(sp.expand(u_squared)),
                        "V_squared_over_W_squared": str(sp.expand(v_squared)),
                        "target": {"left": "y^2", "right": str(sp.expand(c1)), "variable": "r"},
                        "first_point": ["r", "r*V/W"],
                        "second_point": [str(h*r+k), f"({L})*r*U/W"],
                        "map": "first_point + second_point in the target group law",
                        "base_cubic_1_roots": list(map(str, roots1)),
                        "base_cubic_2_roots": list(map(str, roots2)),
                        "isomorphism_h_k_L": list(map(str, (h, k, L))),
                        "nonzero_two_torsion": proof}
    return {"status": "unresolved", **scope, "centers_found": len(found)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quartic", type=Path, help="Quartic monomials [i,j,k,coefficient] as JSON")
    group = parser.add_mutually_exclusive_group()
    group.add_argument("--height", type=int, default=2)
    group.add_argument("--all-rational-centers", action="store_true",
                       help="Use exact elimination for the involution stage")
    args = parser.parse_args()
    try:
        report = recover(json.loads(args.quartic.read_text(encoding="utf-8")), args.height,
                         all_rational_centers=args.all_rational_centers)
    except (ValueError, OSError) as error:
        parser.exit(2, f"Input error: {error}\n")
    print(json.dumps(report, indent=2, default=str))


if __name__ == "__main__":
    main()
