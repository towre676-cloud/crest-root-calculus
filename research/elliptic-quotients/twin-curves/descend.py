"""Construct quartic-field models and verify the maps to the archived sextics."""

import argparse
import json
from pathlib import Path

import sympy as sp


HERE = Path(__file__).resolve().parent
X = sp.Symbol("x")
F_POLY = "x^4-x^3-3*x^2+2*x+4"
L_POLY = "x^24+x^18+27*x^12+163*x^6+64"
H_MINUS_171 = (
    1, 694282057876537344, 472103267541360574464,
    8391550371275812148084736, -1311901521779155773721411584,
)


def poly(expression):
    return sp.Poly(expression, X, domain=sp.QQ)


def column(a, degree):
    return sp.Matrix([a.nth(j) for j in range(degree)])


def powers(a, modulus):
    result = [poly(1)]
    for _ in range(modulus.degree() - 1):
        result.append((a * result[-1]).rem(modulus))
    return result


def multiplication_polynomial(a, modulus):
    n = modulus.degree()
    columns = [column((a * poly(X**j)).rem(modulus), n) for j in range(n)]
    return poly(sp.Matrix.hstack(*columns).charpoly(X).as_expr())


def from_dihedral(u, v, modulus):
    """Return the coefficients of the descended sextic from constant upward."""
    delta2 = (u**2 - 4*v**3).rem(modulus)
    if not u or not v or not delta2:
        raise ArithmeticError("The archived pair must satisfy uv(u^2-4v^3) != 0.")
    a = u - v**2
    b = 3*delta2 + 3*u**2 - 2*u*v**2
    c = delta2*(15*u + v**2)
    d = 2*delta2*(5*delta2 + 5*u**2 + 2*u*v**2)
    coefficients = [delta2**3*a, delta2**2*b, delta2*c, d, c, b, a]
    return [e.rem(modulus) for e in coefficients]


def check_change(coefficients, s1, s2, embed, L):
    """Check the integral coefficient identity for the geometric isomorphism."""
    v = (s1*s2).rem(L)
    s, t = [(x**3).rem(L) for x in (s1, s2)]
    delta = s - t

    def times(a, b):
        product = [poly(0) for _ in range(len(a) + len(b) - 1)]
        for i, x in enumerate(a):
            for j, y in enumerate(b):
                product[i+j] = (product[i+j] + x*y).rem(L)
        return product

    P, Q = [delta, poly(1)], [-delta, poly(1)]
    P_power, Q_power = [[poly(1)]], [[poly(1)]]
    for _ in range(6):
        P_power.append(times(P_power[-1], P))
        Q_power.append(times(Q_power[-1], Q))
    powers_and_factors = [(6, s**2), (4, -s*v**2), (2, t*v**2), (0, -t**2)]
    rhs = [poly(0) for _ in range(7)]
    for k, factor in powers_and_factors:
        for j, value in enumerate(times(P_power[k], Q_power[6-k])):
            rhs[j] = (rhs[j] + factor*value).rem(L)
    if [(delta*embed(c)).rem(L) for c in coefficients] != rhs:
        raise ArithmeticError("The degree-24 change of variables failed.")


def elliptic_j_pair(u, v, disc, F, class_polynomial=H_MINUS_171):
    """Verify the unordered quotient pair using its quadratic over F."""
    inv = 256 * disc.invert(F)
    trace = (inv * (u**2 - 2*v**3 - 9*u*v + 54*v**2 - 27*u)).rem(F)
    product = (inv**2 * (v**2 - 3*u + 9*v)**3).rem(F)
    remainder = (poly(0), poly(0))
    for coeff in class_polynomial:
        a, b = remainder
        remainder = ((poly(coeff)-b*product).rem(F), (a+b*trace).rem(F))
    if any(remainder) or not (trace**2-4*product).rem(F):
        raise ArithmeticError("The two distinct quotient invariants must be roots of the ring class polynomial.")
    return [str(a.as_expr()) for a in (trace, product)]


def compute():
    F, L = poly(F_POLY), poly(L_POLY)
    if not F.is_irreducible or not L.is_irreducible:
        raise ArithmeticError("The supplied field polynomials must be irreducible.")
    dihedral = json.loads((HERE / "dihedral_-19.json").read_text())
    reduced = json.loads((HERE / "equations_-19_reduced.json").read_text())
    moduli = json.loads((HERE / "moduli_-19.json").read_text())
    rows = {}
    absolute_polynomials = set()
    common_j_pair = None
    for index, (name, data) in enumerate(dihedral.items()):
        u, v = [poly(data[key]["expr_in_theta"]).rem(F) for key in ("u", "v")]
        v_powers = powers(v, F)
        coords = sp.Matrix.hstack(*(column(z, 4) for z in v_powers))
        if coords.det() == 0:
            raise ArithmeticError("The v-invariant does not generate the quartic field.")
        v_poly = multiplication_polynomial(v, F)
        if v_poly != poly(reduced["v_minpolys"][index]):
            raise ArithmeticError("The v minimal polynomial disagrees with the archive.")

        s1, s2 = [poly(reduced[name][key]).rem(L) for key in ("s1", "s2")]
        v_in_L = (s1*s2).rem(L)
        theta_as_polynomial_of_v = coords.inv() * column(poly(X), 4)
        theta_in_L = sum((a*b for a, b in zip(theta_as_polynomial_of_v, powers(v_in_L, L))),
                         poly(0)).rem(L)

        def embed(a):
            return a.compose(theta_in_L).rem(L)

        if embed(F) or embed(v) != v_in_L or embed(u) != (s1**3 + s2**3).rem(L):
            raise ArithmeticError("Could not embed the archived dihedral data in the original field.")
        smoothness = (v**2 - 4*u + 18*v - 27).rem(F)
        if not smoothness:
            raise ArithmeticError("The original bielliptic model is singular.")
        j_pair = elliptic_j_pair(u, v, smoothness, F)
        if common_j_pair is not None and j_pair != common_j_pair:
            raise ArithmeticError("The two curves must have the same unordered quotient pair.")
        common_j_pair = j_pair
        # On x^6-s1*x^4+s2*x^2-1, the Clebsch-normalized I2 is
        # 16(v+15) and I10 is 64 times the square of the cubic discriminant.
        # The corresponding absolute ratio is eight times the stored PARI convention.
        disc = (64*smoothness**2).rem(F)
        abs_ratio = ((16*v+240)**5 * disc.invert(F)).rem(F)
        recorded_ratio = poly(moduli[name]["J2^5/J10"]["expr_in_theta"])
        if abs_ratio != (8*recorded_ratio).rem(F):
            raise ArithmeticError("The independent absolute invariant differs from the archive.")
        abs_poly = multiplication_polynomial(abs_ratio, F)
        if not abs_poly.is_irreducible or abs_poly in absolute_polynomials:
            raise ArithmeticError("The moduli ratio must generate F and distinguish the orbits.")
        absolute_polynomials.add(abs_poly)
        coefficients = from_dihedral(u, v, F)
        check_change(coefficients, s1, s2, embed, L)
        rows[name] = {
            "coefficients_constant_first": [str(z.as_expr()) for z in coefficients],
            "theta_in_eta": str(theta_in_L.as_expr()),
            "cubic_discriminant": str(smoothness.as_expr()),
            "v_minimal_polynomial": str(v_poly.as_expr()),
            "clebsch_I2_power_5_over_I10_minimal_polynomial": str(abs_poly.as_expr()),
        }
    return {"F": F_POLY, "L": L_POLY, "H_minus_171_descending": list(H_MINUS_171),
            "elliptic_j_sum_and_product_in_F": common_j_pair, "curves": rows}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", action="store_true", help="Check the stored exact models.")
    args = parser.parse_args()
    result = compute()
    if args.verify and result != json.loads((HERE / "descent_-19.json").read_text()):
        raise ArithmeticError("The exact descent differs from its receipt.")
    print(json.dumps(result, indent=2))


if __name__ == "__main__":
    main()
