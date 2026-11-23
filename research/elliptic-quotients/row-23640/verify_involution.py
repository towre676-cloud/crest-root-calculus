"""Verify the degree-two elliptic quotient and the 7/13 Prym certificate for row 23640.

Run with Python 3.11+ and SymPy: python verify_involution.py
If given a path, the input JSON must contain this row's original 'terms' list.
No external databases, point counts or finite-field packages are used.
"""
import argparse
from itertools import product
import json
from pathlib import Path

import sympy as sp


def require(condition, message):
    if not condition:
        raise ArithmeticError(message)


def verify_equations(terms):
    x, y, z, u, v, t = sp.symbols("x y z u v t")
    source = sum(coef*x**i*y**j*z**k for i, j, k, coef in terms)
    transform = sp.Matrix([[1, 0, 1], [0, 0, 1], [0, 1, -1]])
    require(abs(transform.det()) == 1, "Coordinate change is not unimodular")
    transformed = source.subs(dict(zip((x, y, z), transform*sp.Matrix([u, v, t]))))
    q = u**2 + 7*u*v + 2*v**2
    b = u*v*(u**2+3*u*v+v**2)
    f = 3*t**4 + q*t**2 + b
    require(sp.expand(transformed-f) == 0, "Input does not give the expected double cover")
    d = u**4 + 2*u**3*v + 17*u**2*v**2 + 16*u*v**3 + 4*v**4
    require(sp.expand((6*t**2+q)**2-d-12*f) == 0, "Quotient identity failed")
    # This involution is [x:y:z] -> [x-2y:-y:z+2y] in the original coordinates.
    require(sp.expand(source.subs({x: x-2*y, y: -y, z: z+2*y}, simultaneous=True)-source) == 0,
            "Original-coordinate involution failed")
    disc = int(sp.discriminant(d.subs(v, 1), u))
    require(disc == 57600, "Branch quartic is not the expected separable polynomial")
    a0, a1, a2, a3, a4 = sp.Poly(d.subs(v, 1), u).all_coeffs()
    I = 12*a0*a4 - 3*a1*a3 + a2**2
    J = 72*a0*a2*a4 + 9*a1*a2*a3 - 27*a0*a3**2 - 27*a1**2*a4 - 2*a2**3
    j = 6912*I**3/(4*I**3-J**2)
    require(j == sp.Rational(13997521, 225), "Unexpected elliptic j-invariant")
    # The Euler identity applies at p = 7 and 13. The three affine charts test
    # all projective geometric points by the unit-ideal criterion.
    gradients = [f.diff(var) for var in (u, v, t)]
    for prime in (7, 13):
        for var in (u, v, t):
            coords = [other for other in (u, v, t) if other != var]
            basis = sp.groebner([g.subs(var, 1) for g in gradients], *coords, modulus=prime)
            require(list(basis) == [1], f"Singular reduction in chart {var}=1 mod {prime}")
    return {"double_cover_discriminant": disc, "I": int(I), "J": int(J), "j": str(j)}


class SmallField:
    """Polynomial-basis arithmetic. Only extensions of degree at most three occur."""

    def __init__(self, p, n):
        if n not in (1, 2, 3):
            raise ValueError("Only degrees 1, 2 and 3 are implemented")
        self.p, self.n, self.order = p, n, p**n
        # For degree two or three, a monic polynomial is irreducible iff it
        # has no base-field root. The tuple stores nonleading coefficients.
        self.modulus = (0,)
        if n > 1:
            self.modulus = next(c for c in product(range(p), repeat=n)
                                if all((a**n + sum(b*a**i for i, b in enumerate(c))) % p
                                       for a in range(p)))
        self.digits = [tuple((a//p**i) % p for i in range(n)) for a in range(self.order)]

    def encode(self, digits):
        return sum((d % self.p)*self.p**i for i, d in enumerate(digits))

    def add(self, a, b):
        return self.encode(x+y for x, y in zip(self.digits[a], self.digits[b]))

    def scale(self, k, a):
        return self.encode(k*x for x in self.digits[a])

    def mul(self, a, b):
        coeffs = [0]*(2*self.n-1)
        for i, x in enumerate(self.digits[a]):
            for j, y in enumerate(self.digits[b]):
                coeffs[i+j] += x*y
        for k in range(2*self.n-2, self.n-1, -1):
            for i, c in enumerate(self.modulus):
                coeffs[k-self.n+i] -= coeffs[k]*c
        return self.encode(coeffs[:self.n])


def count_curves(p, degree):
    field = SmallField(p, degree)
    add, scale, mul = field.add, field.scale, field.mul
    roots = [[] for _ in range(field.order)]
    for r in range(field.order):
        roots[mul(r, r)].append(r)
    counts = [0, 0]
    inv6 = pow(6, -1, p)
    # [u:v:t]=[0:0:1] is not on C. Normalize [u:v] on the base to [a:1]
    # or [1:0]. Each quadratic solution W supplies t²=(W-Q)/6.
    for a in range(field.order+1):
        if a == field.order:
            q, b = 1, 0
        else:
            a2 = mul(a, a)
            q = add(add(a2, scale(7, a)), 2 % p)
            b = mul(a, add(add(a2, scale(3, a)), 1))
        d = add(mul(q, q), scale(-12, b))
        counts[1] += len(roots[d])
        for w in roots[d]:
            square = scale(inv6, add(w, scale(-1, q)))
            counts[0] += len(roots[square])
    return counts


def verify_reductions():
    T, U = sp.symbols("T U")
    details = {}
    prym = {}
    for p in (7, 13):
        counts = [count_curves(p, n) for n in (1, 2, 3)]
        s1, s2, s3 = [p**n+1-pair[0] for n, pair in enumerate(counts, start=1)]
        a1 = -s1
        a2 = sp.Rational(s1**2-s2, 2)
        a3 = -sp.Rational(s1**3-3*s1*s2+2*s3, 6)
        charpoly = T**6+a1*T**5+a2*T**4+a3*T**3+p*a2*T**2+p**2*a1*T+p**3
        trace = p + 1 - counts[0][1]
        e_poly = T**2-trace*T+p
        # Verify all three quotient counts against its characteristic polynomial.
        require(counts[1][1] == p**2+1-(trace**2-2*p), "E count over degree 2 disagrees")
        require(counts[2][1] == p**3+1-(trace**3-3*p*trace), "E count over degree 3 disagrees")
        h, rem = sp.div(charpoly, e_poly, T)
        require(rem == 0, "E polynomial does not divide the Jacobian polynomial")
        prym[p] = h
        details[str(p)] = {
            "C_counts": [r[0] for r in counts],
            "E_counts": [r[1] for r in counts],
            "Jacobian": str(sp.expand(charpoly)),
            "elliptic_quotient": str(e_poly),
            "Prym": str(sp.factor(h)),
        }
    f1, f2 = T**2+4*T+13, T**2-2*T+13
    require(sp.expand(prym[13]-f1*f2) == 0, "Unexpected two Prym factors at 13")
    require(sp.discriminant(f1, T) == -36 and sp.discriminant(f2, T) == -48,
            "Ordinary elliptic endomorphism fields disagree with certificate")
    powered = sp.Poly(sp.resultant(prym[7], U-T**2, T), U, domain=sp.QQ)
    require(powered.as_expr() == U**4+46*U**2+2401 and powered.is_irreducible,
            "Degree-two Frobenius at 7 does not certify simplicity")
    details["7"]["Prym_over_F49"] = str(powered.as_expr())
    return details


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quartic_json", nargs="?", type=Path)
    args = parser.parse_args()
    path = args.quartic_json or Path(__file__).with_name("quartic.json")
    terms = json.loads(path.read_text(encoding="utf-8"))["terms"]
    print(json.dumps({"row": 23640, "equations": verify_equations(terms),
                      "reductions": verify_reductions()}, indent=2))


if __name__ == "__main__":
    main()
