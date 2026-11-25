"""Compute the kernel of a sum of elliptic pullbacks on finite torsion."""

from fractions import Fraction

import flint

from absence import good_model
from jac import Jacobian, _all_elements, odd_model


class Field:
    def __init__(self, p, r):
        self.p, self.r = p, r
        self.F = flint.fq_default_ctx(p, r)
        self.R = flint.fq_default_poly_ctx(self.F)
        self.x = self.R.gen()

    def q(self, value):
        value = Fraction(value)
        return self.F(value.numerator % self.p)/self.F(value.denominator % self.p)

    def poly(self, coefficients):
        return self.R([self.q(c) for c in coefficients])


class Elliptic:
    def __init__(self, field, g2, g3):
        self.field = field
        self.g2, self.g3 = field.q(g2), field.q(g3)

    def rhs(self, x):
        return 4*x**3-self.g2*x-self.g3

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        x1, y1 = P
        x2, y2 = Q
        if x1 == x2:
            if (y1+y2).is_zero():
                return None
            lam = (12*x1*x1-self.g2)/(2*y1)
        else:
            lam = (y2-y1)/(x2-x1)
        x3 = lam*lam/4-x1-x2
        return x3, -(y1+lam*(x3-x1))

    def mul(self, n, P):
        result = None
        while n:
            if n & 1:
                result = self.add(result, P)
            P = self.add(P, P)
            n >>= 1
        return result

    def affine_points(self):
        for x in _all_elements(self.field.F, self.field.p, self.field.r):
            rhs = self.rhs(x)
            if rhs.is_zero():
                yield x, rhs
            elif rhs.is_square():
                y = rhs.sqrt()
                yield x, y
                yield x, -y

    def torsion(self, n):
        """Find all E[2] or E[4] by factoring the doubling equations."""
        if n not in (2, 4):
            raise ValueError("Supported torsion levels are 2 and 4.")
        x = self.field.x
        cubic = self.rhs(x)
        roots = [r for r, _ in cubic.roots()]
        points = [None, *((r, self.field.F(0)) for r in roots)]
        if n == 4:
            for r in roots:
                doubling = (12*x*x-self.g2)**2-16*(2*x+r)*cubic
                for t, _ in doubling.roots():
                    rhs = self.rhs(t)
                    if not rhs.is_zero() and rhs.is_square():
                        y = rhs.sqrt()
                        points.extend(((t, y), (t, -y)))
        if len({point_key(P) for P in points}) != n*n or any(self.mul(n, P) is not None for P in points):
            raise ValueError(f"Full E[{n}] is unavailable over the candidate field.")
        return points


def point_key(P):
    return None if P is None else tuple(map(str, P))


class Pullback:
    def __init__(self, field, jacobian, alpha, f, record):
        self.field, self.jacobian, self.alpha = field, jacobian, alpha
        self.N, self.D = field.poly(record["N"]), field.poly(record["D"])
        self.f = field.poly(f)
        self.differential = field.poly([record["a"], record["b"]])
        self.degree = record["degree"]
        self.W = self.N.derivative()*self.D-self.N*self.D.derivative()
        if (not self.N.gcd(self.D).is_one() or max(self.N.degree(), self.D.degree()) != self.degree
                or self.differential.is_zero() or self.W.is_zero()):
            raise ValueError("Reduction does not preserve the horizontal map and differential.")
        g2, g3 = field.q(record["g2"]), field.q(record["g3"])
        curve = 4*self.N**3-g2*self.N*self.D**2-g3*self.D**3
        if g2**3 == 27*g3**2 or self.f*self.W**2 != self.differential**2*self.D*curve:
            raise ValueError("Reduction does not give the certified map on smooth curves.")
        self.cached = {}

    def fibre(self, P):
        key = point_key(P)
        if key not in self.cached:
            self.cached[key] = self._fibre(P)
        return self.cached[key]

    def _fibre(self, P):
        if P is None:
            return None
        X, Y = P
        u = self.N-X*self.D
        if u.degree() != self.degree or not u.is_squarefree():
            return None
        u = u.monic()
        g, inverse, _ = self.W.xgcd(u)
        if not g.is_one():
            return None
        v = (Y*self.differential*self.D**2*inverse) % u
        if not ((v*v-self.f) % u).is_zero():
            raise ArithmeticError("The computed fibre is off the source curve.")
        if self.alpha is not None:
            x = self.field.x
            z = self.alpha*x+1
            u = sum((c*z**i*x**(self.degree-i) for i, c in enumerate(u.coeffs())), self.field.R(0))
            if u.degree() != self.degree:
                return None
            u = u.monic()
            v = sum((c*z**i*x**(3-i) for i, c in enumerate(v.coeffs())), self.field.R(0)) % u
        if not ((v*v-self.jacobian.f) % u).is_zero():
            raise ArithmeticError("The transformed fibre is off the odd model.")
        return self.jacobian._reduce(u, v)

    def image(self, E, P):
        for R in E.affine_points():
            A, B = self.fibre(E.add(P, R)), self.fibre(R)
            if A is not None and B is not None:
                return self.jacobian.add(A, self.jacobian.neg(B))
        raise ValueError("No suitable fibre translate exists over the candidate field.")


def certify_kernel(f, maps, p, r, n):
    if not good_model(f, p):
        raise ValueError("Use a certified good reduction of the source.")
    if any((m["g2"]**3-27*m["g3"]**2).numerator % p == 0 for m in maps):
        raise ValueError("Use a good reduction of both elliptic targets.")
    values = [c for m in maps for c in (*m["N"], *m["D"], m["g2"], m["g3"], m["a"], m["b"])]
    if any(Fraction(c).denominator % p == 0 for c in values):
        raise ValueError("A map coefficient has a pole at the candidate prime.")
    field = Field(p, r)
    g, alpha = odd_model(f, p)
    J = Jacobian(g, p, r)
    elliptic, groups, images = [], [], []
    for record in maps:
        E = Elliptic(field, record["g2"], record["g3"])
        points = E.torsion(n)
        pullback = Pullback(field, J, alpha, f, record)
        image = {point_key(P): J.zero if P is None else pullback.image(E, P) for P in points}
        for P in points:
            for Q in points:
                actual = image[point_key(E.add(P, Q))]
                expected = J.add(image[point_key(P)], image[point_key(Q)])
                if J.key(actual) != J.key(expected):
                    raise ArithmeticError("Pullback values fail the full additivity check.")
        elliptic.append(E)
        groups.append(points)
        images.append(image)
    kernel = [(P, Q) for P in groups[0] for Q in groups[1]
              if J.is_zero(J.add(images[0][point_key(P)], images[1][point_key(Q)]))]
    return field, elliptic, groups, kernel
