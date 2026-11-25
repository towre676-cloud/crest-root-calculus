"""Finite-field core: Jacobian, elliptic targets, pullback f^* and pushforward f_* for certified maps.

A certified map (Proposition 6.1) is X = N/D, Y = y (N'D - ND') / ((a + b x) D^2) onto
Y^2 = 4X^3 - g2 X - g3, with coefficients in Q(i) (Gaussian rationals).  Reduction is taken at
a good prime p; when i occurs we require p = 1 mod 4 and fix i -> sqrt(-1) in F_p.

Points of a reduced Mumford divisor that are not F_q-rational lie in F_{q^2}; they are handled in
the explicit quadratic extension F_q[t]/(u(t)) so that no field embedding has to be chosen.
"""
from __future__ import annotations

import itertools
import sys
from fractions import Fraction

import flint
import sympy

sys.path.insert(0, "../exact-completeness")
from jac import Jacobian, odd_model, weil_poly, jac_order  # noqa: E402


# ------------------------------------------------------------------ Gaussian rationals
def gauss(s):
    """Parse 'a + b*w' style strings (w = i) into (Fraction, Fraction)."""
    w = sympy.symbols("w")
    e = sympy.expand(sympy.sympify(str(s).replace("^", "**"), locals={"w": w}))
    re_, im_ = e.coeff(w, 0), e.coeff(w, 1)
    fr = lambda z: Fraction(int(sympy.fraction(z)[0]), int(sympy.fraction(z)[1]))
    return (fr(re_), fr(im_))


def poly_gauss(s):
    x, w = sympy.symbols("x w")
    P = sympy.Poly(sympy.expand(sympy.sympify(str(s).replace("^", "**"), locals={"x": x, "w": w})), x)
    out = []
    for c in P.all_coeffs()[::-1]:
        out.append(gauss(c))
    return out


def nf_root_mod_p(nf, p):
    """A simple root of the integer polynomial nf (low to high) modulo p, or None."""
    for t in range(p):
        if sum(c * pow(t, k, p) for k, c in enumerate(nf)) % p == 0:
            d = sum(k * c * pow(t, k - 1, p) for k, c in enumerate(nf) if k) % p
            if d:
                return t
    return None


class Field:
    """F_{p^r} with a chosen reduction theta -> root of nf mod p for coefficients in Q(theta).

    Coefficients are Fractions or tuples (c_0, ..., c_{n-1}) in the power basis of theta.
    Gaussian rationals are the case nf = t^2 + 1 (theta = i)."""

    def __init__(self, p, r, need_i=False, nf=None):
        self.p, self.r = p, r
        self.F = flint.fq_default_ctx(p, r)
        self.R = flint.fq_default_poly_ctx(self.F)
        self.x = self.R.gen()
        if need_i and nf is None:
            nf = [1, 0, 1]
        self.nf = nf
        self.theta = None
        self.i = None
        if nf is not None:
            root = nf_root_mod_p(nf, p)
            if root is None:
                raise ValueError("field polynomial has no simple root mod p")
            self.theta = self.F(root)
            if nf == [1, 0, 1]:
                self.i = self.theta

    def _fr(self, a):
        a = Fraction(a)
        return self.F(a.numerator % self.p) / self.F(a.denominator % self.p)

    def q(self, v):
        if isinstance(v, tuple):
            out = self.F(0)
            pw = self.F(1)
            for k, c in enumerate(v):
                if Fraction(c):
                    if k and self.theta is None:
                        raise ValueError("A nonconstant number-field coefficient needs a field polynomial.")
                    out += self._fr(c) * pw
                pw = pw * self.theta if self.theta is not None else pw
            return out
        return self._fr(v)

    def poly(self, cs):
        return self.R([self.q(c) for c in cs])


# ------------------------------------------------------------------ quadratic extension F_q[t]/(u)
class QE:
    """Element a + b t of F_q[t]/(t^2 - s1 t + s0)."""
    __slots__ = ("a", "b", "s1", "s0")

    def __init__(self, a, b, s1, s0):
        self.a, self.b, self.s1, self.s0 = a, b, s1, s0

    def _c(self, o):
        if isinstance(o, QE):
            return o
        return QE(self.a * 0 + o, self.a * 0, self.s1, self.s0)

    def __add__(self, o):
        o = self._c(o)
        return QE(self.a + o.a, self.b + o.b, self.s1, self.s0)

    __radd__ = __add__

    def __neg__(self):
        return QE(-self.a, -self.b, self.s1, self.s0)

    def __sub__(self, o):
        return self + (-self._c(o))

    def __rsub__(self, o):
        return self._c(o) - self

    def __mul__(self, o):
        o = self._c(o)
        # t^2 = s1 t - s0
        bb = self.b * o.b
        return QE(self.a * o.a - bb * self.s0, self.a * o.b + self.b * o.a + bb * self.s1, self.s1, self.s0)

    __rmul__ = __mul__

    def conj(self):
        # t -> s1 - t
        return QE(self.a + self.b * self.s1, -self.b, self.s1, self.s0)

    def norm(self):
        n = self * self.conj()
        return n.a

    def inv(self):
        n = self.norm()
        c = self.conj()
        return QE(c.a / n, c.b / n, self.s1, self.s0)

    def __truediv__(self, o):
        o = self._c(o)
        return self * o.inv()

    def __rtruediv__(self, o):
        return self._c(o) * self.inv()

    def __pow__(self, n):
        r = self._c(1)
        b = self
        while n:
            if n & 1:
                r = r * b
            b = b * b
            n >>= 1
        return r

    def is_zero(self):
        return self.a.is_zero() and self.b.is_zero()

    def __eq__(self, o):
        return (self - o).is_zero()

    def base(self):
        assert self.b.is_zero()
        return self.a


def is_zero(z):
    return z.is_zero()


# ------------------------------------------------------------------ elliptic curves Y^2 = 4X^3 - g2 X - g3
class Ell:
    def __init__(self, K, g2, g3):
        self.K, self.g2, self.g3 = K, K.q(g2), K.q(g3)

    def on(self, P):
        return P is None or is_zero(P[1] * P[1] - (4 * P[0] ** 3 - self.g2 * P[0] - self.g3))

    def add(self, P, Q):
        if P is None:
            return Q
        if Q is None:
            return P
        (x1, y1), (x2, y2) = P, Q
        if is_zero(x1 - x2):
            if is_zero(y1 + y2):
                return None
            lam = (12 * x1 * x1 - self.g2) / (2 * y1)
        else:
            lam = (y2 - y1) / (x2 - x1)
        x3 = lam * lam / 4 - x1 - x2
        return (x3, -(y1 + lam * (x3 - x1)))

    def neg(self, P):
        return None if P is None else (P[0], -P[1])

    def mul(self, n, P):
        if n < 0:
            return self.mul(-n, self.neg(P))
        R = None
        while n:
            if n & 1:
                R = self.add(R, P)
            P = self.add(P, P)
            n >>= 1
        return R

    def eq(self, P, Q):
        if P is None or Q is None:
            return P is None and Q is None
        return is_zero(P[0] - Q[0]) and is_zero(P[1] - Q[1])

    def key(self, P):
        return None if P is None else (str(P[0]), str(P[1]))

    def points(self, limit=None):
        F, p, r = self.K.F, self.K.p, self.K.r
        g = F.gen()
        count = 0
        for t in itertools.product(range(p), repeat=r):
            X = sum((c * g ** i for i, c in enumerate(t) if c), F(0))
            v = 4 * X ** 3 - self.g2 * X - self.g3
            if v.is_zero():
                yield (X, v)
            elif v.is_square():
                s = v.sqrt()
                yield (X, s)
                yield (X, -s)
            count += 1
            if limit and count > limit:
                return


# ------------------------------------------------------------------ certified map on the reduced curve
class CertMap:
    def __init__(self, K, J, alpha, f, N, D, a, b, g2, g3):
        self.K, self.J, self.alpha = K, J, K.F(alpha)
        self.f = K.poly(f)
        self.N, self.D = K.poly(N), K.poly(D)
        self.ab = K.poly([a, b])
        self.W = self.N.derivative() * self.D - self.N * self.D.derivative()
        self.E = Ell(K, g2, g3)
        self.deg = max(self.N.degree(), self.D.degree())
        if self.ab.is_zero() or self.W.is_zero() or not self.N.gcd(self.D).is_one():
            raise ValueError("The map has degenerate denominators or differential.")
        # identity (6.1) on the reduction
        lhs = self.f * self.W * self.W
        rhs = self.ab * self.ab * self.D * (4 * self.N ** 3 - self.E.g2 * self.N * self.D ** 2 - self.E.g3 * self.D ** 3)
        assert (lhs - rhs).is_zero(), "identity (6.1) fails mod p"
        self.WA = self._apply_point(self.alpha, None)  # image of the odd-model base point (alpha, 0)

    # value of the map at an original-model point (x0, y0), coordinates in F_q or QE
    def _apply_point(self, x0, y0):
        Nv, Dv = _ev(self.N, x0), _ev(self.D, x0)
        if is_zero(Dv):
            return None
        X = Nv / Dv
        if y0 is None:          # Weierstrass point: y = 0, image is 2-torsion with Y = 0
            return (X, X * 0)
        Wv = _ev(self.W, x0)
        abv = _ev(self.ab, x0)
        den = abv * Dv * Dv
        if is_zero(den):
            raise ZeroDivisionError("ramified point")
        return (X, y0 * Wv / den)

    # ----- pullback E -> J (class of f^{-1}(P) minus that of f^{-1}(R) translate)
    def fibre(self, P):
        X0, Y0 = P
        u = self.N - X0 * self.D
        if u.degree() != self.deg:
            return None
        u = u.monic()
        if not u.is_squarefree():
            return None
        if not u.gcd(self.W).is_one() or not u.gcd(self.D).is_one() or not u.gcd(self.ab).is_one():
            return None
        if u(self.alpha).is_zero():
            return None
        inv = self.W.xgcd(u)[1]
        y = (Y0 * self.ab * self.D * self.D * inv) % u
        R = self.K.R
        xs = R.gen()
        d = u.degree()
        up = R(0)
        for i, c in enumerate(u.coeffs()):
            up += c * (self.alpha * xs + 1) ** i * xs ** (d - i)
        if up.degree() != d:
            return None
        up = up.monic()
        vp = R(0)
        for i, c in enumerate(y.coeffs()):
            vp += c * (self.alpha * xs + 1) ** i * xs ** (3 - i)
        vp = vp % up
        assert ((vp * vp - self.J.f) % up).is_zero()
        return self.J._reduce(up, vp)

    def pullback(self, P, aux):
        """f^*(P) via f^{-1}(P + R) - f^{-1}(R) for the first R in aux that is generic for both."""
        if P is None:
            return self.J.zero
        for R in aux:
            S = self.E.add(P, R)
            if S is None:
                continue
            A, B = self.fibre(S), self.fibre(R)
            if A is not None and B is not None:
                return self.J.add(A, self.J.neg(B))
        raise RuntimeError("no generic translate in aux list")

    # ----- pushforward J -> E of a reduced divisor class [D - deg(D) infinity']
    def _push_raw(self, Dv):
        u, v = Dv
        if u.degree() == 0:
            return None
        cs = u.coeffs()
        total = None
        if u.degree() == 1:
            xp = -cs[0]
            pts = [(xp, v(xp))]
            imgs = [self._odd_to_image(xp, yp) for xp, yp in pts]
        else:
            fac = u.factor()[1]
            if len(fac) == 2 or (len(fac) == 1 and fac[0][1] == 2):
                roots = [r for r, m in u.roots() for _ in range(m)]
                imgs = [self._odd_to_image(r, v(r)) for r in roots]
            else:
                s1, s0 = -cs[1], cs[0]   # u = t^2 - s1 t + s0
                t = QE(self.K.F(0), self.K.F(1), s1, s0)
                cv = v.coeffs()
                zero = self.K.F(0)
                vt = QE(cv[0] if len(cv) > 0 else zero, cv[1] if len(cv) > 1 else zero, s1, s0)
                P1 = self._odd_to_image(t, vt)
                P2 = None if P1 is None else (P1[0].conj(), P1[1].conj())
                S = self.E.add(P1, P2)
                if S is None:
                    imgs = [None]
                else:
                    imgs = [(S[0].base(), S[1].base())]
        for I in imgs:
            total = self.E.add(total, I)
        # subtract deg * f(base point)
        return self.E.add(total, self.E.neg(self.E.mul(u.degree(), self.WA)))

    def _odd_to_image(self, xp, yp):
        if is_zero(xp):
            raise ZeroDivisionError("point over x = infinity")
        x0 = self.alpha + 1 / xp
        y0 = yp / (xp * xp * xp)
        return self._apply_point(x0, y0)

    def pushforward(self, Dcls, aux_classes):
        try:
            return self._push_raw(Dcls)
        except ZeroDivisionError:
            pass
        for G in aux_classes:
            try:
                A = self._push_raw(self.J.add(Dcls, G))
                B = self._push_raw(G)
                return self.E.add(A, self.E.neg(B))
            except ZeroDivisionError:
                continue
        raise RuntimeError("pushforward: no generic translate")


def _ev(poly, z):
    if not isinstance(z, QE):
        return poly(z)
    r = z._c(0)
    for c in reversed(poly.coeffs()):
        r = r * z + c
    return r
