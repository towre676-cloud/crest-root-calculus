"""Class providers: global morphisms C -> E, known through exact reductions of f^* and f_*.

A provider has a target elliptic curve E (over F_q, format Y^2 = 4X^3 - g2 X - g3), a degree, and
the maps pull: E(F_q) -> J(F_q), push: J(F_q) -> E(F_q).  Certified maps are providers.  So are
sums g = f_1 + rho^ o f_2 with rho the reduction of a global isogeny E_1 -> E_2:
    g^* = f_1^* + f_2^* rho,     g_* = f_{1*} + rho^ f_{2*},
both global, hence g defines a class in NS(J_Qbar) whose reduction is computed exactly.
"""
from __future__ import annotations

from fractions import Fraction

import cypari2
import sympy

from ffcore import Ell, _ev

pari = cypari2.Pari()


class MapProvider:
    def __init__(self, R, i):
        self.R, self.i = R, i
        self.E = R.maps[i].E
        self.deg = R.maps[i].deg

    def pull(self, P):
        return self.R.pull(self.i, P)

    def push(self, D):
        return self.R.push(self.i, D)

    def aux(self, n=40):
        return self.R.aux_points(self.i, n)


def short_ab(g2g3):
    g2, g3 = g2g3
    return -Fraction(g2) / 4, -Fraction(g3) / 4


def pq(fr):
    fr = Fraction(fr)
    return pari(f"{fr.numerator}/{fr.denominator}")


class IsogenyPair:
    """Reduction of a global cyclic isogeny rho: E_s -> E_d and its dual, both over Q.

    PARI's ellisomat on the short model of E_s returns, for the k-isogenous class member L,
    the map E_s -> L and the dual map L -> E_s.  L is identified with the target model E_d by
    (x, y) -> (u^2 x, u^3 y); u^2 is rational, u is taken in F_q.  The sign of the dual is
    fixed by the self-test rho^ rho = [k] on points.
    """

    def __init__(self, K, src, dst, k):
        self.K, self.k = K, k
        As, Bs = short_ab(src)
        Ad, Bd = short_ab(dst)
        E = pari.ellinit([0, 0, 0, pq(As), pq(Bs)])
        jd = pari.ellinit([0, 0, 0, pq(Ad), pq(Bd)]).j()
        L, M = pari.ellisomat(E)[:2]
        idx = None
        for t in range(len(L)):
            if int(M[0, t]) == k and pari.ellinit(L[t][0]).j() == jd:
                idx = t
                break
        if idx is None:
            raise ValueError(f"no rational {k}-isogeny to the target j")
        Lc, fwd, bwd = L[idx][0], L[idx][1], L[idx][2]
        AL, BL = Fraction(str(Lc[0])), Fraction(str(Lc[1]))
        if AL != 0 and BL != 0:
            u2 = (Bd * AL) / (BL * Ad)
        elif BL == 0:
            raise NotImplementedError("j = 1728 identification: use u^4 = Ad/AL")
        else:
            raise NotImplementedError("j = 0 identification: use u^6 = Bd/BL")
        if u2*u2*AL != Ad or u2**3*BL != Bd:
            raise ArithmeticError("The target identification fails over Q.")
        self.u2 = K.q(u2)
        if not self.u2.is_square():
            raise ValueError("u^2 is not a square in this field")
        self.u = self.u2.sqrt()
        self.fx, self.fy, self.fk = self._maps(fwd)
        self.bx, self.by, self.bk = self._maps(bwd)
        self.src, self.dst = Ell(K, src[0], src[1]), Ell(K, dst[0], dst[1])
        self.sign = 1

    def _poly(self, g):
        g = pari.subst(g, "y", 1)
        deg = int(pari.poldegree(g, "x"))
        cs = [Fraction(str(pari.polcoef(g, i, "x"))) for i in range(deg + 1)]
        return self.K.poly(cs)

    def _maps(self, m):
        return self._poly(m[0]), self._poly(m[1]), self._poly(m[2])

    @staticmethod
    def _apply(P, xm, ym, km):
        X, Y = P
        x, y = X, Y / 2
        kv = km(x)
        if kv.is_zero():
            return None
        return (xm(x) / (kv * kv), y * ym(x) / (kv * kv * kv))

    def forward(self, P):
        if P is None:
            return None
        Q = self._apply(P, self.fx, self.fy, self.fk)
        if Q is None:
            return None
        x, y = Q
        out = (self.u2 * x, 2 * self.u2 * self.u * y)
        assert self.dst.on(out), "forward image off the target model"
        return out

    def backward(self, Q):
        if Q is None:
            return None
        X, Y = Q
        x, y = X / self.u2, (Y / 2) / (self.u2 * self.u)
        P = self._apply((x, 2 * y), self.bx, self.by, self.bk)
        if P is None:
            return None
        out = (P[0], 2 * P[1])
        assert self.src.on(out), "dual image off the source model"
        return out if self.sign == 1 else self.src.neg(out)

    def calibrate(self, points):
        for P in points:
            Q = self.backward(self.forward(P))
            kP = self.src.mul(self.k, P)
            if self.src.eq(Q, kP):
                continue
            if self.src.eq(Q, self.src.neg(kP)):
                self.sign = -self.sign
                continue
            raise AssertionError("rho^ rho is not [+-k]")
        for P in points:
            assert self.src.eq(self.backward(self.forward(P)), self.src.mul(self.k, P))


class SumProvider:
    """g = f_1 + rho^ o f_2 : C -> E_1 for providers f_1 (target E_1), f_2 (target E_2) and rho: E_1 -> E_2."""

    def __init__(self, f1, f2, rho, rho_dual):
        self.f1, self.f2, self.rho, self.rhod = f1, f2, rho, rho_dual
        self.E = f1.E
        self.deg = None

    def pull(self, P):
        J = self.f1.R.J
        return J.add(self.f1.pull(P), self.f2.pull(self.rho.forward(P)))

    def push(self, D):
        return self.E.add(self.f1.push(D), self.rho.backward(self.f2.push(D)))

    def aux(self, n=40):
        return self.f1.aux(n)


class IsoPair:
    """Isomorphism rho(X, Y) = (u^2 X, u^3 Y) between (g2, g3) models, optionally composed with the
    CM automorphism [i](X, Y) = (-X, i Y) on the target (j = 1728).  backward is the inverse."""

    def __init__(self, K, src, dst, u2, cm_i=False):
        self.K = K
        self.src, self.dst = Ell(K, src[0], src[1]), Ell(K, dst[0], dst[1])
        if u2*u2*src[0] != dst[0] or u2**3*src[1] != dst[1] or (cm_i and dst[1]):
            raise ArithmeticError("The target isomorphism or CM action does not verify.")
        self.u2 = K.q(u2)
        if not self.u2.is_square():
            raise ValueError("The target isomorphism requires a larger residue field.")
        self.u = self.u2.sqrt()
        self.cm = cm_i
        self.k = 1

    def forward(self, P):
        if P is None:
            return None
        Q = (self.u2 * P[0], self.u2 * self.u * P[1])
        if self.cm:
            Q = (-Q[0], self.K.i * Q[1])
        assert self.dst.on(Q)
        return Q

    def backward(self, Q):
        if Q is None:
            return None
        if self.cm:
            Q = (-Q[0], -self.K.i * Q[1])
        P = (Q[0] / self.u2, Q[1] / (self.u2 * self.u))
        assert self.src.on(P)
        return P
