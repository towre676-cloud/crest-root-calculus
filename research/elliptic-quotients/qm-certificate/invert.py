"""Numerical Jacobi inversion and the morphism P -> alpha(P - W0) for a genus-two curve.

Abel-Jacobi values are integrals of (dx/y, x dx/y) from a Weierstrass point W0 = (r0, 0) along
straight paths, with the branch of y continued from r0.  A divisor P1 + P2 with
AJ(P1) + AJ(P2) = z mod Lambda is found by Newton's method in (x1, x2) with sheet signs, reducing
the residual modulo the period lattice at every step.  Random restarts supply the global search;
the result is accepted only below a residual threshold and refined at high precision.
"""
from __future__ import annotations

import random
import sys

import mpmath as mp

sys.path.insert(0, "../period-reconstruction")
from periods import period_matrix, roots_of, nearest_point  # noqa: E402


class Curve:
    def __init__(self, f, dps=40):
        mp.mp.dps = 70
        self.f = f
        self.Pi, self.E, self.S = period_matrix(f)
        ft = list(f)
        while ft[-1] == 0:
            ft.pop()
        self.roots = roots_of([mp.mpf(c) for c in ft])
        self.lc = mp.mpf(f[6]) if len(f) > 6 and f[6] else mp.mpf(f[5])
        self.r0 = min(self.roots, key=lambda z: abs(mp.im(z)))
        P = self.Pi
        self.Pr = mp.matrix(4, 4)
        for c in range(4):
            self.Pr[0, c], self.Pr[1, c] = mp.re(P[0, c]), mp.im(P[0, c])
            self.Pr[2, c], self.Pr[3, c] = mp.re(P[1, c]), mp.im(P[1, c])
        self.Pri = self.Pr ** -1
        mp.mp.dps = dps

    def fx(self, x):
        s = mp.mpf(0)
        for c in reversed(self.f):
            s = s * x + c
        return s

    def aj(self, xs):
        """(int_{r0}^{xs} dx/y, x dx/y) and the endpoint value of y on the continued branch."""
        r0 = self.r0
        rest = [z for z in self.roots if abs(z - r0) > mp.mpf(10) ** -20]
        mids = [nearest_point(r0, xs, z) for z in rest]
        sq = mp.sqrt(xs - r0)

        def gsqrt(x):
            q = mp.sqrt(self.lc)
            for z, m in zip(rest, mids):
                q *= mp.sqrt(m - z) * mp.sqrt((x - z) / (m - z))
            return q
        vals = []
        for m in (0, 1):
            vals.append(mp.quad(lambda s: (r0 + (xs - r0) * s * s) ** m * 2 * sq / gsqrt(r0 + (xs - r0) * s * s),
                                [0, mp.mpf(1) / 2, 1]))
        yend = sq * gsqrt(xs)
        return vals, yend

    def reduce(self, w):
        v = mp.matrix([mp.re(w[0]), mp.im(w[0]), mp.re(w[1]), mp.im(w[1])])
        c = self.Pri * v
        n = [mp.nint(c[i]) for i in range(4)]
        shift = [sum(self.Pi[r, k] * n[k] for k in range(4)) for r in range(2)]
        return [w[0] - shift[0], w[1] - shift[1]]

    def point_aj(self, x, y):
        vals, yend = self.aj(x)
        s = 1 if abs(yend - y) < abs(yend + y) else -1
        return [s * vals[0], s * vals[1]]

    def invert(self, z, rng, restarts=60, iters=40, tol=None):
        tol = tol or mp.mpf(10) ** (-(mp.mp.dps - 8))
        scale = max(abs(r) for r in self.roots) + 1
        best = None
        for _ in range(restarts):
            xs = [mp.mpc(rng.uniform(-1, 1), rng.uniform(-1, 1)) * scale for _ in range(2)]
            sg = [rng.choice([1, -1]) for _ in range(2)]
            ok = False
            for it in range(iters):
                try:
                    A = [self.aj(x) for x in xs]
                except (ZeroDivisionError, ValueError):
                    break
                w = [sg[0] * A[0][0][k] + sg[1] * A[1][0][k] - z[k] for k in range(2)]
                w = self.reduce(w)
                res = abs(w[0]) + abs(w[1])
                if res < tol:
                    ok = True
                    break
                y = [sg[i] * A[i][1] for i in range(2)]
                # Jacobian columns d/dx_i of (dx/y, x dx/y) integrals
                a, b = 1 / y[0], 1 / y[1]
                c, d = xs[0] / y[0], xs[1] / y[1]
                det = a * d - b * c
                if abs(det) < mp.mpf(10) ** -30:
                    break
                dx0 = (d * w[0] - b * w[1]) / det
                dx1 = (-c * w[0] + a * w[1]) / det
                step = max(abs(dx0), abs(dx1))
                damp = 1 if step < scale else scale / step
                xs = [xs[0] - damp * dx0, xs[1] - damp * dx1]
            if ok:
                ys = [sg[i] * self.aj(xs[i])[1] for i in range(2)]
                return xs, ys, res
            if best is None:
                best = None
        return None


# ------------------------------------------------------------------ fast double-precision search
import numpy as np

_GL_X, _GL_W = np.polynomial.legendre.leggauss(24)
_NSUB = 10
_S = np.concatenate([(k + (_GL_X + 1) / 2) / _NSUB for k in range(_NSUB)])
_WS = np.concatenate([_GL_W / 2 / _NSUB for k in range(_NSUB)])


class FastCurve:
    """Double-precision mirror of Curve for the global Newton search."""

    def __init__(self, C):
        self.C = C
        self.roots = np.array([complex(r) for r in C.roots])
        self.r0 = complex(C.r0)
        self.lc = float(C.lc)
        self.Pi = np.array([[complex(C.Pi[r, c]) for c in range(4)] for r in range(2)])
        self.Pr = np.array([[float(C.Pr[r, c]) for c in range(4)] for r in range(4)])
        self.Pri = np.linalg.inv(self.Pr)
        self.rest = np.array([z for z in self.roots if abs(z - self.r0) > 1e-9])
        self.span = max(abs(a - b) for a in self.roots for b in self.roots)

    def clearance(self, xs):
        r0 = self.r0
        d = xs - r0
        t = np.clip(np.real((self.rest - r0) * np.conj(d)) / abs(d) ** 2, 0, 1)
        return np.min(abs(r0 + t * d - self.rest))

    def aj(self, xs):
        r0 = self.r0
        d = xs - r0
        t = np.clip(np.real((self.rest - r0) * np.conj(d)) / abs(d) ** 2, 0, 1)
        mids = r0 + t * d
        if np.min(abs(mids - self.rest)) < 0.05 * self.span:
            raise ValueError("path too close to a branch point")
        x = r0 + d * _S ** 2
        q = np.sqrt(self.lc + 0j) * np.ones_like(x)
        for z, m in zip(self.rest, mids):
            q = q * np.sqrt(m - z) * np.sqrt((x - z) / (m - z))
        sq = np.sqrt(d)
        integ0 = 2 * sq / q
        v0 = np.sum(_WS * integ0)
        v1 = np.sum(_WS * x * integ0)
        xe = xs
        qe = np.sqrt(self.lc + 0j)
        for z, m in zip(self.rest, mids):
            qe = qe * np.sqrt(m - z) * np.sqrt((xe - z) / (m - z))
        return np.array([v0, v1]), sq * qe

    def reduce(self, w):
        v = np.array([w[0].real, w[0].imag, w[1].real, w[1].imag])
        n = np.round(self.Pri @ v)
        return w - self.Pi @ n

    def invert(self, z, rng, restarts=400, iters=60, scale=None):
        scale = scale or (max(abs(self.roots)) + 1)
        for _ in range(restarts):
            xs = np.array([complex(rng.uniform(-1, 1), rng.uniform(-1, 1)) * scale for _ in range(2)])
            sg = np.array([rng.choice([1, -1]) for _ in range(2)])
            for it in range(iters):
                try:
                    A0, y0 = fast_aj_routed(self, self.router, xs[0])
                    A1, y1 = fast_aj_routed(self, self.router, xs[1])
                except (ValueError, ZeroDivisionError):
                    break
                w = self.reduce(sg[0] * A0 + sg[1] * A1 - z)
                if np.max(abs(w)) < 1e-12:
                    return xs, np.array([sg[0] * y0, sg[1] * y1])
                y = [sg[0] * y0, sg[1] * y1]
                Jm = np.array([[1 / y[0], 1 / y[1]], [xs[0] / y[0], xs[1] / y[1]]])
                try:
                    dx = np.linalg.solve(Jm, w)
                except np.linalg.LinAlgError:
                    break
                st = max(abs(dx))
                if not np.isfinite(st):
                    break
                damp = 1 if st < 0.3 * scale else 0.3 * scale / st
                xs = xs - damp * dx
        return None


def refine(C, xs, ys, z, steps=8, Rt=None):
    """High-precision Newton refinement from a double-precision solution."""
    xs = [mp.mpc(complex(v)) for v in xs]
    ys = [mp.mpc(complex(v)) for v in ys]
    for _ in range(steps):
        A = [mp_aj_routed(C, Rt, x) for x in xs]
        sg = [1 if abs(A[i][1] - ys[i]) < abs(A[i][1] + ys[i]) else -1 for i in range(2)]
        w = [sg[0] * A[0][0][k] + sg[1] * A[1][0][k] - z[k] for k in range(2)]
        w = C.reduce(w)
        y = [sg[i] * A[i][1] for i in range(2)]
        a, b, c, d = 1 / y[0], 1 / y[1], xs[0] / y[0], xs[1] / y[1]
        det = a * d - b * c
        dx0 = (d * w[0] - b * w[1]) / det
        dx1 = (-c * w[0] + a * w[1]) / det
        xs = [xs[0] - dx0, xs[1] - dx1]
        ys = [y[0], y[1]]
        res = abs(w[0]) + abs(w[1])
    ys = [sg[i] * mp_aj_routed(C, Rt, xs[i])[1] for i in range(2)]
    return xs, ys, res


def solve_divisor(C, F, z_mp, rng, tries=6):
    """Fast search plus high-precision verification; returns (xs, ys) accurate to the working precision."""
    z = np.array([complex(z_mp[0]), complex(z_mp[1])])
    for _ in range(tries):
        r = F.invert(z, rng)
        if r is None:
            continue
        try:
            xs, ys = refine(C, r[0], r[1], z_mp, steps=4, Rt=F.router)[:2]
            A = [mp_aj_routed(C, F.router, x) for x in xs]
        except (ValueError, ZeroDivisionError):
            continue
        sg = [1 if abs(A[i][1] - ys[i]) < abs(A[i][1] + ys[i]) else -1 for i in range(2)]
        w = C.reduce([sg[0] * A[0][0][k] + sg[1] * A[1][0][k] - z_mp[k] for k in range(2)])
        if abs(w[0]) + abs(w[1]) < mp.mpf(10) ** (-(mp.mp.dps - 10)):
            return xs, ys
    return None


# ------------------------------------------------------------------ routed Abel-Jacobi (r0 -> hub -> x)
def _clear_seg(a, b, pts):
    d = b - a
    out = 1e300
    for z in pts:
        t = min(max(((z - a) * d.conjugate()).real / abs(d) ** 2, 0), 1)
        out = min(out, abs(a + t * d - z))
    return out


class Router:
    """Choose for each endpoint a path r0 -> x or r0 -> h -> x with clearance >= 5% of the span."""

    def __init__(self, F):
        self.F = F
        R = max(abs(F.roots - F.r0)) + 1.0
        self.hubs = [F.r0 + R * np.exp(2j * np.pi * k / 24) * rr for rr in (0.6, 1.1, 1.7) for k in range(24)]
        self.hubs = [h for h in self.hubs if _clear_seg(F.r0, h, F.rest) > 0.05 * F.span]

    def path(self, x):
        F = self.F
        rest_all = list(F.roots)
        rest_x = [z for z in rest_all]
        if _clear_seg(F.r0, x, F.rest) > 0.05 * F.span:
            return None
        best, bh = -1, None
        for h in self.hubs:
            c = min(_clear_seg(h, x, rest_x), _clear_seg(F.r0, h, F.rest))
            if c > best:
                best, bh = c, h
        if best < 0.03 * F.span:
            raise ValueError("no clear route")
        return bh


def _ysqrt_factory(lc, rest, a, b, sqrt_fn, conj=False):
    d = b - a
    mids = []
    for z in rest:
        t = ((z - a) * (d.conjugate() if hasattr(d, "conjugate") else mp.conj(d))).real / abs(d) ** 2
        t = min(max(t, 0), 1)
        mids.append(a + t * d)

    def ys(x):
        q = sqrt_fn(lc)
        for z, m in zip(rest, mids):
            q = q * sqrt_fn(m - z) * sqrt_fn((x - z) / (m - z))
        return q
    return ys


def fast_aj_routed(F, Rt, x):
    h = Rt.path(x)
    if h is None:
        return F.aj(x)
    A1, yh = F.aj(h)
    # leg h -> x, branch y = sqrt(f) continuous from yh
    ys = _ysqrt_factory(F.lc + 0j, list(F.roots), h, x, np.sqrt)
    eps = yh / ys(h)
    s = _S
    xx = h + (x - h) * s
    yy = eps * ys(xx)
    v0 = np.sum(_WS * (x - h) / yy)
    v1 = np.sum(_WS * xx * (x - h) / yy)
    return A1 + np.array([v0, v1]), eps * ys(x)


def mp_aj_routed(C, Rt, x):
    h = Rt.path(complex(x))
    if h is None:
        return C.aj(x)
    h = mp.mpc(h)
    A1, yh = C.aj(h)
    roots = list(C.roots)
    ys = _ysqrt_factory(C.lc, roots, h, x, mp.sqrt)
    eps = yh / ys(h)
    vals = []
    for m in (0, 1):
        vals.append(A1[m] + mp.quad(lambda s: (h + (x - h) * s) ** m * (x - h) / (eps * ys(h + (x - h) * s)), [0, 0.5, 1]))
    return vals, eps * ys(x)


def make(f, dps=40):
    C = Curve(f, dps=dps)
    F = FastCurve(C)
    F.router = Router(F)
    return C, F
