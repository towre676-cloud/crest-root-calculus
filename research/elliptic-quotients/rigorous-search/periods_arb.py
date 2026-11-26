"""Rigorous period matrix of a genus-two curve in ball arithmetic (Arb via python-flint).

Same homology as period-reconstruction/periods.py: after the exact Mobius change
x = z0 + 1/xt (z0 a Gaussian rational), a branch point e1 and the star of segments to the other
roots give loops gamma_j with gamma_j . gamma_k = +1 for j before k in angular order from a cut,
and periods 2 int_{e1}^{ej} omega_+.  Here every quantity is a ball:

  * roots: rigorous isolation by acb_poly.roots;
  * integrals: acb.integral with the integrand declared analytic only where every square root
    stays off its branch cut (acb.sqrt(analytic=True) returns a non-finite ball otherwise), so a
    wrong branch can make the computation fail but cannot produce a wrong enclosure;
  * the discrete choices (base root, cut, angular order, the global sign of y_+) are made from
    midpoints and are accepted only after strict ball comparisons or disjointness tests.

The output matrix of (dx/y, x dx/y) periods on the original model is an enclosure of the true
period matrix for the stated homology basis, with intersection form E (exact integers).
"""
from __future__ import annotations

from fractions import Fraction

from flint import acb, acb_mat, acb_poly, arb, ctx


def set_prec(bits):
    ctx.prec = bits


def gauss(re, im):
    return acb(arb(re), arb(im))


Z0_CHOICES = [(Fraction(a, 7), Fraction(b, 11)) for a in (-5, -2, 1, 3, 6) for b in (2, 3, 5, 8, 13)]


def moved_poly(f, z0):
    """Enclose the exact Gaussian-rational coefficients of x^6 f(z0 + 1/x)."""
    from math import comb
    zr, zi = z0
    z = acb(arb(zr.numerator) / zr.denominator, arb(zi.numerator) / zi.denominator)
    g = [acb(0)] * 7
    for i, c in enumerate(f):
        for k in range(i + 1):
            g[k] += acb(c) * comb(i, k) * z ** (i - k)
    out = [acb(0)] * 7
    for k in range(7):
        out[6 - k] = g[k]
    return out, z


def mid(z):
    return complex(float(z.real.mid()), float(z.imag.mid()))


def nearest(a, b, p):
    d = b - a
    t = ((p - a) * d.conjugate()).real / abs(d) ** 2
    t = min(max(t, 0.0), 1.0)
    return a + t * d


def clearance(e, i):
    m = float("inf")
    for j in range(6):
        if j == i:
            continue
        for k in range(6):
            if k in (i, j):
                continue
            m = min(m, abs(nearest(e[i], e[j], e[k]) - e[k]))
    return m


def star_periods_arb(fco):
    """fco: list of 7 acb coefficients (low to high) of the moved sextic."""
    P = acb_poly(fco)
    roots = P.roots(tol=arb(2) ** (-(ctx.prec - 20)), maxprec=8 * ctx.prec)
    return star_periods_from_roots(roots, fco[6])


def star_periods_from_roots(roots, lc):
    """Periods of the sextic lc times the product over six supplied branch points."""
    if lc.contains(0) or len(roots) != 6 or any(
        x.overlaps(y) for i, x in enumerate(roots) for y in roots[:i]
    ):
        raise ArithmeticError("The moved sextic does not have six isolated root balls.")
    em = [mid(r) for r in roots]
    import cmath
    best = None
    for i in range(6):
        angs = sorted(cmath.phase(em[j] - em[i]) for j in range(6) if j != i)
        gaps = [((angs[(t + 1) % 5] - angs[t]) % (2 * cmath.pi)) for t in range(5)]
        g = max(gaps)
        cl = clearance(em, i)
        if best is None or cl > best[0]:
            t = gaps.index(g)
            best = (cl, i, angs[t] + g / 2, g)
    _, i1, cut, _ = best
    e1, e1m = roots[i1], em[i1]
    others = [(roots[j], em[j]) for j in range(6) if j != i1]

    def theta(zm):
        th = cmath.phase(zm - e1m)
        while th <= cut:
            th += 2 * cmath.pi
        while th > cut + 2 * cmath.pi:
            th -= 2 * cmath.pi
        return th

    others.sort(key=lambda t: theta(t[1]))
    # Rotate the proposed cut to the negative real axis using an exact dyadic.
    r = -cmath.rect(1, -cut)
    rotation = acb(r.real, r.imag)
    rays = [rotation * (z - e1) for z, _ in others]
    if any(w.imag.contains(0) and not w.real > 0 for w in rays):
        raise ArithmeticError("A ray might meet the angular cut.")
    angles = [w.arg() for w in rays]
    if not all(a < b for a, b in zip(angles, angles[1:])):
        raise ArithmeticError("Cannot certify the cyclic order of the rays.")
    common = (1 / rotation).sqrt()
    fpe1 = lc
    for z, _ in others:
        fpe1 = fpe1 * (e1 - z)
    rho0 = fpe1.sqrt()
    periods = []
    pi = arb.pi()
    for (ej, ejm), ray in zip(others, rays):
        sq = ray.sqrt(analytic=True) * common
        rest = [(z, zm) for z, zm in others if z is not ej]
        mids = [nearest(e1m, ejm, zm) for z, zm in rest]
        midsb = [acb(m.real, m.imag) for m in mids]
        for (z, _), m in zip(rest, midsb):
            if (m - z).contains(0) or not all(
                ((end - z) / (m - z)).real > 0 for end in (e1, ej)
            ):
                raise ArithmeticError("Cannot certify a square-root branch along the segment.")
        pref = (-(ej - e1) * lc).sqrt()

        def Qraw(t, analytic):
            x = e1 + (ej - e1) * (t / 2).sin() ** 2
            q = pref
            for (z, _), m in zip(rest, midsb):
                q = q * (m - z).sqrt() * ((x - z) / (m - z)).sqrt(analytic=analytic)
            return q
        q0 = Qraw(acb(0), False)
        positive, negative = q0.overlaps(rho0), q0.overlaps(-rho0)
        if positive == negative:
            raise ArithmeticError("Cannot determine the relative signs of the lifted segments.")
        sign = 1 if positive else -1
        row = []
        for m in (0, 1):
            def integrand(t, analytic, m=m):
                angle = pi * t
                x = e1 + (ej - e1) * (angle / 2).sin() ** 2
                return pi * x ** m * sq / (sign * Qraw(angle, analytic))
            val = acb(0)
            for k in range(8):
                a, b = arb(k) / 8, arb(k + 1) / 8
                val += acb.integral(integrand, acb(a), acb(b), rel_tol=arb(2) ** (-(ctx.prec - 20)),
                                    abs_tol=arb(2) ** (-(ctx.prec - 20)), deg_limit=4 * ctx.prec,
                                    eval_limit=10 ** 7, depth_limit=4 * ctx.prec)
            if not val.is_finite():
                raise ArithmeticError("The period integral was not enclosed.")
            row.append(2 * val)
        periods.append(row)
    I = [[0 if j == k else (1 if j < k else -1) for k in range(5)] for j in range(5)]
    relation = [periods[0][m] - periods[1][m] + periods[2][m] - periods[3][m] + periods[4][m] for m in range(2)]
    if not all(z.contains(0) for z in relation):
        raise ArithmeticError("The alternating homology relation failed.")
    return {"periods": periods, "intersection": I}


def period_matrix_arb(f, bits=300):
    """Rigorous 2x4 period matrix of (dx/y, x dx/y) on y^2 = f(x), basis gamma_1..gamma_4, form E."""
    if len(f) not in (6, 7) or not f[-1]:
        raise ValueError("Use a degree-five or degree-six model.")
    set_prec(bits)
    best = None
    for z0 in Z0_CHOICES:
        g, z = moved_poly(f, z0)
        if g[6].contains(0):
            continue
        em = [mid(r) for r in acb_poly(g).roots(tol=arb(2) ** -60, maxprec=4000)]
        span = max(abs(a - b) for a in em for b in em)
        cl = max(clearance(em, i) for i in range(6)) / span
        if best is None or cl > best[0]:
            best = (cl, z0)
    if best is None:
        raise ArithmeticError("No admissible Möbius center was found.")
    g, z = moved_poly(f, best[1])
    S = star_periods_arb(g)
    P = []
    for w0, w1 in S["periods"]:
        P.append([-w1, -(z * w1 + w0)])     # dx/y = -xt dxt/yt,  x dx/y = -(z0 xt + 1) dxt/yt
    S["z0"] = best[1]
    S["Pi"] = [[P[c][r] for c in range(4)] for r in range(2)]
    S["E"] = [row[:4] for row in S["intersection"][:4]]
    P, E = acb_mat(S["Pi"]), acb_mat(S["E"])
    riemann = P * E.inv() * P.transpose()
    if not all(riemann[i, j].contains(0) for i in range(2) for j in range(2)):
        raise ArithmeticError("Riemann’s first bilinear relation failed.")
    return S
