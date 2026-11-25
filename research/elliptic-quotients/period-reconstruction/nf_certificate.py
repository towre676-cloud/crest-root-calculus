"""Propose an algebraic coefficient field and verify the map identity in that field.

The normalized differential supplies a candidate primitive element.
Recognition can abstain if the chosen field or the numerical bounds are
insufficient. The accepted identity is checked over an exact number field.
"""
import cypari2
import mpmath as mp

pari = cypari2.Pari()
pari.allocatemem(128 * 10 ** 6)
pari("default(realprecision, 90)")


def _p(z, digits=60):
    return pari(f"{mp.nstr(mp.re(z), digits)} + I*({mp.nstr(mp.im(z), digits)})")


def minimal_polynomial(z, maxdeg=8, tol=mp.mpf(10) ** -40, bound=10 ** 25):
    zp = _p(z)
    if abs(mp.im(z)) < tol:
        q = pari.bestappr(pari(mp.nstr(mp.re(z), 60)), 10 ** 15)
        if abs(mp.mpf(str(pari.numerator(q))) / mp.mpf(str(pari.denominator(q))) - mp.re(z)) < tol:
            return pari("t") - q
    for d in range(2, maxdeg + 1):
        P = pari.algdep(zp, d)
        if pari.poldegree(P) != d:
            continue
        if max(abs(int(pari.polcoef(P, i))) for i in range(d + 1)) > bound:
            continue
        if abs(complex(pari.subst(P, "x", zp))) < float(tol) * 10 ** 5:
            # irreducibility
            fac = pari.factor(P)
            for i in range(len(fac[0])):
                Q = fac[0][i]
                if abs(complex(pari.subst(Q, "x", zp))) < 1e-30:
                    return pari.subst(Q, "x", pari("t")) / pari.pollead(Q)
    return None


def in_field(c, theta, deg, bound=10 ** 18):
    if deg == 1:
        if abs(mp.im(c)) > mp.mpf(10) ** -30:
            return None
        q = pari.bestappr(pari(mp.nstr(mp.re(c), 60)), bound)
        if abs(mp.mpf(str(pari.numerator(q))) / mp.mpf(str(pari.denominator(q))) - mp.re(c)) > mp.mpf(10) ** -30:
            return None
        return [q]
    vec = [c] + [theta ** k for k in range(deg)]
    rel = pari.lindep([_p(v) for v in vec], 45)
    rel = [int(rel[i]) for i in range(deg + 1)]
    if rel[0] == 0 or max(abs(r) for r in rel) > bound:
        return None
    coeffs = [pari(f"{-rel[k + 1]}/{rel[0]}") for k in range(deg)]
    val = sum(mp.mpf(str(pari.numerator(a))) / mp.mpf(str(pari.denominator(a))) * theta ** k for k, a in enumerate(coeffs))
    if abs(val - c) > mp.mpf(10) ** -30:
        return None
    return coeffs


def certify_nf(f, N, D, g2, g3, lam):
    ca, cb = lam
    gen = cb if abs(ca - 1) < mp.mpf(10) ** -40 else ca
    mpol = minimal_polynomial(gen) if abs(gen) > mp.mpf(10) ** -40 else pari("t")
    if mpol is None:
        return {"verified": False, "reason": "lambda not algebraic within bounds"}
    deg = int(pari.poldegree(mpol))
    theta = gen if deg > 1 else mp.mpf(1)
    rec = []
    for c in list(N) + list(D) + [g2, g3, ca, cb]:
        r = in_field(c, theta, deg)
        if r is None:
            return {"verified": False, "reason": f"coefficient not in Q(lambda), deg {deg}"}
        rec.append(r)
    mod = f"Mod(t, {mpol})" if deg > 1 else None

    def el(r):
        s = " + ".join(f"({a})*t^{k}" for k, a in enumerate(r))
        return pari(s.replace("t", mod)) if mod else pari(s.replace("*t^0", ""))

    d = len(N) - 1
    x = pari("x")
    Np = sum((el(rec[i]) * x ** i for i in range(d + 1)), pari(0))
    Dp = sum((el(rec[d + 1 + i]) * x ** i for i in range(d + 1)), pari(0))
    G2, G3, A, B = (el(r) for r in rec[-4:])
    fx = sum((int(c) * x ** i for i, c in enumerate(f)), pari(0))
    if Dp == 0 or A + B*x == 0 or pari.poldegree(pari.gcd(Np, Dp), "x") != 0:
        return {"verified": False, "reason": "degenerate or unreduced map"}
    W = pari.deriv(Np, "x") * Dp - Np * pari.deriv(Dp, "x")
    if W == 0:
        return {"verified": False, "reason": "constant horizontal map"}
    lhs = fx * W ** 2
    rhs = (A + B * x) ** 2 * Dp * (4 * Np ** 3 - G2 * Np * Dp ** 2 - G3 * Dp ** 3)
    ok = pari.simplify(pari.lift(lhs - rhs)) == 0 if mod else (lhs - rhs) == 0
    nonsing = pari.simplify(G2 ** 3 - 27 * G3 ** 2) != 0
    degree = max(int(pari.poldegree(pari.lift(Np) if mod else Np, "x")), int(pari.poldegree(pari.lift(Dp) if mod else Dp, "x")))
    return {"verified": bool(ok and nonsing), "field_polynomial": str(mpol), "N": str(Np), "D": str(Dp),
            "g2": str(G2), "g3": str(G3), "differential": f"({A}) dx/y + ({B}) x dx/y", "map_degree": degree}
