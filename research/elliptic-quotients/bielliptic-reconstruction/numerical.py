"""Equation-only reconstruction of the optimal elliptic spectrum of a bielliptic genus-two curve.

Input: integer coefficients of f, deg f in {5, 6}, for C : y^2 = f(x).  Output, for every
elliptic involution sigma defined over Q:

  * the involution, verified exactly as an automorphism of the branch locus;
  * j(E1), j(E2) for the two quotients C/sigma and C/(sigma iota), exactly, from binary-quartic
    invariants of the branch data on the quotient line;
  * a non-CM certificate (non-integral j);
  * the minimal cyclic isogeny degree k in {1} U {odd primes <= KMAX}, certified exactly by
    j1 = j2 or Phi_k(j1, j2) = 0 together with Phi_l(j1, j2) != 0 for smaller primes l;
  * the permutation pi = gamma^{-1} rho on E1[2], with the gluing gamma read off canonically
    from Weierstrass labels and rho read off the period lattices at high precision;
  * the complete optimal degree function d(a, b) = 2 (a^2 + k b^2) / |K_(a,b)|.

Epistemic labels: steps 1-3 are exact.  Step 4 is a high-precision numerical identification
of a discrete invariant, reported with its separation margin.
"""
from __future__ import annotations

import json
import sys
from fractions import Fraction
from math import gcd, isqrt, pi as PI

import cypari2

pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9)
pari.set_real_precision_bits(600)
pari("default(realprecision, 180)")
ODD_PRIMES = [3, 5, 7, 11, 13]
PREC = 640


def P(s):
    return pari(s)


def frac(g):
    return Fraction(int(pari.numerator(g)), int(pari.denominator(g)))


# ---------------------------------------------------------------- step 1: involutions
def pairings(items):
    if not items:
        yield []
        return
    a = items[0]
    for b in items[1:]:
        rest = [r for r in items[1:] if r != b]
        for tail in pairings(rest):
            yield [(a, b)] + tail


def find_involutions(f):
    """Rational Mobius involutions sigma(x) = (A x + B)/(C x - A) permuting the roots of f."""
    fx = pari.Pol(list(reversed(f)))
    roots = list(pari.polroots(fx, precision=PREC))
    if len(roots) == 5:
        roots.append(None)  # infinity
    found = []
    for Pr in pairings(list(range(6))):
        rows = []
        for i, j in Pr[:2]:
            a, b = roots[i], roots[j]
            if a is None:
                rows.append([1, 0, -b])          # sigma(inf) = b  <=> A = C b
            elif b is None:
                rows.append([1, 0, -a])
            else:
                rows.append([a + b, 1, -a * b])  # A(a+b) + B - C ab = 0
        A = rows[0][1] * rows[1][2] - rows[0][2] * rows[1][1]
        B = rows[0][2] * rows[1][0] - rows[0][0] * rows[1][2]
        C = rows[0][0] * rows[1][1] - rows[0][1] * rows[1][0]
        piv = max((A, B, C), key=lambda t: abs(t))
        A, B, C = A / piv, B / piv, C / piv
        cand = []
        for t in (A, B, C):
            q = pari.bestappr(pari.real(t), 10 ** 12)
            if abs(pari.real(t) - q) > 1e-60 or abs(pari.imag(t)) > 1e-60:
                cand = None
                break
            cand.append(frac(q))
        if cand is None:
            continue
        A, B, C = cand
        if A * A + B * C == 0:
            continue
        if verify_involution(f, A, B, C) and (A, B, C) not in found:
            found.append((A, B, C))
    # normalise projectively and deduplicate
    out = []
    for A, B, C in found:
        s = next(t for t in (C, A, B) if t != 0)
        key = (A / s, B / s, C / s)
        if key not in out:
            out.append(key)
    return out


def verify_involution(f, A, B, C):
    """Exact: (Cx - A)^6 f((Ax+B)/(Cx-A)) is a constant multiple of f (as sextic forms)."""
    x = P("x")
    num = P(f"({A.numerator}/{A.denominator})*x + ({B.numerator}/{B.denominator})")
    den = P(f"({C.numerator}/{C.denominator})*x - ({A.numerator}/{A.denominator})")
    g = sum((int(c) * num ** i * den ** (6 - i) for i, c in enumerate(f)), P(0))
    fx = sum((int(c) * x ** i for i, c in enumerate(f)), P(0))
    # compare as binary sextics: pad degree
    return pari.poldegree(g) <= 6 and all(
        pari.polcoef(g, i) * _lc6(f) == pari.polcoef(fx, i) * pari.polcoef(g, 6 if f[6] else _top(f)) or False
        for i in range(7)) or _proportional(g, fx)


def _lc6(f):
    return f[6] if f[6] else f[5]


def _top(f):
    return 6 if f[6] else 5


def _proportional(g, fx):
    cg = [pari.polcoef(g, i) for i in range(7)]
    cf = [pari.polcoef(fx, i) for i in range(7)]
    i0 = next(i for i in range(7) if cf[i] != 0)
    if cg[i0] == 0:
        return False
    lam = cg[i0] / cf[i0]
    return all(cg[i] == lam * cf[i] for i in range(7))


# ---------------------------------------------------------------- step 2: quotients and j
def quartic_j(q):
    """j-invariant of the genus-one double cover branched at the roots of a binary quartic."""
    a, b, c, d, e = [pari.polcoef(q, 4 - i) for i in range(5)]
    I = 12 * a * e - 3 * b * d + c * c
    J = 72 * a * c * e + 9 * b * c * d - 27 * a * d * d - 27 * e * b * b - 2 * c ** 3
    return 6912 * I ** 3 / (4 * I ** 3 - J ** 2)


def quotient_data(f, inv):
    A, B, C = inv
    q = lambda t: P(f"{t.numerator}/{t.denominator}")
    x, u = P("x"), P("u")
    sig_num, sig_den = q(A) * x + q(B), q(C) * x - q(A)
    # choose lambda rational, not a branch point, not fixed, with sigma(lambda) finite
    fx = sum((int(c) * x ** i for i, c in enumerate(f)), P(0))
    lam = None
    for t in range(2, 200):
        L = P(t)
        if pari.subst(fx, "x", L) == 0 or pari.subst(sig_den, "x", L) == 0:
            continue
        if pari.subst(sig_num, "x", L) / pari.subst(sig_den, "x", L) == L:
            continue
        lam = L
        break
    # phi = 1/(x - lam) + 1/(sigma(x) - lam) = N/D
    phi = 1 / (x - lam) + sig_den / (sig_num - lam * sig_den)
    N, D = pari.numerator(phi), pari.denominator(phi)
    # orbit values of the branch points (and of infinity if deg f = 5)
    F = pari.polresultant(fx, D * u - N, "x")
    if pari.poldegree(fx) == 5:  # root at infinity: phi(inf) = 1/(A/C - lam) or 0
        val = 0 if C == 0 else 1 / (q(A) / q(C) - lam)
        F = F * (u - val)
    F = pari.subst(F, "u", P("u"))
    sq = pari.gcd(F, pari.deriv(F, "u"))
    cubic = F / sq
    cubic = cubic / pari.pollead(cubic)
    assert pari.poldegree(cubic) == 3, "orbit polynomial must be a cubic"
    # fixed points of sigma: C x^2 - 2 A x - B = 0 (one at infinity if C == 0)
    fix = q(C) * x ** 2 - 2 * q(A) * x - q(B)
    Pfix = pari.polresultant(fix, D * u - N, "x")
    if C == 0:
        Pfix = Pfix * (u - 0)
    Pfix = Pfix / pari.pollead(Pfix)
    assert pari.poldegree(Pfix) == 2
    return cubic, Pfix


def j_pair(cubic, Pfix):
    fac = pari.factor(Pfix)
    u = P("u")
    if len(fac[0]) == 2 or fac[1][0] == 2:
        rts = [-pari.polcoef(fac[0][i], 0) / pari.polcoef(fac[0][i], 1) for i in range(len(fac[0]))]
        js = [quartic_j(cubic * (u - r)) for r in rts]
        field = None
    else:
        t = P("Mod(t, %s)" % pari.subst(Pfix, "u", P("t")))
        j1 = quartic_j(cubic * (u - t))
        js = [j1, pari.subst(pari.lift(j1), "t", P("t"))]  # conjugate computed below
        js = [j1, P("Mod(%s, %s)" % (pari.subst(pari.lift(j1), "t", -P("t") - pari.polcoef(Pfix, 1)), pari.subst(Pfix, "u", P("t"))))]
        field = str(pari.subst(Pfix, "u", P("t")))
    return js, field


def is_integral(j):
    if pari.type(j) == "t_POLMOD" and pari.poldegree(pari.lift(j)) <= 0:
        j = pari.polcoef(pari.lift(j), 0)
    if pari.type(j) == "t_POLMOD":
        cp = pari.charpoly(j)
        return all(pari.denominator(pari.polcoef(cp, i)) == 1 for i in range(pari.poldegree(cp) + 1))
    return pari.denominator(j) == 1


def isogeny_degree(j1, j2):
    """Minimal k in {1} U ODD_PRIMES with a cyclic k-isogeny, exact; also reports l = 2."""
    if j1 == j2:
        return 1, {}
    tested = {}
    for l in [2] + ODD_PRIMES:
        Phi = pari.polmodular(l)
        val = pari.subst(pari.subst(Phi, "x", j1), "y", j2)
        tested[l] = (val == 0)
        if val == 0:
            return l, tested
    return None, tested


# ---------------------------------------------------------------- step 4: rho on 2-torsion
def cubic_model_roots(cubic, Pval):
    """2-torsion abscissae e_i = 1/(u_i - P) of the model with P sent to infinity, labelled by u_i."""
    us = list(pari.polroots(cubic, precision=PREC))
    return us, [1 / (ui - Pval) for ui in us]


def ell_from_roots(es):
    s1 = es[0] + es[1] + es[2]
    s2 = es[0] * es[1] + es[0] * es[2] + es[1] * es[2]
    s3 = es[0] * es[1] * es[2]
    return pari.ellinit([0, -s1, 0, s2, -s3], precision=PREC)


def gauss_reduce(w1, w2):
    while True:
        if abs(w1) > abs(w2):
            w1, w2 = w2, w1
        m = pari.round(pari.real(w2 / w1))
        if m == 0:
            return w1, w2
        w2 = w2 - m * w1


def _coords(z, b1, b2):
    """Real coordinates of z in the basis (b1, b2)."""
    det = pari.imag(pari.conj(b1) * b2)
    s = pari.imag(pari.conj(z) * b2) / det
    t = pari.imag(pari.conj(b1) * z) / det
    return s, t


def _homothety(w1, w2, l1, l2, tol):
    """All c with c<w1, w2> = <l1, l2> (as lattices), found among short-vector ratios."""
    W = gauss_reduce(w1, w2)
    L = gauss_reduce(l1, l2)
    shorts = [L[0], L[1], L[0] + L[1], L[0] - L[1]]
    out = []
    for v in shorts:
        for sg in (1, -1):
            c = sg * v / W[0]
            ok = True
            M = []
            for w in (w1, w2):
                s, t = _coords(c * w, l1, l2)
                rs, rt = pari.round(s), pari.round(t)
                if abs(s - rs) > tol or abs(t - rt) > tol:
                    ok = False
                    break
                M.append((int(rs), int(rt)))
            if ok and abs(M[0][0] * M[1][1] - M[0][1] * M[1][0]) == 1:
                out.append(c)
    return out


def rho_permutation(cubic, Pfix, k):
    Ps = list(pari.polroots(Pfix, precision=PREC))
    us = list(pari.polroots(cubic, precision=PREC))
    tol = pari("1e-100")
    results = []
    for P1, P2 in [(Ps[0], Ps[1])]:
        e = [1 / (ui - P1) for ui in us]
        f2 = [1 / (ui - P2) for ui in us]
        E1, E2 = ell_from_roots(e), ell_from_roots(f2)
        w1, w2 = E1.omega()[0], E1.omega()[1]
        h1, h2 = E2.omega()[0], E2.omega()[1]
        for a in range(1, k + 1):
            if k % a:
                continue
            d = k // a
            for b in range(d):
                l1, l2 = a * h1 + b * h2, d * h2
                for c in _homothety(w1, w2, l1, l2, tol):
                    perm, margin = [None] * 3, 1e300
                    for z in (w1 / 2, w2 / 2, (w1 + w2) / 2):
                        x1 = pari.ellztopoint(E1, z, precision=PREC)[0]
                        d1 = sorted((abs(x1 - e[t]), t) for t in range(3))
                        x2 = pari.ellztopoint(E2, c * z, precision=PREC)[0]
                        d2 = sorted((abs(x2 - f2[t]), t) for t in range(3))
                        margin = min(margin, float(d1[1][0] / max(d1[0][0], pari("1e-170"))),
                                     float(d2[1][0] / max(d2[0][0], pari("1e-170"))))
                        perm[d1[0][1]] = d2[0][1]
                    if None not in perm and len(set(perm)) == 3:
                        results.append((perm, margin, (a, b, d)))
    if not results:
        return None
    perms = {tuple(r[0]) for r in results}
    if len(perms) != 1:
        raise RuntimeError(f"inconsistent 2-torsion actions {perms}")
    best = max(results, key=lambda r: r[1])
    return list(best[0]), best[1], len(results)


# ---------------------------------------------------------------- step 5: spectrum
PERM_MATRIX = {"identity": ((1, 0), (0, 1)), "transposition": ((0, 1), (1, 0)), "3-cycle": ((0, 1), (1, 1))}


def kernel_order(a, b, k, G):
    """|{P in E1 : (aP, b rho P) in graph(gamma)}|, rho = diag(1, k), gamma = G on 2-torsion (k odd)."""
    g = gcd(a, k) if a else k
    N = 2 * g
    count = 0
    for x in range(N):
        for y in range(N):
            if (2 * a * x) % N or (2 * a * y) % N or (2 * b * x) % N or (2 * b * k * y) % N:
                continue
            v1 = ((2 * a * x // N) % 2, (2 * a * y // N) % 2)
            v2 = ((2 * b * x // N) % 2, (2 * b * k * y // N) % 2)
            Gv = ((G[0][0] * v1[0] + G[0][1] * v1[1]) % 2, (G[1][0] * v1[0] + G[1][1] * v1[1]) % 2)
            if Gv == v2:
                count += 1
    return count


def degree(a, b, k, G):
    n = 2 * (a * a + k * b * b)
    K = kernel_order(a, b, k, G)
    assert n % K == 0
    return n // K


def spectrum(k, G, bound):
    """Exact multiplicities of optimal degrees <= bound.

    |K| divides 4 gcd(a, k)^2 <= 4k^2 and d = 2(a^2 + k b^2)/|K|, so d <= bound forces
    a^2 <= 2 k^2 bound and b^2 <= 2 k bound; the window below is therefore complete.
    """
    from collections import Counter
    cnt = Counter()
    Ra = isqrt(2 * k * k * bound) + 1
    Rb = isqrt(2 * k * bound) + 1
    for a in range(0, Ra + 1):
        for b in range(-Rb, Rb + 1):
            if a == 0 and b != 1:
                continue
            if gcd(a, b) != 1:
                continue
            g = gcd(a, k)
            if 2 * (a * a + k * b * b) > bound * 4 * g * g:
                continue
            d = degree(a, b, k, G)
            if d <= bound:
                cnt[d] += 1
    return cnt


def leading_constant(k, G):
    """c with N(X) ~ c X, from the PR #17 counting theorem (r = 2, n = 2k)."""
    n = 2 * k
    total = 0
    for a in range(n):
        for b in range(n):
            if gcd(gcd(a, b), n) != 1:
                continue
            total += kernel_order(a, b, k, G)
    euler = 1.0
    for l in set(_pf(n)):
        euler *= 1 - l ** -2
    area = PI / (2 * k ** 0.5)
    zeta2 = PI ** 2 / 6
    return area * total / (2 * n * n * zeta2 * euler), total


def _pf(n):
    out, d = [], 2
    while d * d <= n:
        while n % d == 0:
            out.append(d)
            n //= d
        d += 1
    if n > 1:
        out.append(n)
    return out


# ---------------------------------------------------------------- driver
def reconstruct(f, bound=200):
    invs = find_involutions(f)
    report = {"f": f, "involutions": [[str(t) for t in inv] for inv in invs], "quotients": []}
    for inv in invs:
        cubic, Pfix = quotient_data(f, inv)
        (j1, j2), field = j_pair(cubic, Pfix)
        rec = {"involution": [str(t) for t in inv], "orbit_cubic": str(cubic), "fixed_image_poly": str(Pfix),
               "j1": str(j1), "j2": str(j2), "field": field,
               "non_cm_certificate": (not is_integral(j1)) and (not is_integral(j2))}
        k, tested = isogeny_degree(j1, j2)
        rec["k"] = k
        rec["modular_tests"] = {str(l): bool(v) for l, v in tested.items()}
        if k is None:
            rec["status"] = "no cyclic isogeny of degree 1 or prime <= 13: multiplicity one or out of range"
        elif k == 2:
            rec["status"] = "k = 2 is outside the odd-k theorem"
        else:
            pm, margin, nhom = rho_permutation(cubic, Pfix, k)
            ctype = cycle_name(pm)
            G = PERM_MATRIX[ctype]
            rec.update({"rho_on_labels": pm, "pi_type": ctype, "identification_margin": margin,
                        "homotheties_found": nhom})
            if ctype == "identity" and k == 1:
                rec["status"] = "inconsistent: identity gluing with k = 1 is not a Jacobian"
            else:
                sp = spectrum(k, G, bound)
                c, tot = leading_constant(k, G)
                rec.update({"spectrum": {str(d): m for d, m in sorted(sp.items())}, "leading_constant": c,
                            "residue_kernel_sum": tot, "status": "reconstructed"})
        report["quotients"].append(rec)
    return report


def cycle_name(pm):
    fixed = sum(1 for i in range(3) if pm[i] == i)
    return {3: "identity", 1: "transposition", 0: "3-cycle"}[fixed]


if __name__ == "__main__":
    f = json.loads(sys.argv[1])
    print(json.dumps(reconstruct(f), indent=1, default=str))
