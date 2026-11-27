"""Optimal degree-four elliptic quotients defined over number fields.

The rational solver in etale_fibres.py keeps only Galois-stable pairings of the
Weierstrass points and only rational points of the pole incidence. This module
drops both restrictions. A pairing is a fixed-point-free involution of the
sextic algebra, and it is located exactly over its field of definition K without
computing a splitting field. The incidence is then solved over K, and every
algebraic solution is returned over its absolute point field L, one
representative for each Galois orbit. Each cover is checked by the exact
identity S^6 g(N/S^2) = f H^2 in L[x].
"""
from __future__ import annotations

from dataclasses import dataclass

import sympy as sp

from reconstruct_four import pari, sextic_chart, x
from etale_fibres import pair_indices

pari.allocatemem(2*10**9)
for _name in ("x", "u", "t", "a", "b", "z"):
    pari(_name)
# Field generators must have the lowest variable priority.
W = pari('W=varlower("w")')
X, U, T, A, B, Z = (pari(v) for v in ("x", "u", "t", "a", "b", "z"))


# ----------------------------------------------------------------------------
# Number fields as absolute polynomials in the lowest-priority variable w.

@dataclass(frozen=True)
class Field:
    modulus: object  # irreducible polynomial in w, or None for Q

    @property
    def degree(self):
        return 1 if self.modulus is None else int(pari.poldegree(self.modulus))


Q = Field(None)


def to_w(poly_in_z):
    return pari.substpol(poly_in_z, Z, W)


def extend(field, factor):
    """Adjoin a root of an irreducible factor in z over the field.

    Returns the absolute field, the image of the old generator and the new root.
    """
    if pari.poldegree(factor, Z) == 1:
        root = -pari.polcoef(factor, 0, Z)/pari.polcoef(factor, 1, Z)
        return field, (None if field.modulus is None else pari.Mod(W, field.modulus)), root
    if field.modulus is None:
        absolute = pari.polredbest(to_w(factor), 1)
        poly, root = absolute[0], absolute[1]
        return Field(poly), None, pari.Mod(pari.lift(root), poly)
    relative = pari.rnfequation(field.modulus, pari.lift(factor), 1)
    poly, old, k = to_w(relative[0]), pari.lift(relative[1]), relative[2]
    reduced = pari.polredbest(poly, 1)
    new_poly, change = reduced[0], pari.lift(reduced[1])
    old = pari.Mod(pari.subst(to_w(old), W, change), new_poly)
    theta = pari.Mod(change, new_poly)
    return Field(new_poly), old, theta - k*old


def reembed(value, field, image, new_field):
    if field.modulus is None or image is None:
        return value if new_field.modulus is None else value*pari.Mod(1, new_field.modulus)
    return pari.subst(pari.lift(value), W, pari.lift(image))*pari.Mod(1, new_field.modulus)


def is_zero(value):
    return pari.lift(value) == 0


def factor_over(field, poly):
    """Irreducible factors of a univariate polynomial in z over the field."""
    if pari.poldegree(poly, Z) <= 0:
        return []
    if field.modulus is None:
        factors = pari.factor(poly)[0]
    else:
        factors = pari.nffactor(field.modulus, poly)[0]
    return [f for f in factors if pari.poldegree(f, Z) > 0]


# ----------------------------------------------------------------------------
# Pairings as fixed-point-free involutions of A = Q[x]/(f).

@dataclass(frozen=True)
class NFPairing:
    field: Field
    involution: object  # sigma(x) in K[x] mod f
    trace: object       # x + sigma(x)
    norm: object        # x sigma(x)
    orbit: int          # number of conjugate pairings
    resolvent: object   # minimal polynomial of the separating invariant


def pari_source(source):
    return pari(str(source.as_expr()).replace("**", "^"))


def compose(poly, inner, f):
    """poly(inner) mod f by Horner's rule, keeping field coefficients reduced."""
    result = 0*X
    for k in range(int(pari.poldegree(poly, X)), -1, -1):
        result = (result*inner + pari.polcoef(poly, k, X)) % f
    return result


def involution_ok(f, sigma):
    if not is_zero(compose(f, sigma, f)):
        return False
    if not is_zero(compose(sigma, sigma, f) - X):
        return False
    return pari.poldegree(pari.gcd(X - sigma, f), X) == 0


WEIGHTS = ([1, 2, 3, 5, 7, 11], [1, -1, 4, -9, 16, -25], [3, 1, 4, 1, 5, 9], [2, -7, 1, 8, -2, 8])


def functional_charpoly(sigma, weights):
    """Characteristic polynomial over Q of a fixed rational functional of the coefficients."""
    element = sum(w*pari.polcoef(sigma, k, X) for k, w in enumerate(weights))
    element = pari(element)
    if element.type() != "t_POLMOD":
        return Z - element
    return pari.substpol(pari.charpoly(element), pari.variable(pari.charpoly(element)), Z)


def pairwise_distinct(pairings):
    """Exact proof that the conjugates of all certified involutions are pairwise distinct.

    Equal involutions over the algebraic closure have equal values of every functional.
    If the product of the characteristic polynomials is squarefree, no two of the
    sum-of-degrees conjugates coincide, within an orbit or across orbits.
    """
    for weights in WEIGHTS:
        product = 1 + 0*Z
        for p in pairings:
            product *= functional_charpoly(p.involution, weights)
        if pari.poldegree(pari.gcd(product, pari.deriv(product, Z)), Z) == 0:
            return list(weights)
    return None


def numerical_pairings(f, digits):
    """All 15 pairings as numerical involution polynomials, with an integral invariant."""
    pari.set_real_precision(digits)
    # Roots of a primitive integral model times its leading coefficient are algebraic integers.
    lead = pari.pollead(f/pari.content(f))
    roots = list(pari('polroots(%s)' % f))
    integral = [lead*r for r in roots]
    derivative = pari.deriv(f, X)
    idempotents = [pari.divrem(f, X - r, X)[0]/pari.subst(derivative, X, r) for r in roots]
    data = []
    for pairs in pair_indices(tuple(range(6))):
        partner = {}
        for i, j in pairs:
            partner[i], partner[j] = j, i
        sigma = sum(roots[partner[i]]*idempotents[i] for i in range(6))
        data.append((pairs, sigma))
    return roots, integral, data


def pairing_orbits(source, digits=400):
    """Galois orbits of pairings, each with an exact involution over its field of definition."""
    f = pari_source(source)
    for attempt in range(6):
        roots, integral, data = numerical_pairings(f, digits)
        # Integral symmetric invariants of the pairs; a power sum that separates all 15 pairings.
        for weight in range(1, 60):
            invariants = [sum((integral[i]*integral[j] + weight*(integral[i] + integral[j]))**3
                              for i, j in pairs) for pairs, _ in data]
            if min(abs(p - q) for k, p in enumerate(invariants) for q in invariants[k+1:]) > 1e-20:
                break
        else:
            raise ArithmeticError("No separating invariant in the tested range.")
        resolvent = 1 + 0*Z
        for v in invariants:
            resolvent *= Z - v
        coeffs = [pari.polcoef(resolvent, k, Z) for k in range(16)]
        rounded = [pari.round(pari.real(c)) for c in coeffs]
        if max(abs(c - r) for c, r in zip(coeffs, rounded)) > pari("1e-50"):
            digits *= 2
            continue
        resolvent = sum(r*Z**k for k, r in enumerate(rounded))
        if pari.poldegree(pari.gcd(resolvent, pari.deriv(resolvent, Z)), Z) > 0:
            digits *= 2
            continue
        orbits, ok = [], True
        for factor in pari.factor(resolvent)[0]:
            members = [k for k, v in enumerate(invariants) if abs(pari.subst(factor, Z, v)) <
                       pari("1e-30")*(1+abs(v))**pari.poldegree(factor, Z)]
            if len(members) != pari.poldegree(factor, Z):
                ok = False
                break
            pairing = exact_pairing(f, factor, [invariants[k] for k in members],
                                    [data[k][1] for k in members], digits)
            if pairing is None:
                ok = False
                break
            orbits.append(pairing)
        if ok and sum(p.orbit for p in orbits) == 15 and pairwise_distinct(orbits) is not None:
            return f, orbits
        digits *= 2
    raise ArithmeticError("Pairing recognition failed at every precision.")


def exact_pairing(f, factor, values, sigmas, digits):
    degree = int(pari.poldegree(factor, Z))
    if degree == 1:
        field = Q
        sigma = sum(pari.bestappr(pari.real(pari.polcoef(sigmas[0], k, X)), pari(10)**(digits//3))*X**k
                    for k in range(6))
    else:
        absolute = pari.polredbest(to_w(factor), 1)
        modulus, root = absolute[0], pari.lift(absolute[1])
        field = Field(modulus)
        # Interpolate in the separating invariant, then rewrite in the reduced generator.
        sigma = 0*X
        for k in range(6):
            interpolant = pari.polinterpolate(values, [pari.polcoef(s, k, X) for s in sigmas], W)
            exact = sum(pari.bestappr(pari.real(pari.polcoef(interpolant, i, W)), pari(10)**(digits//3))*W**i
                        for i in range(degree))
            sigma += pari.Mod(pari.subst(exact, W, root), modulus)*X**k
    if not involution_ok(f, sigma):
        return None
    norm = (X*sigma) % f
    return NFPairing(field, sigma, (X + sigma) % f, norm, degree, factor)


# ----------------------------------------------------------------------------
# The incidence over K, following map_conditions in etale_fibres.py.

def coordinates(vectors, target, field):
    """Solve sum c_i vectors[i] = target for c in K."""
    matrix = pari.matrix(6, len(vectors), [pari.polcoef(v, r, X) for r in range(6) for v in vectors])
    rhs = pari.Col([pari.polcoef(target, r, X) for r in range(6)])
    solution = pari.matinverseimage(matrix, rhs)
    if len(solution) == 0:
        raise ArithmeticError("Element outside the paired subalgebra.")
    return [solution[i] for i in range(len(vectors))]


@dataclass(frozen=True)
class Algebra:
    generator: object
    modulus: object   # cubic in u over K
    trace: object     # trace function in u
    norm: object


def paired_algebra(f, pairing):
    one = 1 + 0*X
    for k in range(6):
        generator = (pairing.trace + k*pairing.norm) % f
        powers = [one, generator, (generator*generator) % f]
        matrix = pari.matrix(6, 3, [pari.polcoef(p, r, X) for r in range(6) for p in powers])
        if pari.matrank(matrix) == 3:
            break
    else:
        raise ArithmeticError("No generator separates the three pairs.")
    cube = (powers[2]*generator) % f
    c = coordinates(powers, cube, pairing.field)
    modulus = U**3 - (c[2]*U**2 + c[1]*U + c[0])
    tr = coordinates(powers, pairing.trace, pairing.field)
    nm = coordinates(powers, pairing.norm, pairing.field)
    return Algebra(generator, modulus, tr[0] + tr[1]*U + tr[2]*U**2, nm[0] + nm[1]*U + nm[2]*U**2)


def conditions(f, algebra, degree):
    s = [1 + 0*X, X + A, X**2 + A*X + B][degree]
    d = s*s
    n1 = (d*algebra.generator) % f
    n2 = (d*algebra.generator**2) % f
    c1, c2 = pari.polcoef(n1, 5, X), pari.polcoef(n2, 5, X)
    if is_zero(c1) and is_zero(c2):
        return None
    n = c2*n1 - c1*n2
    value = c2*U - c1*U**2
    q = X**2 - algebra.trace*X + algebra.norm
    quotient, remainder = pari.divrem(n - value*d, q, X)
    p0, p1, p2 = (pari.polcoef(quotient, k, X) for k in range(3))
    disc = pari.divrem(p1**2 - 4*p0*p2, algebra.modulus, U)[1]
    equations = [pari.polcoef(disc, k, U) for k in range(3)]
    charpoly = pari.polresultant(algebra.modulus, T - value, U)
    # The degeneracy locus as a list of factors, each of which must stay nonzero.
    unit = [pari.polresultant(s, f, X), pari.polresultant(s, n, X), pari.poldisc(charpoly, T)]
    if degree < 2:
        unit.append(pari.polcoef(n, 4, X))
    return s, n, value, [e for e in equations if not is_zero(e)], unit


def strip(poly, unit, variables):
    """Remove from poly every factor shared with a degeneracy factor (the saturation)."""
    for factor in unit:
        while True:
            g = pari.gcd(poly, factor)
            main = next((v for v in variables if pari.poldegree(g, v) > 0), None)
            if main is None:
                break
            poly = pari.divrem(poly, g, main)[0]
    return poly


def degenerate(unit, field, image, new, variables, point):
    return any(is_zero(reembed_eval(u, field, image, new, variables, point)) for u in unit)


def evaluate(value, variables, point):
    for v, p in zip(variables, point):
        value = pari.subst(value, v, p)
    return value


def solve_incidence(field, equations, unit, variables):
    """Algebraic points of V(equations) outside V(unit), one per Galois orbit over K."""
    if any(is_zero(u) for u in unit):
        return []  # the whole chart lies in the degeneracy locus
    if not variables:
        return [(field, None, ())] if not equations else []
    if not equations:
        raise ArithmeticError("Positive-dimensional incidence.")
    common = equations[0]
    for e in equations[1:]:
        common = pari.gcd(common, e)
    if len(variables) == 1:
        # A univariate common divisor is the finite solution set itself, not a family.
        (a,) = variables
        saturated = strip(common, unit, variables)
        points = []
        for factor in factor_over(field, pari.substpol(saturated, a, Z)):
            new, image, root = extend(field, factor)
            point = (root,)
            if any(not is_zero(reembed_eval(e, field, image, new, variables, point)) for e in equations):
                continue
            if not degenerate(unit, field, image, new, variables, point):
                points.append((new, image, point))
        return points
    residual = strip(common, unit, variables)
    if any(pari.poldegree(residual, v) > 0 for v in variables):
        raise ArithmeticError("Positive-dimensional incidence outside the degeneracy locus.")
    reduced = [pari.divrem(e, common, variables[0])[0] if pari.poldegree(common, variables[0]) > 0
               else (pari.divrem(e, common, variables[-1])[0] if pari.poldegree(common, variables[-1]) > 0 else e)
               for e in equations]
    points = []
    a, b = variables
    # Every equation cuts the projection, so spurious pairwise intersections drop out.
    first, projection = reduced[0], 0*A
    for other in reduced[1:]:
        projection = pari.gcd(projection, pari.polresultant(first, other, b))
    if is_zero(projection):
        raise ArithmeticError("Projection vanished identically.")
    for factor in factor_over(field, pari.substpol(projection, a, Z)):
        # Over the relative extension K[z]/(factor), discard fibres lying in the degeneracy locus
        # before any absolute field is formed.
        one = pari.Mod(1, factor)
        fibre = 0*B
        for e in reduced:
            fibre = pari.gcd(fibre, pari.substpol(e, a, Z)*one)
        fibre = strip(fibre, [pari.substpol(u, a, Z)*one for u in unit], [b])
        if pari.poldegree(fibre, b) <= 0:
            continue
        mid, image, alpha = extend(field, factor)
        specialised = [pari.subst(reembed(e, field, image, mid), a, alpha) for e in reduced]
        g = 0*B
        for e in specialised:
            g = pari.gcd(g, e)
        if pari.poldegree(g, b) <= 0:
            continue
        for second in factor_over(mid, pari.substpol(g, b, Z)):
            top, image2, beta = extend(mid, second)
            alpha2 = reembed(alpha, mid, image2, top)
            full = compose_images(field, image, mid, image2, top)
            point = (alpha2, beta)
            if any(not is_zero(reembed_eval(e, field, full, top, variables, point)) for e in equations):
                continue
            if degenerate(unit, field, full, top, variables, point):
                continue
            points.append((top, full, point))
    return points


def compose_images(field, image, mid, image2, top):
    if field.modulus is None:
        return None
    if image2 is None or mid.modulus is None:
        return image if mid == top else reembed(image, mid, image2, top)
    return reembed(image, mid, image2, top)


def reembed_eval(value, field, image, new, variables, point):
    return evaluate(reembed(value, field, image, new), variables, point)


# ----------------------------------------------------------------------------
# Covers over the point field.

@dataclass(frozen=True)
class NFCover:
    field: Field
    numerator: object
    half_poles: object
    lift: object
    cubic: object
    pairing: NFPairing
    orbit: int


def monic_sqrt(poly):
    """Square root of a monic polynomial of even degree, or None."""
    n = int(pari.poldegree(poly, X))
    if n % 2:
        return None
    k = n//2
    h = X**k
    for i in range(1, k+1):
        target = pari.polcoef(poly, n - i, X)
        current = pari.polcoef(h*h, n - i, X)
        h += (target - current)/2*X**(k - i)
    return h if is_zero(h*h - poly) else None


def build_cover(f, field, s, n, g):
    """Exact cover (X,Y) = (N/S^2, H y/S^3) onto Y^2 = short cubic over the point field."""
    d = s*s
    if pari.poldegree(pari.gcd(s, f), X) > 0 or pari.poldegree(pari.gcd(n, d), X) > 0:
        return None
    if max(pari.poldegree(n, X), pari.poldegree(d, X)) != 4:
        return None
    gc = [pari.polcoef(g, k, T) for k in range(4)]
    cleared = sum(gc[i]*n**i*d**(3-i) for i in range(4))
    quotient, remainder = pari.divrem(cleared, f, X)
    if not is_zero(remainder) or is_zero(quotient):
        return None
    scalar = pari.pollead(quotient, X)
    h = monic_sqrt(quotient/scalar)
    if h is None:
        return None
    shift = gc[2]/3
    short = scalar**3*pari.subst(g, T, T/scalar - shift)
    cover = (scalar*(n + shift*d), s, scalar**2*h, short)
    numerator, half, lift, cubic = cover
    sc = [pari.polcoef(cubic, k, T) for k in range(4)]
    assert is_zero(sc[3] - 1) and is_zero(sc[2])
    assert pari.poldegree(pari.gcd(cubic, pari.deriv(cubic, T)), T) == 0
    numerator, lift, cubic = normalize_scale(numerator, lift, cubic, field)
    sc = [pari.polcoef(cubic, k, T) for k in range(4)]
    identity = sum(sc[i]*numerator**i*(half*half)**(3-i) for i in range(4)) - f*lift*lift
    assert is_zero(identity)
    return numerator, half, lift, cubic


def rational_content(value):
    """Rational content of a field element in the power basis."""
    lifted = pari.lift(value)
    return abs(pari.content(lifted)) if lifted.type() == "t_POL" else abs(lifted)


def normalize_scale(numerator, lift, cubic, field):
    """Rescale Y^2 = X^3 + a4 X + a6 by a rational u, minimal at every prime of the contents.

    Over Q this is the global minimal short model away from 2 and 3; elsewhere it is the
    smallest rational rescaling, so the identity stays exact over the point field.
    """
    a4, a6 = pari.polcoef(cubic, 1, T), pari.polcoef(cubic, 0, T)
    contents = [c for c in (rational_content(a4), rational_content(a6)) if c != 0]
    primes = set()
    for c in contents:
        for part in (pari.numerator(c), pari.denominator(c)):
            primes.update(int(p) for p in pari.factor(part)[0]) if part != 1 else None
    u = pari(1)
    for p in sorted(primes):
        e = min(pari.floor(pari.valuation(c, p)/w) for c, w in zip(
            (rational_content(a4), rational_content(a6)), (4, 6)) if c != 0)
        u *= pari(p)**e
    numerator, lift, cubic = numerator/u**2, lift/u**3, pari.subst(cubic, T, u**2*T)/u**6
    if field.modulus is not None and field.degree <= 8:
        numerator, lift, cubic = ideal_scale(numerator, lift, cubic, field)
    return numerator, lift, cubic


def ideal_scale(numerator, lift, cubic, field):
    """Divide out the largest principal ideal I with I^4 | a4 and I^6 | a6 (fractionally)."""
    bnf = pari.bnfinit(field.modulus, 1)
    a4, a6 = (pari.lift(pari.polcoef(cubic, k, T)) for k in (1, 0))
    nonzero = [(c, w) for c, w in ((a4, 4), (a6, 6)) if c != 0]
    primes = set()
    for c, _ in nonzero:
        factorization = pari.idealfactor(bnf, c)
        primes.update(str(P) for P in factorization[0])
    ideal = pari(1)
    for P in (pari(p) for p in sorted(primes)):
        e = min(pari.floor(pari.idealval(bnf, c, P)/w) for c, w in nonzero)
        ideal = pari.idealmul(bnf, ideal, pari.idealpow(bnf, P, e))
    found = pari.bnfisprincipal(bnf, ideal)
    if any(int(k) for k in found[0]):
        return numerator, lift, cubic
    u = pari.Mod(pari.nfbasistoalg(bnf, found[1]), field.modulus)
    return numerator/u**2, lift/u**3, pari.subst(cubic, T, u**2*T)/u**6


def constant_on_pairs(f, numerator, half, sigma, field):
    """Exact test that X = N/S^2 takes equal values on x and sigma(x) at every root of f.

    The three branch values are distinct, so the fibre partition of the six Weierstrass
    points is then exactly the certified pairing.
    """
    f = f*pari.Mod(1, field.modulus) if field.modulus is not None else f
    moved_n, moved_s = compose(numerator, sigma, f), compose(half, sigma, f)
    return is_zero((numerator*moved_s**2 - moved_n*half**2) % f)


def differential(cover):
    """The pullback of dX/Y, as the linear polynomial w with w dx/y (exact division)."""
    numerator, half, lift = cover.numerator, cover.half_poles, cover.lift
    quotient, remainder = pari.divrem(pari.deriv(numerator, X)*half - 2*numerator*pari.deriv(half, X), lift, X)
    assert is_zero(remainder) and pari.poldegree(quotient, X) <= 1 and not is_zero(quotient)
    return quotient


def j_invariant(cubic):
    a4, a6 = pari.polcoef(cubic, 1, T), pari.polcoef(cubic, 0, T)
    return 1728*4*a4**3/(4*a4**3 + 27*a6**2)


def recover(source, log=None):
    """One representative of every Galois orbit of geometric optimal degree-four maps."""
    source, center = sextic_chart(source)
    f, orbits = pairing_orbits(source)
    output = []
    for pairing in orbits:
        algebra = paired_algebra(f, pairing)
        for degree in range(3):
            if log:
                log(f"pairing orbit {pairing.orbit} over {pairing.field.modulus}, pole degree {degree}")
            data = conditions(f, algebra, degree)
            if data is None:
                continue
            s, n, value, equations, unit = data
            variables = [(), (A,), (A, B)][degree]
            for field, image, point in solve_incidence(pairing.field, equations, unit, list(variables)):
                sub = lambda e: evaluate(reembed(e, pairing.field, image, field), variables, point)
                g = sub(pari.polresultant(algebra.modulus, T - value, U))
                cover = build_cover(sub(f), field, sub(s), sub(n), g)
                if cover is None:
                    raise ArithmeticError("A saturated point failed the exact cover identity.")
                sigma = reembed(pairing.involution, pairing.field, image, field)
                assert constant_on_pairs(f, cover[0], cover[1], sigma, field)
                output.append(NFCover(field, *cover, pairing=pairing, orbit=field.degree))
    return source, center, f, orbits, output


def text(value):
    return str(pari.lift(value)) if value is not None else None


def summary(source, center, f, orbits, covers):
    return {
        "reciprocal_chart_center": center,
        "pairing_orbits": [
            {"orbit_size": p.orbit, "field": text(p.field.modulus) or "Q", "involution": text(p.involution)}
            for p in orbits
        ],
        "geometric_maps": sum(c.orbit for c in covers),
        "maps": [
            {
                "point_field": text(c.field.modulus) or "Q",
                "conjugates": c.orbit,
                "pairing_field": text(c.pairing.field.modulus) or "Q",
                "half_poles": text(c.half_poles),
                "numerator": text(c.numerator),
                "lift": text(c.lift),
                "cubic": text(c.cubic),
                "j": text(j_invariant(c.cubic)),
                "differential": text(differential(c)),
            }
            for c in covers
        ],
    }


if __name__ == "__main__":
    import argparse
    import json

    parser = argparse.ArgumentParser(description=__doc__.splitlines()[0])
    parser.add_argument("sextic", type=json.loads,
                        help='Ascending coefficients of a rational quintic or sextic, e.g. [1,0,-2] or ["1/2",0,3].')
    args = parser.parse_args()
    poly = sp.Poly.from_list([sp.Rational(c) for c in args.sextic[::-1]], x, domain=sp.QQ)
    print(json.dumps(summary(*recover(poly)), indent=2))
