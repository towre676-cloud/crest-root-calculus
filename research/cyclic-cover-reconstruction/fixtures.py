"""Reproducible source equations with independent planted-cover provenance."""

from math import gcd
import random

import sympy as sp

from cyclic import X, coefficients, power_reduce


def poly(expression):
    return sp.Poly(expression, X, domain=sp.QQ)


def integral(polynomial):
    result = polynomial.integrate()
    return result - result.nth(0)


def normalization_example(k, U, G, V):
    branch = G.compose(U).exquo(V**k)
    if gcd(branch.degree(), k) != 1 or G.gcd(G.diff()).degree() != 0:
        raise ArithmeticError("A test fixture has inadmissible branch data.")
    if G.compose(U) != branch*V**k or not U.diff().rem(V**(k-1)).is_zero:
        raise ArithmeticError("A test fixture fails its planted identity.")
    return branch


def picard():
    return poly(X**4-2*X**3+sp.Rational(7,12)*X**2+sp.Rational(5,12)*X+sp.Rational(5,18))


def saturated(k, a=0, distance=1):
    d = 2*k-1
    V = poly((X-a)*(X-a-distance))
    U = integral(d*V**(k-1))
    first, second = U.eval(a), U.eval(a+distance)
    m = next(m for m in range(2, k+1) if gcd(k, m) == 1)
    extra = [max(first, second)+i for i in range(1, m-1)]
    G = poly(sp.prod(X-r for r in [first, second, *extra]))
    return U, V, G


def nonsaturated(k, a=0, distance=1):
    V = poly((X-a+distance)*(X-a)*(X-a-distance))
    start, finish = a, a+distance
    moments = [(poly(X**i)*V**(k-1)).integrate() for i in (0, 1)]
    weights = [value.eval(finish)-value.eval(start) for value in moments]
    rho = weights[1]/weights[0]
    d = 3*(k-1)+2
    U = integral(d*V**(k-1)*poly(X-rho))
    G = poly((X-U.eval(a))*(X-U.eval(a-distance)))
    if len({U.eval(j) for j in (a-distance, a, a+distance, rho)}) != 3:
        raise ArithmeticError("The nonsaturated fixture lost its third critical value.")
    return U, V, G


def mixed(k, a=0, distance=1):
    # Over one target root U has a k-fold preimage. Over the other
    # target root it has a double preimage retained in the branch.
    V = poly(X-a)
    d = k+1
    U = integral(d*V**(k-1)*poly(X-a-distance))
    first, second = U.eval(a), U.eval(a+distance)
    m = next(m for m in range(2, k+1) if gcd(k, m) == 1)
    G = poly(sp.prod(X-r for r in [first, second, *[max(first, second)+i for i in range(1, m-1)]]))
    F = normalization_example(k, U, G, V)
    if F.gcd(F.diff()).degree() == 0:
        raise ArithmeticError("Expected a mixed branch exponent.")
    return U, V, G


def known_map(U, V, G):
    return {"inner": coefficients(U), "outer": coefficients(G), "vertical": coefficients(V)}


def build_panel(seed):
    rng = random.Random(seed)
    plans = [(3, "picard-quintic"), (3, "complement"), (5, "complement"),
             (7, "complement"), (3, "mixed"), (4, "mixed"), (5, "mixed"),
             (3, "saturated"), (4, "saturated"), (5, "saturated"),
             (7, "saturated"), (9, "saturated"), (11, "saturated")]
    fixtures = []
    for index, (k, family) in enumerate(plans):
        if family == "picard-quintic":
            U, V, G = [poly(z) for z in (X**5-5*X**4/2+5*X**3/3, X*(X-1), X*(X-sp.Rational(1,6)))]
        else:
            origin = rng.choice((-2, -1, 0, 1, 2))
            distance = rng.choice((1, 2)) if family == "mixed" else 1
            U, V, G = {"complement": nonsaturated, "mixed": mixed, "saturated": saturated}[family](k, origin, distance)
        if index % 3 == 2:
            G = G.mul_ground(sp.Rational(-2, 3))
        F = normalization_example(k, U, G, V)
        stratum = [U.degree(), G.degree(), V.degree()]
        fixtures.append((k, F, stratum, family, known_map(U, V, G), None))
        if family in {"picard-quintic", "complement", "mixed"}:
            offset = rng.choice((-2, -1, 1, 2))
            # Change one coefficient and retain its actual branch pattern.
            power = 1 if family == "picard-quintic" else 0
            changed = F + poly(offset*X**power)
            if power_reduce(changed, k)[1] != changed:
                raise ArithmeticError("A perturbation ceased to be power-free.")
            fixtures.append((k, changed, stratum, f"{family}-perturbed", None, [power, offset]))
    rng.shuffle(fixtures)
    public, provenance = {}, {}
    for i, (k, F, triple, family, planting, alteration) in enumerate(fixtures, 1):
        identifier = f"c{i:02}"
        public[identifier] = {"exponent": k, "coefficients": coefficients(F), "stratum": triple}
        provenance[identifier] = {"family": family, "planting": planting, "alteration": alteration}
    return public, provenance
