"""Check the quotient identities, torsion gluing, and class counts."""

from fractions import Fraction
from itertools import product
from math import gcd

import sympy as sp

from spectrum import (companion_degree, cyclic_kernel, describe_parameter,
                      diagonal_kernel, enumerate_spectrum, factor_formula,
                      family_degree, kernel_size, map_recipe, residue_weights)


def check_identities():
    x, t, u = sp.symbols("x t u")
    cubic = (u-1)*(u-t**2)*(u-t**4)
    curve = cubic.subs(u, x**2)
    h = x-t**2/x
    # q0 and q1, followed by the two diagonal quotients.
    for horizontal, multiplier, target in (
        (x**2, 1, cubic),
        (t**4/x**2, sp.I*t**3/x**3, cubic),
        (h, (x-sp.I*t)/x**2, u*(u-2*sp.I*t)*(u**2-(t**2-1)**2)),
        (-h, -(x+sp.I*t)/x**2, u*(u-2*sp.I*t)*(u**2-(t**2-1)**2)),
    ):
        assert sp.cancel(curve*multiplier**2-target.subs(u, horizontal)) == 0

    for horizontal, multiplier, expected in (
        (x**2, 1, 2*x),
        (t**4/x**2, sp.I*t**3/x**3, 2*sp.I*t),
        (h, (x-sp.I*t)/x**2, x+sp.I*t),
        (-h, -(x+sp.I*t)/x**2, x-sp.I*t),
    ):
        assert sp.cancel(sp.diff(horizontal, x)/multiplier-expected) == 0

    roots = (sp.Integer(1), t**2, t**4)
    w = t**4+t**2+1
    b = t**4*(t**2+2)/w
    d = t**2*(2*t**2+1)/w
    s = -t**4*(t**4-1)**2/w**2
    for root, next_root in zip(roots, roots[1:]+roots[:1]):
        assert sp.cancel(d+s/(root-b)-next_root) == 0
    curve = cubic.subs(u, x**2+b)
    eta_squared = sp.prod(d-root for root in roots)
    assert sp.cancel(cubic.subs(u, d+s/x**2)-eta_squared*curve/x**6) == 0


def check_gluing():
    kernel = diagonal_kernel()
    weights = residue_weights(2, 2, kernel)
    assert weights == {(1, 0): 1, (0, 1): 1, (1, 1): 2}
    assert residue_weights(2, 2, cyclic_kernel()) == {(1, 0): 1, (0, 1): 1, (1, 1): 1}
    for a, b in product(range(-12, 13), repeat=2):
        if gcd(a, b) == 1:
            assert family_degree(a, b) == 2*(a*a+b*b)//kernel_size((a, b), 2, kernel)
            degree, names, (r, s) = map_recipe(a, b)
            if names == ("g_plus", "g_minus"):
                assert (r+s, r-s) == (a, b) and degree == 2*(r*r+s*s)
            else:
                assert (r, s) == (a, b) and degree == 2*(a*a+b*b)
    # Residues of a proper cyclic graph in E[3]^3 give a different correction.
    graph = {((j, 0), (0, j), (0, 0)) for j in range(3)}
    assert set(residue_weights(3, 3, graph).values()) == {1}
    graph = {((j, 0), (j, 0), (0, 0)) for j in range(3)}
    weights = residue_weights(3, 3, graph)
    assert sorted(vector for vector, size in weights.items() if size == 3) == [(1, 1, 0), (2, 2, 0)]


def check_spectrum():
    counts = enumerate_spectrum(10_000)
    companion = enumerate_spectrum(10_000, companion_degree)
    for degree in range(1, 10_001):
        primes = sp.factorint(degree//2) if degree >= 2 and degree % 2 == 0 else {}
        assert counts[degree] == factor_formula(degree, primes)
        odd_part = degree//4 if degree % 4 == 0 else degree//2
        expected = 0
        if degree % 2 == 0 and odd_part % 2:
            expected = factor_formula(2*odd_part, sp.factorint(odd_part))//2
        assert companion[degree] == expected
    assert {d: counts[d] for d in (2, 4, 10, 18, 26, 34, 50, 130)} == {
        2: 4, 4: 0, 10: 8, 18: 0, 26: 8, 34: 8, 50: 8, 130: 16,
    }
    assert describe_parameter("2")[2] == Fraction(148176, 25)
    for p, q in product(range(2, 24), repeat=2):
        if p > q and gcd(p, q) == 1:
            assert describe_parameter(Fraction(p, q))[2].denominator > 1
    print(f"Verified {sum(counts.values())} and {sum(companion.values())} classes through degree 10000.")


def check_point_counts():
    checked = 0
    for prime in sp.primerange(11, 150):
        if pow(-1, (prime-1)//2, prime) != 1 or pow(35, (prime-1)//2, prime) != 1:
            continue
        squares = {x*x % prime for x in range(1, prime)}

        def affine_count(values):
            return sum(1 if value == 0 else 2 if value in squares else 0 for value in values)

        original = 2+affine_count((x*x-1)*(x*x-4)*(x*x-16) % prime for x in range(prime))
        companion = 2+affine_count((x*x+175)*(x*x+28)*(x*x-560) % prime for x in range(prime))
        elliptic = 1+affine_count((x-1)*(x-4)*(x-16) % prime for x in range(prime))
        assert original == companion == 2*elliptic-prime-1
        checked += 1
    assert checked >= 5
    print(f"Verified matching point counts at {checked} completely split primes.")


if __name__ == "__main__":
    check_identities()
    check_gluing()
    check_spectrum()
    check_point_counts()
    print("Verified six symbolic maps, the four differentials and both finite kernels.")
