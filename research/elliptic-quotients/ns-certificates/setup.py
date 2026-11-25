"""Reduction of a curve and its certified maps at a good prime."""
import random
from fractions import Fraction
from pathlib import Path
import sys

import cypari2
import sympy

from ffcore import Field, CertMap
sys.path.insert(0, str(Path(__file__).resolve().parent.parent / "exact-completeness"))
from jac import Jacobian, odd_model, weil_poly, jac_order
from absence import good_model

pari = cypari2.Pari()


def _den(v):
    if isinstance(v, tuple):
        from math import lcm
        return lcm(*(Fraction(c).denominator for c in v))
    return Fraction(v).denominator


def good_primes(f, maps, need_i, start=11, stop=2000):
    for p in sympy.primerange(start, stop):
        if not good_model(f, p) or (need_i and p % 4 != 1):
            continue
        try:
            odd_model(f, p)
            field = Field(p, 1, need_i=need_i)
            check_targets(field, maps)
        except (ValueError, ZeroDivisionError):
            continue
        yield p


def check_targets(field, maps):
    for m in maps:
        values = [*m["N"], *m["D"], m["a"], m["b"], m["g2"], m["g3"]]
        if any(_den(c) % field.p == 0 for c in values):
            raise ValueError("A map coefficient has a pole at the chosen prime.")
        g2, g3 = field.q(m["g2"]), field.q(m["g3"])
        if g2**3 == 27*g3**2:
            raise ValueError("The elliptic target has bad reduction.")


class Reduction:
    def __init__(self, f, maps, p, r=1, seed=5, force_i=False, nf=None):
        need_i = force_i or any(m["gaussian"] for m in maps)
        if not good_model(f, p):
            raise ValueError("The source has bad reduction.")
        self.K = Field(p, r, need_i=need_i, nf=nf)
        check_targets(self.K, maps)
        g, alpha = odd_model(f, p)
        self.J = Jacobian(g, p, r)
        self.P = weil_poly(g, p)
        self.order = jac_order(self.P, p, r)
        self.maps = []
        for m in maps:
            candidate = CertMap(self.K, self.J, alpha, f, m["N"], m["D"], m["a"], m["b"], m["g2"], m["g3"])
            if candidate.deg != m["degree"]:
                raise ValueError("Map reduction changes its degree.")
            self.maps.append(candidate)
        self.rng = random.Random(seed)
        self.aux = {}
        self.auxJ = [self.J.random_divisor(self.rng) for _ in range(12)]

    def aux_points(self, i, n=30):
        if i not in self.aux:
            E = self.maps[i].E
            pts = []
            for P in E.points():
                pts.append(P)
                if len(pts) >= 400:
                    break
            self.rng.shuffle(pts)
            self.aux[i] = pts[:n]
        return self.aux[i]

    def pull(self, i, P):
        return self.maps[i].pullback(P, self.aux_points(i))

    def push(self, i, Dc):
        return self.maps[i].pushforward(Dc, self.auxJ)

    def ell_order(self, i):
        E = self.maps[i].E
        return sum(1 for _ in E.points()) + 1
