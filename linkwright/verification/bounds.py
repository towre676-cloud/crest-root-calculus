"""Outward-rounded fixed-point intervals and separating-axis lower bounds (millimetres)."""

from dataclasses import dataclass
from fractions import Fraction
from math import isqrt

S = 100_000_000


def ceil_div(n, d):
    return -((-n) // d)


@dataclass(frozen=True, slots=True)
class Box:
    lo: int
    hi: int

    def __post_init__(self):
        if self.lo > self.hi:
            raise ValueError("empty interval")

    @classmethod
    def enclose(cls, lo, hi=None):
        a, b = Fraction(str(lo)) * S, Fraction(str(lo if hi is None else hi)) * S
        return cls(a.numerator // a.denominator, ceil_div(b.numerator, b.denominator))

    def __add__(self, other):
        return Box(self.lo + other.lo, self.hi + other.hi)

    def __sub__(self, other):
        return Box(self.lo - other.hi, self.hi - other.lo)

    def __neg__(self):
        return Box(-self.hi, -self.lo)

    def __mul__(self, other):
        prod = [x * y for x in (self.lo, self.hi) for y in (other.lo, other.hi)]
        return Box(min(prod) // S, ceil_div(max(prod), S))

    def __truediv__(self, other):
        if other.lo <= 0:
            raise ValueError("expected a strictly positive divisor")
        pairs = [(n * S, d) for n in (self.lo, self.hi) for d in (other.lo, other.hi)]
        return Box(min(n // d for n, d in pairs), max(ceil_div(n, d) for n, d in pairs))

    def square(self):
        near = max(self.lo, -self.hi, 0)
        far = max(abs(self.lo), abs(self.hi))
        return Box(near * near // S, ceil_div(far * far, S))

    def sqrt(self):
        if self.lo < 0:
            raise ValueError("radical enclosure crosses zero; refine inputs")
        low, high = isqrt(self.lo * S), isqrt(self.hi * S)
        return Box(low, high + int(high * high < self.hi * S))

    def rounded(self):
        return [self.lo / S, self.hi / S]


ZERO = Box(0, 0)
ONE = Box(S, S)
TWO = Box(2 * S, 2 * S)
ORIGIN = (ZERO, ZERO)


def circle_arcs(divisions):
    """Four rotated copies of tan(theta/2) in [0,1], with exact closed endpoints."""
    for quad in range(4):
        for i in range(divisions):
            t, u = Fraction(i, divisions), Fraction(i + 1, divisions)
            cos = Box.enclose((1 - u * u) / (1 + u * u), (1 - t * t) / (1 + t * t))
            sin = Box.enclose(2 * t / (1 + t * t), 2 * u / (1 + u * u))
            yield [(cos, sin), (-sin, cos), (-cos, -sin), (sin, -cos)][quad]


def plus_pose(g, a, b, c, direction):
    ax, ay = (a * d for d in direction)
    dx, dy = g - ax, -ay
    d2 = dx.square() + dy.square()
    along = (b.square() - c.square() + d2) / (TWO * d2)
    radical = ((b + c).square() - d2) * (d2 - (b - c).square())
    normal = radical.sqrt() / (TWO * d2)
    joint = (ax + along * dx - normal * dy, ay + along * dy + normal * dx)
    return {"O": ORIGIN, "D": (g, ZERO), "A": (ax, ay), "B": joint}


def direction_candidates(p, q, r, s):
    def mid(v):  # twice the nominal midpoint; only directions matter
        return v[0].lo + v[0].hi, v[1].lo + v[1].hi

    pm, qm, rm, sm = map(mid, (p, q, r, s))
    for a, b in ((pm, qm), (rm, sm)):
        yield b[1] - a[1], a[0] - b[0]
    for a in (pm, qm):
        for b in (rm, sm):
            yield a[0] - b[0], a[1] - b[1]
    yield 1, 0
    yield 0, 1


def projection(axis, pt):
    lo = hi = 0
    for n, b in zip(axis, pt):
        a, c = n * b.lo, n * b.hi
        lo += min(a, c)
        hi += max(a, c)
    return lo, hi


def segment_gap(p, q, r, s):
    """Conservative centreline gap in units of 1/S, including all endpoint boxes."""
    best = 0
    for axis in direction_candidates(p, q, r, s):
        length2 = sum(x * x for x in axis)
        if not length2:
            continue
        pl, ph = projection(axis, p)
        ql, qh = projection(axis, q)
        rl, rh = projection(axis, r)
        sl, sh = projection(axis, s)
        gap = max(min(rl, sl) - max(ph, qh), min(pl, ql) - max(rh, sh))
        if gap > 0:
            length = isqrt(length2)
            bound = gap // (length + int(length * length != length2))
            best = max(best, bound)
    return best
