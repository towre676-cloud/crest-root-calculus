"""Exact six-bar configuration examples. Python 3.10+, standard library only."""

import argparse
import json
from collections import Counter
from fractions import Fraction as Q
from itertools import combinations, product
from pathlib import Path


def determinant(rows):
    a, b, c = rows
    return (
        a[0] * (b[1] * c[2] - b[2] * c[1])
        - a[1] * (b[0] * c[2] - b[2] * c[0])
        + a[2] * (b[0] * c[1] - b[1] * c[0])
    )


def discriminant(form):
    a, b, c = form
    return b * b - 4 * a * c


def resultant(f, g):
    a, b, c = f
    d, e, h = g
    return (a * h - c * d) ** 2 - (a * e - b * d) * (b * h - c * e)


ALIGNED = ((-3, 0, 5), (7, 0, -1), (12, 0, -6))
ROTATED = ((0, 1, 0), (22, -49, 22), (120, -161, -120))


def invariants(forms):
    return {
        "determinant": determinant(forms),
        "discriminants": [discriminant(f) for f in forms],
        "resultants": [resultant(f, g) for f, g in combinations(forms, 2)],
    }


def form(f, u, v):
    a, b, c = f
    return a * u * u + b * u * v + c * v * v


def heron(crank, ground, coupler, rocker, u, v):
    # Squared-distance numerator, with denominator S=u²+v².
    gx, gy = ground
    s = u * u + v * v
    dist = (gx * gx + gy * gy + crank * crank) * s - 2 * crank * (
        gx * (v * v - u * u) + 2 * gy * u * v
    )
    return ((coupler + rocker) ** 2 * s - dist) * (dist - (coupler - rocker) ** 2 * s)


def check_mechanisms():
    # Five distinct affine values certify equality of each quartic polynomial.
    for u, v in [(u, 1) for u in range(5)] + [(1, 0)]:
        a, b, c = [form(f, u, v) for f in ALIGNED]
        assert heron(4, (6, 0), 6, 2, u, v) == 12**2 * a * b
        assert heron(4, (6, 0), 7, 1, u, v) == 8**2 * a * c
        a, b, c = [form(f, u, v) for f in ROTATED]
        assert heron(28, (0, 21), 24, 11, u, v) == 336**2 * a * b
        assert heron(28, (8, 15), 32, 7, u, v) == 112**2 * a * c
        assert (
            120**2 * (b + 49 * a) ** 2 - 22**2 * (c + 161 * a) ** 2
            == 4 * (22 * 120) ** 2 * a * a
        )


class Field:
    """Small extension fields in an integer-coded power basis, degrees 1 to 3."""

    def __init__(self, p, n):
        self.p, self.n, self.q = p, n, p**n
        self.modulus = [0] if n == 1 else self.pick_modulus()
        self.digits = [tuple(i // p**k % p for k in range(n)) for i in range(self.q)]
        self.weights = [p**i for i in range(n)]
        self.squares = Counter(self.mul(x, x) for x in range(self.q))

    def pick_modulus(self):
        # In degree 2 or 3 absence of linear factors is an irreducibility test.
        for b, c in product(range(self.p), range(1, self.p)):
            if all((x**self.n + b * x + c) % self.p for x in range(self.p)):
                return [c, b] + [0] * (self.n - 2)
        raise RuntimeError("No irreducible modulus found.")

    def linear(self, terms):
        return sum(
            (sum(s * self.digits[x][k] for s, x in terms) % self.p) * w
            for k, w in enumerate(self.weights)
        )

    def mul(self, x, y):
        if self.n == 1:
            return x * y % self.p
        a, b = self.digits[x], self.digits[y]
        v = [0] * (2 * self.n - 1)
        for i, c in enumerate(a):
            for j, d in enumerate(b):
                v[i + j] += c * d
        for h in range(len(v) - 1, self.n - 1, -1):
            for j, c in enumerate(self.modulus):
                v[h - self.n + j] -= v[h] * c
        return sum(v[i] % self.p * w for i, w in enumerate(self.weights))


def count_covers(forms, f):
    q = f.q
    covers = [0, 0, 0]
    total = 0
    shared = []
    for t in range(q + 1):
        vals = (
            [a % f.p for a, _, _ in forms]
            if t == q
            else [f.linear([(a, f.mul(t, t)), (b, t), (c, 1)]) for a, b, c in forms]
        )
        a, b, c = vals
        sizes = [f.squares[f.mul(s, r)] for s, r in ((a, b), (a, c), (b, c))]
        covers = [n + k for n, k in zip(covers, sizes)]
        if a == 0:
            # At the node, divide the product of the two radicals by A.
            # Its finite limiting values are the square roots of BC.
            total += sizes[2]
            shared.append({"base": "infinity" if t == q else t, "points": sizes[2]})
        else:
            total += sizes[0] * sizes[1]
    assert total == sum(covers) - 2 * (q + 1)
    return total, covers, shared


def count_plane(f):
    # Independently enumerate the projective quartic, using a histogram for Y.
    def sq(x):
        return f.mul(x, x)

    def left(y, z):
        val = f.linear([(1, sq(y)), (49, sq(z))])
        return f.linear([(120**2, sq(val))])

    def right(x, z):
        val = f.linear([(1, sq(x)), (161, sq(z))])
        return f.linear([(22**2, sq(val)), (4 * (22 * 120) ** 2, sq(sq(z)))])

    counts = Counter(left(y, 1) for y in range(f.q))
    return (
        sum(counts[right(x, 1)] for x in range(f.q))
        + sum(left(y, 0) == right(1, 0) for y in range(f.q))
        + int(left(1, 0) == right(0, 0))
    )


def poly_product(a, b):
    out = [0] * (len(a) + len(b) - 1)
    for i, c in enumerate(a):
        for j, d in enumerate(b):
            out[i + j] += c * d
    return out


def local_check(forms, p, extension_degrees, plane=False):
    info = invariants(forms)
    assert all(x % p for x in info["discriminants"] + info["resultants"]) and p != 2
    output = []
    for n in extension_degrees:
        f = Field(p, n)
        count, quotients, common = count_covers(forms, f)
        row = {
            "degree": n,
            "q": f.q,
            "C": count,
            "quotient_counts": quotients,
            "common_branch_fibres": common,
        }
        if plane:
            row["direct_plane"] = count_plane(f)
            assert count == row["direct_plane"]
        output.append(row)
    traces = [p + 1 - k for k in output[0]["quotient_counts"]]
    lp = [1]
    for a in traces:
        lp = poly_product(lp, [1, -a, p])
    if extension_degrees == [1, 2, 3]:
        s1, s2, s3 = [r["q"] + 1 - r["C"] for r in output]
        c1 = -s1
        c2 = (s1 * s1 - s2) // 2
        c3 = -(s1**3 - 3 * s1 * s2 + 2 * s3) // 6
        assert lp == [1, c1, c2, c3, p * c2, p * p * c1, p**3]
    return {
        "prime": p,
        "traces": traces,
        "factor_product_coefficients": lp,
        "counts": output,
    }


D = Q(2017, 2401)


def cyclic_modular(p):
    """Evaluate the written quotient map at every projective point over Fp."""

    def red(x):
        return x.numerator * pow(x.denominator, -1, p) % p

    d = red(D)
    delta = d * d % p
    k = 6 * delta % p
    l = (1 + 2 * delta) % p
    assert all(v % p for v in (2, 3, 7, d, 1 - d * d, 4 * d * d - 1, k - l, k - 2 * l))

    def quartic(x, y, z):
        return (k * (x**4 + y**4 + z**4) - l * (x * x + y * y + z * z) ** 2) % p

    points = (
        [(x, y, 1) for x, y in product(range(p), repeat=2)]
        + [(x, 1, 0) for x in range(p)]
        + [(1, 0, 0)]
    )
    on = [P for P in points if not quartic(*P)]
    images = Counter()
    inf = 0
    for x, y, z in on:
        s1 = (x + y + z) % p
        if not s1:
            inf += 1
            continue
        u = (x * y + x * z + y * z) * pow(s1, -2, p) % p
        w = x * y * z * pow(s1, -3, p) % p
        v = (x - y) * (y - z) * (z - x) * pow(s1, -3, p) % p
        assert (
            24 * delta * w - l * (1 - 2 * u) ** 2 + k * (1 - 4 * u + 2 * u * u)
        ) % p == 0
        assert (v * v - (u * u - 4 * u**3 - 4 * w - 27 * w * w + 18 * u * w)) % p == 0
        images[u, v] += 1
    assert images and set(images.values()) == {3}
    return {
        "prime": p,
        "plane_points": len(on),
        "points_on_s1_zero": inf,
        "rational_image_points_on_chart": len(images),
        "rational_preimages_per_image": 3,
    }


def cyclic_geometry():
    # Write an element of Q(sqrt(3)) as (a,b), a+b sqrt(3).
    def times(x, y):
        return (x[0] * y[0] + 3 * x[1] * y[1], x[0] * y[1] + x[1] * y[0])

    assert Q(47, 49) ** 2 + Q(24, 49) ** 2 / 3 == 1
    assert D == 2 * Q(47, 49) ** 2 - 1
    assert 147**2 + 245**2 + 147 * 245 == 343**2
    for u, v in [(u, 1) for u in range(5)] + [(1, 0)]:
        s = u * u + v * v
        cx = v * v - u * u
        sy = 2 * u * v
        h0 = (cx - D * s, 0)
        hp = (-Q(1, 2) * cx - D * s, Q(sy, 2))
        hm = (hp[0], -hp[1])
        for sign, h in [(1, hp), (-1, hm)]:
            # Squared crank-to-ground distance numerator of the unit-scale loop.
            dist = (2 * 49**2 * s - 49**2 * cx, -sign * 49**2 * sy)
            f = times(
                ((47 + 24) ** 2 * s - dist[0], -dist[1]),
                (dist[0] - (47 - 24) ** 2 * s, dist[1]),
            )
            comb = times(h0, h)
            assert f == tuple(-((2 * 49**2) ** 2) * c for c in comb)
        squares = [times(h, h) for h in (h0, hp, hm)]
        total = (sum(z[0] for z in (h0, hp, hm)), sum(z[1] for z in (h0, hp, hm)))
        assert tuple(6 * D * D * sum(q[i] for q in squares) for i in (0, 1)) == tuple(
            (1 + 2 * D * D) * z for z in times(total, total)
        )
    # The cyclic quotient chart has an origin candidate over Q(sqrt(-3)).
    u = Q(1, 2)
    k = 6 * D**2
    l = 1 + 2 * D**2
    w = (l * (1 - 2 * u) ** 2 - k * (1 - 4 * u + 2 * u * u)) / (24 * D**2)
    discriminant = u * u - 4 * u**3 - 4 * w - 27 * w * w + 18 * u * w
    assert w == Q(1, 8) and discriminant == Q(-3, 64)
    # At theta=0, the distance in each unscaled triangle is 49.
    assert 47 - 24 < 49 < 47 + 24


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument(
        "--write", action="store_true", help="Write receipt.json beside this file."
    )
    args = parser.parse_args()

    check_mechanisms()
    cyclic_geometry()
    aligned = local_check(ALIGNED, 13, [1, 2, 3])
    small = [local_check(ROTATED, p, [1, 2, 3], plane=True) for p in (7, 13)]
    obstructions = [local_check(ROTATED, p, [1], plane=True) for p in (23, 71)]
    for i, j in combinations(range(3), 2):
        assert any(
            (row["traces"][i] == 0) != (row["traces"][j] == 0) for row in obstructions
        )

    result = {
        "aligned_hyperelliptic": {
            "forms": ALIGNED,
            **invariants(ALIGNED),
            "local": aligned,
        },
        "rotated_plane_quartic": {
            "forms": ROTATED,
            **invariants(ROTATED),
            "local": small,
            "ordinary_supersingular_witnesses": obstructions,
        },
        "cyclic_plane_quartic": {
            "d": str(D),
            "first_loop": [147, 147, 141, 72],
            "second_loop": [245, 245, 235, 120],
            "rocker_pivot_distance": 343,
            "quotient_point": {"U": "1/2", "W": "1/8", "V_squared": "-3/64"},
            "quotient_chart_replays": [cyclic_modular(p) for p in (31, 37, 43, 61)],
        },
    }
    data = json.dumps(result, indent=2) + "\n"
    if args.write:
        Path(__file__).with_name("receipt.json").write_text(data)
    print(data, end="")


if __name__ == "__main__":
    main()
