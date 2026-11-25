"""Certified bounded elliptic shortlist from a period enclosure."""
from itertools import product
from math import gcd, isqrt

from flint import arb_mat

from enumeration import ldl, short_vectors


PAIRS = [(0, 1), (0, 2), (0, 3), (1, 2), (1, 3), (2, 3)]


def real_frame(periods):
    return arb_mat([row for differential in periods
                    for row in ([z.real for z in differential], [z.imag for z in differential])])


def frame(data):
    P = real_frame(data["Pi"])
    J0 = arb_mat([[0, -1, 0, 0], [1, 0, 0, 0], [0, 0, 0, -1], [0, 0, 1, 0]])
    J = P.inv() * J0 * P
    G = arb_mat(data["E"]) * J
    G = (G + G.transpose()) / 2
    if G[0, 0] < 0:
        G = -G
    integral_lower_bound(G)
    return J, G


def integer(value):
    result = value.unique_fmpz()
    if result is None:
        raise ArithmeticError("The enclosure did not determine a unique integer.")
    return int(result)


def integral_lower_bound(G, scale=2 ** 32):
    """Build a positive definite integer A with A <= scale * G."""
    n = G.nrows()
    A = [[0] * n for _ in range(n)]
    radii = [0] * n
    for i in range(n):
        for j in range(i + 1):
            interval = scale * G[i, j]
            A[i][j] = A[j][i] = integer(interval.mid().floor())
            error = integer((interval - A[i][j]).abs_upper().ceil())
            radii[i] += error
            if i != j:
                radii[j] += error
    for i, radius in enumerate(radii):
        A[i][i] -= radius
    ldl(A)
    return A


def normalize(pi):
    return pi if next(x for x in pi if x) > 0 else tuple(-x for x in pi)


def dyadic_bounds(ball, scale):
    return integer((ball * scale).lower().floor()), integer((ball * scale).upper().ceil())


def elliptic_shortlist(data, J, G, X):
    """Superset of all true planes. Ball containment does not imply membership."""
    scale = 2 ** 32
    form = integral_lower_bound(G, scale)
    bound = isqrt((2 * X * scale) ** 2 // 3)
    count, planes = 0, {}
    resolution = 2 ** 45
    for u in short_vectors(form, bound):
        count += 1
        column = arb_mat([[x] for x in u])
        Ju = J * column
        length = (column.transpose() * G * column)[0, 0]
        if not length > 0:
            raise ArithmeticError("The witness vector length was not enclosed away from zero.")
        ratios = [(u[i] * Ju[j, 0] - u[j] * Ju[i, 0]) / length for i, j in PAIRS]
        limits = [dyadic_bounds(r, resolution) for r in ratios]
        # Prefer the coordinate that most often rejects a trial degree. This
        # only changes the order in which the exact inequalities are tested.
        first = max(range(6), key=lambda i: abs(limits[i][0] + limits[i][1]))
        low, high = limits[first]
        for d in range(1, X + 1):
            if -((-d * low) // resolution) > d * high // resolution:
                continue
            ranges = [range(-((-d * a) // resolution), (d * b) // resolution + 1)
                      for a, b in limits]
            if any(not values for values in ranges):
                continue
            for pi in product(*ranges):
                if gcd(*pi) != 1 or pi[0] * pi[5] - pi[1] * pi[4] + pi[2] * pi[3]:
                    continue
                e = abs(sum(data["E"][i][j] * c for (i, j), c in zip(PAIRS, pi)))
                if e == d and all((d * r - c).contains(0) for r, c in zip(ratios, pi)):
                    planes[normalize(pi)] = d
    return planes, count


def endomorphism_form(G, J, penalty):
    """Ambient positive form agreeing with the Rosati norm on complex-linear matrices."""
    inverse = G.inv()
    indices = [(i, j) for i in range(4) for j in range(4)]
    Q = arb_mat([[G[a, c] * inverse[d, b] for c, d in indices] for a, b in indices])
    C = arb_mat([[(J[a, c] if b == d else 0) - (J[d, b] if a == c else 0)
                  for c, d in indices] for a, b in indices])
    form = Q + penalty * C.transpose() * C
    return (form + form.transpose()) / 2
