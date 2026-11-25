"""Exact degree counts for two quadratic integer rings and curve-model checks."""
import argparse
from collections import Counter
from itertools import product
import json
from math import isqrt
from pathlib import Path
import sys

import sympy as sp


HERE = Path(__file__).resolve().parent
NS = HERE.parent / "ns-certificates"
sys.path.insert(0, str(NS))
from isotropic import spectrum_from_lattice


def gaussian_formula(degree):
    """Closed geometric multiplicity for a non-product Gaussian two-gluing."""
    odd = degree
    while odd % 2 == 0:
        odd //= 2
    two = degree // odd
    if two >= 8 or (two == 1 and odd % 4 == 1):
        return 0
    value = 2
    for p, e in sp.factorint(odd).items():
        value *= p ** (e - 1) * (p - (1 if p % 4 == 1 else -1))
    return 2 * value if two == 2 and odd % 4 == 1 else value


def multiply(x, y, t, u):
    a, b = x
    c, d = y
    return a * c - u * b * d, a * d + b * c + t * b * d


def norm(v, t, u):
    a, b = v
    return a*a + t*a*b + u*b*b


def primitive(x, y, t, u):
    """Euclidean algorithm in the Gaussian or discriminant -7 ring."""
    while y != (0, 0):
        c, d = y
        den = norm(y, t, u)
        a, b = multiply(x, (c + t*d, -d), t, u)
        q1 = (2*b + den) // (2*den)
        q0 = (2*a + t*b - t*q1*den + den) // (2*den)
        l, r = multiply((q0, q1), y, t, u)
        x, y = y, (x[0] - l, x[1] - r)
    return norm(x, t, u) == 1


GAUSSIAN = (0, 3, 2, 1)
EIGENLINE_SWAP = (0, 2, 1, 3)
SUM_SWAP = (0, 1, 3, 2)


def gluing_table(t, u, gamma):
    residues = tuple(product(range(2), repeat=2))
    if t == 0:
        # Basis 1, i in O/2. The maximal ideal is the kernel of a+b mod 2.
        label = lambda v: (v[0] % 2) + 2*(v[1] % 2)
        admissible = lambda x, y: any((a+b) % 2 for a, b in (x, y))
    else:
        # Evaluate the basis element at the two roots 0 and 1 modulo 2.
        label = lambda v: (v[0] % 2) + 2*((v[0]+v[1]) % 2)
        admissible = lambda x, y: label(x) | label(y) == 3
    return {
        (x, y): sum(label(multiply(y, z, t, u)) == gamma[label(multiply(x, z, t, u))]
                    for z in residues)
        for x, y in product(residues, repeat=2) if admissible(x, y)
    }


def enumeration(bound, t, u, table, units):
    # A coordinate has norm <= floor(bound*max(k)/2) and
    # 2*norm(a,b) >= a*a+b*b in either of the two rings.
    cap = bound * max(table.values()) // 2
    lim = isqrt(2*cap)
    points = []
    for v in product(range(-lim, lim+1), repeat=2):
        size = norm(v, t, u)
        if size <= cap:
            points.append((size, v, tuple(c % 2 for c in v)))
    points.sort()
    counts = Counter()
    for a, x, ox in points:
        for b, y, oy in points:
            if a+b > cap:
                break
            k = table.get((ox, oy))
            if not k or 2*(a+b) > k*bound or not primitive(x, y, t, u):
                continue
            degree, remainder = divmod(2*(a+b), k)
            if remainder:
                raise ArithmeticError("Nonintegral polarization degree.")
            counts[degree] += 1
    if any(count % units for count in counts.values()):
        raise ArithmeticError("The action of units is not free on the counted vectors.")
    return {degree: count//units for degree, count in sorted(counts.items())}


def verify_models():
    """Exact models for the discriminant -7 pair, including both degree-two maps."""
    u, x = sp.symbols("u x")
    delta = sp.sqrt(-7)
    source = ((7+delta)/2, (7-delta)/2, -7)
    h = u**3 - 35*u + 98
    if sp.expand(sp.prod(u-a for a in source) - h) != 0:
        raise ArithmeticError("Incorrect elliptic branch points.")
    if sp.Rational(1728) * 4*(-35)**3 / (4*(-35)**3 + 27*98**2) != -3375:
        raise ArithmeticError("Incorrect j-invariant.")
    records = []
    for perm in (EIGENLINE_SWAP, SUM_SWAP):
        target = [source[perm[i+1]-1] for i in range(3)]
        q = sp.simplify((target[1]-target[0])/(target[2]-target[0]))
        a, b, c = source
        s = sp.simplify((q*(a-c)*b - (a-b)*c)/((a-b)-q*(a-c)))
        r = [sp.simplify(v+s) for v in source]
        if 0 in r or len(set(r)) != 3:
            raise ArithmeticError("Singular sextic.")
        lam = sp.simplify((target[1]-target[0])/(1/r[1]-1/r[0]))
        mu = sp.simplify(target[0]-lam/r[0])
        if any(sp.simplify(mu+lam/a-b) != 0 for a, b in zip(r, target)):
            raise ArithmeticError("The two-torsion labels did not match.")
        coefficient = sp.simplify(h.subs(u, mu))
        f = sp.expand(sp.prod(x*x-v for v in r))
        check = sp.together(x**6*h.subs(u, mu+lam/x**2)).expand()
        if coefficient == 0 or sp.simplify(check-coefficient*f) != 0 or sp.simplify(h.subs(u, x*x-s)-f) != 0:
            raise ArithmeticError("The elliptic maps failed exact substitution.")
        records.append({"shift": str(s), "lambda": str(lam), "mu": str(mu),
                        "eta_squared": str(coefficient), "sextic": str(f)})
    return records


def verify_nonprincipal_transfer():
    """Compare the torsion intersections from O and (2, 1+sqrt(-5))."""
    table = gluing_table(0, 5, (0, 2, 1, 3))
    reps = tuple(product(range(2), repeat=2))

    def label(v):
        return v[0] % 2 + 2*(v[1] % 2)

    def ideal_action(u, torsion):
        # lambda*u is represented in the ideal basis (2, lambda).
        a, b = u
        c, d = torsion
        xi, zeta = -3*b, a+b
        v = (2*xi + zeta, zeta)
        # Multiply by 2P for P=(c+d*lambda/2)/2, then read E[2].
        x, y = multiply(v, (2*c+d, d), 0, 5)
        if x % 2 or y % 2:
            raise ArithmeticError("The map did not land in E[2].")
        return x//2, y//2

    # lambda=1+sqrt(-5) gives an isogeny E_a -> E of norm-degree 3.
    torsion_images = {(c - 2*d, c + d) for c, d in reps}
    if len({label(p) for p in torsion_images}) != 4:
        raise ArithmeticError("The torsion transfer is not invertible modulo two.")
    gamma = (0, 2, 1, 3)
    for (a, b), expected in table.items():
        intersection = sum(label(ideal_action(b, p)) == gamma[label(ideal_action(a, p))]
                           for p in reps)
        if intersection != expected:
            raise ArithmeticError("The finite correction changed with the ideal class.")
    return len(table)


def verify(bound=150):
    table = gluing_table(0, 1, GAUSSIAN)
    if Counter(table.values()) != {1: 6, 2: 6}:
        raise ArithmeticError("Unexpected Gaussian kernel distribution.")
    formula = {d: m for d in range(1, bound+1) if (m := gaussian_formula(d))}
    if formula != enumeration(bound, 0, 1, table, 4):
        raise ArithmeticError("Gaussian divisor formula disagrees with exact primitive vectors.")
    for suffix in ("102", "120"):
        stored = json.loads((NS / f"cert_cm_{suffix}.json").read_text(encoding="utf-8"))
        lattice, _ = spectrum_from_lattice(stored["NS_gram"], stored["theta"], bound)
        if lattice != formula:
            raise ArithmeticError(f"The Gaussian formula disagrees with NS certificate {suffix}.")
    tables = [gluing_table(1, 2, gamma) for gamma in (EIGENLINE_SWAP, SUM_SWAP)]
    if [Counter(t.values()) for t in tables] != [{1: 6, 2: 3}, {1: 4, 2: 5}]:
        raise ArithmeticError("Unexpected discriminant -7 kernel distribution.")
    split_counts = [enumeration(bound, 1, 2, t, 2) for t in tables]
    if [d.get(3, 0) for d in split_counts] != [0, 8] or any(min(d) != 2 for d in split_counts):
        raise ArithmeticError("The split spectra did not distinguish two genus-two examples.")
    report = {"bound": bound, "nonprincipal_class_residues": verify_nonprincipal_transfer(),
              "gaussian": {"classes": sum(formula.values()), "distinct_degrees": len(formula),
                           "initial_counts": {d: m for d, m in formula.items() if d <= 12}},
              "discriminant_minus_7": [
                  {"model": model, "classes": sum(counts.values()),
                   "initial_counts": {d: m for d, m in counts.items() if d <= 12}}
                  for model, counts in zip(verify_models(), split_counts)]}
    print(json.dumps(report, indent=2))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--bound", type=int, default=150)
    args = parser.parse_args()
    if args.bound < 12:
        parser.error("Use a bound of at least 12 to include the distinguishing controls.")
    verify(args.bound)
