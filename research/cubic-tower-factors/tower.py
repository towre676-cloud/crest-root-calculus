"""Subset and sign factors for the Crest–Root nonabelian cubic tower."""
from itertools import combinations
from math import isqrt, prod

import numpy as np
import sympy as sp

from field import Field


def check_stage(p, indices):
    if not indices or any(i < 0 for i in indices) or len(set(indices)) != len(indices):
        raise ValueError("Supply distinct nonnegative tower indices")
    if p < 5 or not sp.isprime(p):
        raise ValueError("Use a prime greater than three")
    z = sp.Symbol("z")
    cubic = sp.Poly(-4*z**3+z**2+z+1, z, modulus=p)
    if cubic.discriminant() % p == 0:
        raise ValueError("The elliptic base has bad reduction")
    earlier = []
    for i in indices:
        branch = sp.Poly((12*i+1)*z**2 + (12*i*i+1)*z + 4*i**3+1, z, modulus=p)
        if (branch.degree() != 2 or not branch.is_sqf or branch.gcd(cubic).degree()
                or any(branch.gcd(d).degree() for d in earlier)):
            raise ValueError(f"Bad branch configuration at index {i} modulo {p}")
        earlier.append(branch)


def cubic_tables(field):
    """Histograms of x^3+c*x with c in {0,1,g}, indexed by logs (zero is -1)."""
    q, order = field.q, field.order
    tables = [np.zeros(q, dtype=np.int64) for _ in range(3)]
    for start in range(0, q, 65536):
        x = field.log[np.arange(start, min(start+65536, q))]
        cubes = np.where(x == -1, -1, 3*x % order)
        images = [cubes, field.add(cubes, x), field.add(cubes, field.shift(x, 1))]
        for table, image in zip(tables, images):
            np.add.at(table, image, 1)
    return tables


def fibre_sizes(field, tables, alpha, beta):
    """Count x^3+alpha*x=beta by scaling x to square-class representatives."""
    odd = alpha % 2
    log_r = (alpha-odd)//2
    target = field.shift(beta, -3*log_r)
    usual = np.where(odd == 0, tables[1][target], tables[2][target])
    return np.where(alpha == -1, tables[0][beta], usual)


def trace_sums(p, indices, extension, max_order):
    """Return ordinary and sigma-twisted H1 traces by nonempty support and the base."""
    check_stage(p, indices)
    q = p**extension
    if (q+1+2*(isqrt(q)+1))*3**len(indices) >= 2**63:
        raise ValueError("Fiber sums exceed the supported 64-bit accumulator bound")
    field = Field(p, extension)
    tables = cubic_tables(field)
    n = len(indices)
    supports = [s for k in range(1, min(n, max_order)+1) for s in combinations(range(n), k)]
    # Matrix rows index ordinary and sign-twisted point counts, respectively.
    sums = np.zeros((2, len(supports)), dtype=np.int64)
    base = np.zeros(2, dtype=np.int64)
    total = np.zeros(2, dtype=np.int64)
    for first in range(0, field.q, 16384):
        encoded = np.arange(first, min(first+16384, field.q))
        z = field.log[encoded]
        z2 = np.where(z == -1, -1, 2*z % field.order)
        z3 = np.where(z == -1, -1, 3*z % field.order)
        rhs = field.add(field.add(field.shift(z3, field.log[-4 % p]), z2), field.add(z, 0))
        for twist in (0, 1):
            # w=eta*w', t=eta*t' with eta^2=g^twist, eta^q=(-1)^twist eta.
            # Opposite w' give opposite cubic roots. Count one sign with weight two.
            valid = (rhs == -1) | (rhs % 2 == twist)
            base_z, values = encoded[valid], rhs[valid]
            weight = np.where(values == -1, 1, 2)
            base[twist] += weight.sum()
            root_w = np.where(values == -1, -1, ((values-twist) % field.order)//2)
            beta = field.shift(root_w, -twist)
            chars = []
            for i in indices:
                a = field.log[field.translate(base_z, i)]
                alpha = field.shift(a, field.log[3] - twist)
                chars.append(fibre_sizes(field, tables, alpha, beta) - 1)
            full = weight.copy()
            for char in chars:
                full *= char+1
            total[twist] += full.sum()
            for j, support in enumerate(supports):
                terms = weight
                for i in support:
                    terms = terms*chars[i]
                sums[twist, j] += terms.sum()

    # At the normalized infinity the single sheet is rational and the two
    # remaining sheets have b^2=-(12i+1)/12 and are interchanged by sigma.
    base += 1
    legendre = [int(sp.legendre_symbol((-(12*i+1)*pow(12, -1, p)) % p, p))**extension
                for i in indices]
    for sign in (0, 1):
        sheet_chars = [1 + (-1)**sign * e for e in legendre]
        total[sign] += prod(c+1 for c in sheet_chars)
        for j, support in enumerate(supports):
            sums[sign, j] += prod(sheet_chars[i] for i in support)
    result = {}
    for j, support in enumerate(supports):
        key = tuple(indices[i] for i in support)
        plain, twisted = (int(-v) for v in sums[:, j])
        if (plain+twisted) % 2:
            raise ArithmeticError("Nonintegral sign trace")
        result[key] = ((plain+twisted)//2, (plain-twisted)//2)
    base_traces = [field.q+1-int(count) for count in base]
    if base_traces[0] != -base_traces[1]:
        raise ArithmeticError("Elliptic inversion trace failed")
    if max_order >= n:
        plus = sum(a for a, b in result.values())
        minus = sum(b for a, b in result.values())
        if list(total) != [q+1-base_traces[0]-plus-minus, q+1+base_traces[0]-plus+minus]:
            raise ArithmeticError("Subset traces failed to reconstruct the two full point counts")
    return base_traces[0], result, [int(v) for v in total]


def weil_polynomial(traces, p, dimension):
    """Low-to-high coefficients of det(1-TF); reciprocity supplies the upper half."""
    if len(traces) < dimension:
        raise ValueError("Insufficient traces to recover the factor")
    coeffs = [1]
    for j in range(1, dimension+1):
        numerator = -sum(coeffs[j-i]*traces[i-1] for i in range(1, j+1))
        if numerator % j:
            raise ArithmeticError("Newton identities failed integrality")
        coeffs.append(numerator//j)
    return coeffs + [coeffs[dimension-i]*p**i for i in range(1, dimension+1)]


def trace_from_factor(coeffs, power):
    traces = []
    for k in range(1, power+1):
        degree = len(coeffs)-1
        linear = k*coeffs[k] if k <= degree else 0
        traces.append(-linear-sum(coeffs[j]*traces[k-j-1] for j in range(1, min(k-1, degree)+1)))
    return traces[-1]
