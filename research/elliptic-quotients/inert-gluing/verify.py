"""Exact local controls for the odd inert CM self-gluing theorem."""
import argparse
from collections import Counter, defaultdict
from dataclasses import dataclass
from fractions import Fraction
from itertools import combinations, product
from math import gcd, isqrt, lcm
import json


ORDERS = {-3: (1, 1), -4: (0, 1), -7: (1, 2)}
CONTROLS = {-3: (5, 11, 17), -4: (3, 7, 11, 19), -7: (3, 5, 13)}


@dataclass(frozen=True)
class ResidueField:
    p: int
    trace_omega: int
    norm_omega: int

    def reduce(self, a):
        return tuple(x % self.p for x in a)

    def add(self, a, b):
        return self.reduce((a[0] + b[0], a[1] + b[1]))

    def neg(self, a):
        return self.reduce((-a[0], -a[1]))

    def conj(self, a):
        return self.reduce((a[0] + self.trace_omega * a[1], -a[1]))

    def mul(self, a, b):
        x, y = a
        z, w = b
        return self.reduce((x*z - self.norm_omega*y*w, x*w + y*z + self.trace_omega*y*w))

    def norm(self, a):
        x, y = a
        return (x*x + self.trace_omega*x*y + self.norm_omega*y*y) % self.p

    def trace(self, a):
        return (2*a[0] + self.trace_omega*a[1]) % self.p

    def inv(self, a):
        inverse_norm = pow(self.norm(a), -1, self.p)
        return self.mul((inverse_norm, 0), self.conj(a))

    def gamma(self, a, b, t):
        return self.add(self.mul(a, t), self.mul(b, self.conj(t)))

    def legendre(self, x):
        value = pow(x % self.p, (self.p-1)//2, self.p)
        return -1 if value == self.p-1 else value


def determinant(a, b):
    return a[0]*b[1] - a[1]*b[0]


def form_gram(field, a):
    """In the basis (p,0),(-Tr(a),1),(-Tr(a*omega),omega)."""
    t, n = field.trace_omega, field.norm_omega
    row = (field.p, -field.trace(a), -field.trace(field.mul(a, (0, 1))))
    base = ((0, 0, 0), (0, 4, 2*t), (0, 2*t, 4*n))
    return [[x*y + b for y, b in zip(row, line)] for x, line in zip(row, base)]


def graph_pairing(field, a, b):
    """Four coefficients from an alternating pairing on the two graph generators."""
    one = (1, 0)
    omega = (0, 1)
    g = field.gamma(a, b, one)
    h = field.gamma(a, b, omega)
    # Pair (one,g) with H(omega,h), using the standard symplectic form on each factor.
    crossed = []
    for c in (one, omega):
        crossed.append(determinant(one, field.mul(c, h))
                       + determinant(g, field.mul(field.conj(c), omega)))
    return tuple(x % field.p for x in (1, determinant(g, h), *crossed))


def monomial_orbit(field, units, a, b):
    # Expand both factor orders and all unit changes on each axis.
    for first, second in ((a, b), (field.neg(field.conj(a)), b)):
        for u, v in product(units, repeat=2):
            yield (field.mul(field.mul(v, field.conj(u)), first),
                   field.mul(field.mul(v, u), second))


def humbert_orbit(field, units, a):
    return {field.mul(unit, c) for unit, c in product(units, (a, field.conj(a)))}


def square_shells(field, bound):
    t, n = field.trace_omega, field.norm_omega
    pairs = defaultdict(list)
    # In these orders 2*N(x+y*omega) >= x*x+y*y.
    for x, y in product(range(-bound, bound+1), repeat=2):
        norm = x*x + t*x*y + n*y*y
        if 4*norm <= bound*bound:
            pairs[norm].append((x, y))
    return [[(u, x, y) for u in range(-d, d+1, 2)
             for x, y in pairs[(d*d-u*u)//4]]
            for d in range(1, bound+1)]


def degree_counts(field, a, shells):
    s = field.trace(a)
    t = field.trace(field.mul(a, (0, 1)))
    result = []
    for shell in shells:
        count = 0
        for u, x, y in shell:
            z, remainder = divmod(u+s*x+t*y, field.p)
            if remainder == 0 and gcd(z, x, y) == 1:
                count += 1
        result.append(count)
    return result


def exact_theta_level(gram):
    """Smallest N for which N*(2*gram)^-1 is even integral."""
    m = [[2*x for x in row] for row in gram]
    def cofactor(i, j):
        # Cyclic row and column orders incorporate the cofactor sign.
        r, s = (i+1) % 3, (i+2) % 3
        t, u = (j+1) % 3, (j+2) % 3
        return m[r][t]*m[s][u] - m[r][u]*m[s][t]

    cofactors = [[cofactor(i, j) for j in range(3)] for i in range(3)]
    determinant = sum(v * c for v, c in zip(m[0], cofactors[0]))
    return lcm(*(Fraction(cofactors[i][j], determinant*(2 if i == j else 1)).denominator
                 for i, j in product(range(3), repeat=2)))


def sturm_index(level):
    factor = Fraction(level, 6)
    remainder, p = level, 2
    while p*p <= remainder:
        if remainder % p == 0:
            factor *= Fraction(p+1, p)
            while remainder % p == 0:
                remainder //= p
        p += 1
    if remainder > 1:
        factor *= Fraction(remainder+1, remainder)
    return int(factor)


def check_case(discriminant, prime, report_forms=False):
    t, n = ORDERS[discriminant]
    field = ResidueField(prime, t, n)
    if (prime < 3 or prime % 2 == 0 or any(prime % d == 0 for d in range(2, isqrt(prime)+1))
            or field.legendre(discriminant) != -1):
        raise ValueError("Choose an odd inert prime for the selected order.")

    units = [field.reduce((x, y)) for x, y in product(range(-1, 2), repeat=2)
             if x*x + t*x*y + n*y*y == 1]
    elements = tuple(product(range(prime), repeat=2))
    by_norm = defaultdict(list)
    for a in elements:
        by_norm[field.norm(a)].append(a)
    lines = [(1, x) for x in range(prime)] + [(0, 1)]
    inverse = {x: field.inv(x) for x in lines}
    form_fibers = defaultdict(set)
    genus_counts = Counter()
    form_records = []
    graphs = 0
    form_genus = {}
    level = -4 * discriminant * prime**2
    for a in elements:
        radius = (1 + field.norm(a)) % prime
        if radius == 0:
            continue
        center_key = min(humbert_orbit(field, units, a))
        sign = field.legendre(radius)
        genus_counts[sign] += 1
        circle = {r for r in elements if field.norm(field.add(r, field.neg(a))) == radius}
        if len(circle) != prime + 1 or (0, 0) in circle:
            raise ArithmeticError("Incorrect set of discounted directions.")
        divisible = sum(field.norm(r) == prime-1 for r in circle)
        if divisible != 1-sign:
            raise ArithmeticError("The degree divisibility statistic disagrees with the genus.")

        gram = form_gram(field, a)
        # On the kernel hyperplane, these two basis vectors give a unimodular block.
        det2 = gram[1][1]*gram[2][2] - gram[1][2]**2
        if det2 % prime == 0 or field.legendre(det2) != field.legendre(-4*discriminant)*sign:
            raise ArithmeticError("Genus squareclass disagrees with the hyperplane determinant.")
        if center_key == a:
            form_genus[a] = sign
            if exact_theta_level(gram) != level:
                raise ArithmeticError("The inverse Gram matrix gives a different theta level.")
            if report_forms:
                form_records.append({"scalar_part": a, "gram": gram, "genus_sign": sign})
        for b in by_norm[radius]:
            graphs += 1
            actual = {field.mul(field.gamma(a, b, x), inverse[x]) for x in lines}
            if actual != circle:
                raise ArithmeticError("The graph intersections disagree with the norm circle.")
            expected = (1, prime-1, field.trace(a), field.trace(field.mul(a, (0, 1))))
            if graph_pairing(field, a, b) != expected:
                raise ArithmeticError("The exact alternating restriction disagrees with the trace condition.")
            representative = min(monomial_orbit(field, units, a, b))
            if min(humbert_orbit(field, units, representative[0])) != center_key:
                raise ArithmeticError("A polarized graph orbit crossed two Humbert classes.")
            form_fibers[center_key].add(representative)

    scalar_graphs = by_norm[prime-1]
    if len(scalar_graphs) != prime+1 or graphs+len(scalar_graphs) != prime*(prime*prime-1):
        raise ArithmeticError("The scalar and nonscalar graphs did not exhaust the determinant coset.")
    scalar_weight = prime**4 + prime**2
    nonscalar_weight = prime**2*(prime+1) + prime**2-prime
    for x, y in scalar_graphs:
        norm = x*x + t*x*y + n*y*y
        first_diagonal, remainder = divmod(1+norm, prime)
        if remainder or first_diagonal*prime-norm != 1:
            raise ArithmeticError("The scalar polarization is not integral and unimodular.")

    if graphs != (prime+1)*(prime*prime-prime-1):
        raise ArithmeticError("The nonscalar graph count disagrees with the determinant coset.")
    if genus_counts != {1: (prime*prime-2*prime-1)//2, -1: (prime*prime-1)//2}:
        raise ArithmeticError("The finite genus partition has the wrong sizes.")
    if len(units) == 2:
        if len(form_fibers) != (prime*prime+prime-1-field.legendre(-1))//4:
            raise ArithmeticError("The Humbert class count disagrees with Burnside's formula.")
        for a, fiber in form_fibers.items():
            count = (prime+1)//2 if a == field.conj(a) else prime+1
            if len(fiber) != count:
                raise ArithmeticError("Unexpected number of polarized surfaces in a Humbert fiber.")

    out = {"discriminant": discriminant, "prime": prime, "nonscalar_graphs": graphs,
           "scalar_graphs": len(scalar_graphs),
           "scalar_to_nonscalar_leading_constant": str(Fraction(scalar_weight, nonscalar_weight)),
           "theta_level": level, "weight_two_sturm_bound_at_theta_level": sturm_index(level),
           "prime_coprime_degree_density_by_genus": {
               sign: str(Fraction(prime*(prime+sign), prime**2+2*prime-1)) for sign in (-1, 1)},
           "Humbert_isometry_classes": len(form_fibers),
           "Humbert_classes_by_genus_sign": dict(sorted(Counter(form_genus.values()).items())),
           "polarized_isomorphism_classes": sum(map(len, form_fibers.values())),
           "surface_fiber_sizes": dict(sorted(Counter(map(len, form_fibers.values())).items()))}
    cutoff = 2*prime
    shells = square_shells(field, cutoff)
    spectra = {a: degree_counts(field, a, shells) for a in form_genus}
    unresolved, last_first_difference = [], 0
    for a, b in combinations(form_genus, 2):
        if form_genus[a] != form_genus[b]:
            continue
        index = next((d for d, (x, y) in enumerate(zip(spectra[a], spectra[b]), 1) if x != y), None)
        if index is None:
            unresolved.append([a, b])
        else:
            last_first_difference = max(last_first_difference, index)
    out["same_genus_spectral_comparison"] = {
        "checked_through_degree": cutoff, "largest_first_difference": last_first_difference,
        "unresolved_pairs": unresolved,
    }
    if report_forms:
        for record in form_records:
            record["primitive_square_counts"] = spectra[record["scalar_part"]]
        out["forms"] = form_records
    if (discriminant, prime) == (-7, 5):
        one, two = (0, 0), (2, 3)
        if (field.legendre(1 + field.norm(one)) != field.legendre(1 + field.norm(two))
                or [degree_counts(field, a, shells)[1] for a in (one, two)] != [2, 0]):
            raise ArithmeticError("The same-genus degree-two control failed.")
        out["same_genus_control"] = {"scalar_parts": [one, two], "degree_two_counts": [2, 0]}
    if (discriminant, prime) == (-4, 19):
        a, b = (0, 3), (0, 7)
        left = spectra[min(humbert_orbit(field, units, a))]
        right = spectra[min(humbert_orbit(field, units, b))]
        if (form_genus[a] != form_genus[b] or left[:25] != right[:25]
                or [left[25], right[25]] != [4, 0]):
            raise ArithmeticError("The level-19 delayed separation failed.")
        out["delayed_separation"] = {"scalar_parts": [a, b], "first_difference": 26,
                                     "degree_26_counts": [4, 0]}
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--order", type=int, choices=tuple(ORDERS))
    parser.add_argument("--prime", type=int)
    parser.add_argument("--forms", action="store_true", help="Include one Gram matrix per isometry class.")
    args = parser.parse_args()
    if args.prime is not None and args.order is None:
        parser.error("--prime requires --order.")
    controls = [(d, p) for d, primes in CONTROLS.items() for p in primes]
    if args.order is not None:
        controls = [(args.order, p) for p in ((args.prime,) if args.prime is not None else CONTROLS[args.order])]
    try:
        results = [check_case(d, p, args.forms) for d, p in controls]
    except ValueError as error:
        parser.error(str(error))
    print(json.dumps(results, indent=2))


if __name__ == "__main__":
    main()
