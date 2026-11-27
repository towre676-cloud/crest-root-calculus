"""C6 rank-zero certificates over number fields with arbitrary ideal class groups.

Put this file next to the previous isogeny_over_field.py module.
"""
import argparse
import json
from collections import deque
from itertools import product
from math import inf, prod

from isogeny_over_field import lift_torsion, pari


def half_classes(order, value):
    """Solutions of 2z = -value modulo order."""
    if order % 2:
        return [(-value * pow(2, -1, order)) % order]
    if value % 2:
        return []
    half = order // 2
    root = (-value // 2) % half
    return [root, root + half]


def global_squareclasses(field, budget):
    primes = list(pari.idealprimedec(field, 3))
    orders = list(map(int, field.bnf_get_cyc()))
    basis = list(field.bnf_get_gen())
    units = [field.bnf_get_tu()[1], *field.bnf_get_fu()]
    unit_products = [pari.Mod(1, field.nf_get_pol())]
    if 2 ** len(units) > budget:
        return None, None
    for unit in units:
        unit_products += [u * unit for u in unit_products]

    one = pari.idealhnf(field, 1)
    prime_classes = [list(map(int, pari.bnfisprincipal(field, pr, 0))) for pr in primes]
    seeds = []
    for parities in product((0, 1), repeat=len(primes)):
        target = [sum(e * c[i] for e, c in zip(parities, prime_classes)) for i in range(len(orders))]
        halves = [half_classes(h, v) for h, v in zip(orders, target)]
        for coordinates in product(*halves):
            if (len(seeds) + 1) * len(unit_products) > budget:
                return None, None
            seeds.append((parities, coordinates))

    reps = []
    for parities, coordinates in seeds:
        ideal = one
        for pr, exponent in zip(primes, parities):
            if exponent:
                ideal = pari.idealmul(field, ideal, pr)
        for gen, exponent in zip(basis, coordinates):
            if exponent:
                ideal = pari.idealmul(field, ideal, pari.idealpow(field, gen, 2 * exponent))
        remainder, principal = pari.bnfisprincipal(field, ideal, 1)
        if any(remainder):
            raise ArithmeticError("Class-group relation failed to produce a principal ideal")
        generator = pari.nfbasistoalg(field, principal)
        if pari.idealhnf(field, generator) != pari.idealhnf(field, ideal):
            raise ArithmeticError("Wrong principal ideal generator")
        reps.extend(generator * unit for unit in unit_products)
    description = dict(cyclic_invariants=orders, prime_classes_above_3=prime_classes,
                       squareclasses=len(reps), parity_class_lifts=len(seeds),
                       unit_squareclasses=len(unit_products))
    return reps, description


class LocalDescent:
    """Test smooth isogeny quartics over a completion by exhaustive residue disks."""

    def __init__(self, field, pr):
        self.field = field
        self.pr = pr
        self.uniformizer = pari.nfbasistoalg(field, pr[1])
        if self.valuation(self.uniformizer) != 1:
            raise ArithmeticError("Wrong prime-ideal uniformizer")
        self.twice_dyadic_index = 2 * self.valuation(pari(2))
        self.digits = None
        self.unit_precision = pari.matrix(1, 2, [pr, self.twice_dyadic_index + 1])

    def residue_digits(self):
        if self.digits is None:
            hnf = pari.idealhnf(self.field, self.pr)
            diag = [int(col[i]) for i, col in enumerate(hnf)]
            self.digits = [pari.nfbasistoalg(self.field, pari.Col(coords))
                           for coords in product(*(range(h) for h in diag))]
            if len(self.digits) != pari.idealnorm(self.field, self.pr):
                raise ArithmeticError("Wrong number of residue representatives")
        return self.digits

    def valuation(self, value):
        if value == 0:
            return inf
        return int(pari.nfeltval(self.field, value, self.pr))

    def is_square(self, value):
        if value == 0:
            return True
        v = self.valuation(value)
        if v % 2:
            return False
        unit = value / self.uniformizer**v
        coords = pari.nfalgtobasis(self.field, unit)
        if any(c.denominator() != 1 for c in coords):
            # Get a globally integral lift even when denominators involve a
            # different prime above the same rational prime.
            integral = pari.idealchinese(self.field, self.unit_precision, [unit])
            lifted = pari.nfbasistoalg(self.field, integral)
            if self.valuation(lifted - unit) <= self.twice_dyadic_index:
                raise ArithmeticError("Unit lift failed its squareclass precision")
            unit = lifted
        return bool(pari.nfislocalpower(self.field, self.pr, unit, 2))

    def covers(self, representative, a, b):
        # Normalize the first coefficient up to a local square; only local
        # integrality is needed. PARI accepts fractions integral at this prime.
        valuation = self.valuation(representative)
        first = representative / self.uniformizer ** (2 * (valuation // 2))
        last = b / first
        if self.valuation(first) < 0 or self.valuation(last) < 0:
            raise ArithmeticError("Nonintegral local model")
        if self.is_square(first) or self.is_square(last):
            return True
        zero = pari(0)
        queue = deque([(first, last, zero, 0), (last, first, zero, 1)])
        while queue:
            lead, constant, t, depth = queue.popleft()
            t2 = t*t
            value = (lead*t2 + a)*t2 + constant
            if self.is_square(value):
                return True
            v = self.valuation(value)
            linear = (4*lead*t2 + 2*a)*t
            derivative_value = self.valuation(linear)
            if v > 2*derivative_value:
                return True       # A simple root with ordinate zero by Hensel.
            lower_bound = min(depth + derivative_value,
                              2*depth + self.valuation(6*lead*t2 + a),
                              3*depth + self.valuation(4*lead*t),
                              4*depth + self.valuation(lead))
            if lower_bound > v + self.twice_dyadic_index:
                continue          # All values have the same nonsquare class.
            place = self.uniformizer**depth
            for digit in self.residue_digits():
                queue.append((lead, constant, t + digit*place, depth + 1))
        return False


def local_partition(local, global_reps):
    """Project onto the distinct squareclasses in one completion."""
    distinct = []
    classes = []
    for e in global_reps:
        for i, d in enumerate(distinct):
            if local.is_square(e/d):
                classes.append(i)
                break
        else:
            classes.append(len(distinct))
            distinct.append(e)
    return distinct, classes


def classify(polynomial, max_classes=None):
    if (not isinstance(polynomial, list) or len(polynomial) < 2 or polynomial[-1] != 1
            or any(type(c) != int for c in polynomial)):
        raise ValueError("Supply a monic integer polynomial as a list in ascending order")
    f = pari.Polrev(polynomial, pari('w'))
    if not pari.polisirreducible(f):
        raise ValueError("Reducible field polynomial")
    field = pari.bnfinit(f, 1)
    if pari.bnfcertify(field) != 1:
        raise ArithmeticError("Unable to certify the class group and units")
    result = dict(field=str(f), class_number=int(field.bnf_get_no()), status="partial")
    classes, summary = global_squareclasses(field, inf if max_classes is None else max_classes)
    if classes is None:
        return dict(result, reason="Squareclass budget exceeded")
    result.update(summary)
    real = [i for i, e in enumerate(classes)
            if all(int(v) > 0 for v in pari.nfeltsign(field, e))]
    tests = {sign: [real.copy(), list(range(len(classes)))] for sign in (-1, 1)}
    for prime in (3, 2):
        for pr in pari.idealprimedec(field, prime):
            local = LocalDescent(field, pr)
            distinct, membership = local_partition(local, classes)
            for sign, images in tests.items():
                a = -3 * sign
                for j, (aj, bj) in enumerate(((a, 3), (-2*a, -3))):
                    active = set(membership[t] for t in images[j])
                    compatible = {i for i in active if local.covers(distinct[i], aj, bj)}
                    images[j] = [t for t in images[j] if membership[t] in compatible]
    sizes = {str(sign): [len(xs) for xs in images] for sign, images in tests.items()}
    for bounds in sizes.values():
        if prod(bounds) < 4 or any(n == 0 or n & (n - 1) for n in bounds):
            raise ArithmeticError("Local images have invalid subgroup sizes")
    result['selmer_sizes'] = sizes
    result['rank_upper_bounds'] = {s: prod(v).bit_length()-3 for s, v in sizes.items()}
    for sign in (-1, 1):
        if prod(sizes[str(sign)]) == 4:
            return dict(result, status="complete", sign=sign,
                        verified_affine_points=lift_torsion(field, sign), rational_infinities=2)
    return dict(result, reason="Both exact 2-isogeny descent rank upper bounds are positive")


if __name__ == "__main__":
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("polynomial")
    cli.add_argument("--max-classes", type=int)
    args = cli.parse_args()
    print(json.dumps(classify(json.loads(args.polynomial), args.max_classes), indent=2))
