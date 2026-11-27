"""Exact Hecke and modular-symbol certificates for the residual theta sums."""

import argparse
from collections import defaultdict
from math import lcm
import json
from pathlib import Path
import sys

from local_square import ThetaDifference, case_data, delay_census, local_square, pari, zero_bound
from neighbours import Closure


# Each entry gives the conductor and eigenvalues as polynomials in a_3.
# For the case of character discriminant 2860, the splitting operator is T(7²).
PACKETS = {
    232: [(232, (-1, 1), {3: 1, 5: 1, 7: 2, 11: 3, 13: -1, 17: 0, 19: 0, 23: 4, 31: 3})],
    760: [
        (152, (-1, 1), {3: 1, 7: 3, 11: 2, 13: 1, 17: -5, 23: -1, 31: 4}),
        (760, (-2, -2, 1), {3: (0, 1), 7: 0, 11: 2, 13: (2, -1),
                            17: (2, 2), 23: (4, -2), 31: 4}),
    ],
    520: [
        (104, (-1, 1), {3: 1, 7: 5, 11: -2, 17: -3, 19: -2, 23: 4, 31: -4}),
        (520, (-4, -2, 1), {3: (0, 1), 7: 0, 11: (4, -1), 17: 2,
                             19: (4, -1), 23: (0, -1), 31: (4, -3)}),
    ],
    2860: [
        (1430, (-4, 1), {3: 2, 7: 4, 17: -6, 19: 0, 23: 8, 31: -4}),
        (286, (3, 1), {3: 2, 7: -3, 17: -6, 19: 0, 23: 1, 31: 10}),
    ],
    28: [(1792, (-2, 1), {3: 2, 5: 2, 11: 0, 13: -2, 17: 6, 19: 6, 23: -4, 31: 8})],
}


def polynomial(value):
    return pari.Polrev(value) if isinstance(value, tuple) else pari.Pol(value)


def substitute(poly, value):
    return pari.subst(poly, pari("x"), value)


def known_packets(cases):
    common_levels = {}
    for _, entry in cases:
        char = entry["theta_character_conductor"]
        # Use the lcm since every case in a character group has the same level primes.
        common_levels[char] = lcm(common_levels.get(char, 1), entry["theta_level"] // 2)
    certified = {}
    for char, level in common_levels.items():
        specs = PACKETS[char]
        splitting_prime = 7 if char == 2860 else 3
        matches = [[] for _ in specs]
        shortest_fingerprint = [2 for _ in specs]
        for conductor in map(int, pari.divisors(level)):
            for form in pari.mfeigenbasis(pari.mfinit([conductor, 2], 0)):
                generator = pari.mfcoef(form, splitting_prime)
                for index, (_, factor, fingerprint) in enumerate(specs):
                    if substitute(polynomial(factor), generator):
                        continue
                    for position, (prime, expected) in enumerate(fingerprint.items(), 1):
                        if pari.mfcoef(form, prime) != substitute(polynomial(expected), generator):
                            shortest_fingerprint[index] = max(shortest_fingerprint[index], position)
                            break
                    else:
                        matches[index].append(conductor)
        for (conductor, factor, fingerprint), found, length in zip(specs, matches, shortest_fingerprint):
            if found != [conductor]:
                raise ArithmeticError(f"Fingerprint did not isolate conductor {conductor}: {found}")
            certified[char, factor] = (conductor, dict(list(fingerprint.items())[:length]), level)
    return certified


class WeightedTheta(ThetaDifference):
    def __init__(self, forms, weights):
        pairs = [(q, int(w)) for q, w in zip(forms, weights, strict=True) if w]
        if not pairs or sum(w for _, w in pairs):
            raise ArithmeticError("A nonzero cuspidal theta combination was expected.")
        super().__init__([q for q, _ in pairs])
        self.weights = [w for _, w in pairs]

    def through(self, bound):
        if bound > self.checked:
            self.coefficients = sum(w*q.qfrep(bound) for q, w in zip(self.forms, self.weights))
            self.checked = bound
        return self.coefficients


def neighbour_operator(seeds, prime, bound):
    closure = Closure([], [])
    labels = [[closure.classify(q) for q in forms] for forms in seeds]
    adjacency = []
    for gram in closure.forms:
        adjacency.append([closure.classify(q) for q in delay_census.neighbours(gram, prime)])
    theta = pari.Mat([pari.Col(q.qfrep(bound)) for q in closure.forms])
    return closure.forms, labels, adjacency, theta


def cyclic_span(vector, adjacency, theta):
    powers = [vector]
    coefficients = [theta*vector]
    while True:
        values = [0]*len(vector)
        for targets, weight in zip(adjacency, vector, strict=True):
            if weight:
                for target in targets:
                    values[target] += weight
        vector = pari.Col(values)
        coefficients.append(theta*vector)
        relations = pari.Mat(coefficients).matker()
        if len(relations):
            relation = relations[0] / relations[0][len(relations[0])-1]
            return pari.Polrev(list(relation)), pari.Mat(powers), pari.Mat(coefficients[:-1])
        powers.append(vector)


def project_block(minimal, factor):
    quotient = minimal / factor
    return pari.lift(pari.Mod(quotient, factor)**-1) * quotient


def coefficient_vector(poly, size):
    return pari.Col([pari.polcoef(poly, n) for n in range(size)])


def verify_fingerprint(theta, original, factor, projector, coefficient_basis, fingerprint, determinant):
    bound = zero_bound(theta.level)
    coefficients = theta.through(max(fingerprint)**2 * bound)
    for prime, eigenvalue in fingerprint.items():
        q_squared = prime*prime
        action = pari.lift(pari.Mod(polynomial(eigenvalue)*projector, original))
        expected = coefficient_basis*coefficient_vector(action, int(pari.poldegree(original)))
        for n in range(1, bound+1):
            term = coefficients[n*q_squared-1] + pari.kronecker(-determinant*n, prime)*coefficients[n-1]
            if n%q_squared == 0:
                term += prime*coefficients[n//q_squared-1]
            if term != expected[n-1]:
                raise ArithmeticError(f"Hecke fingerprint not certified at prime {prime} and norm {n}.")
    if not pari.polresultant(factor, polynomial(fingerprint[3])**2 - 16):
        raise ArithmeticError("The fingerprint may include a unary theta component.")


def central_symbol(conductor, split, factor, fingerprint, character):
    symbols = pari.msinit(conductor, 2, -1)
    new = pari.msnew(symbols)
    generator = pari.mshecke(symbols, split, new)
    identity = pari.matid(len(generator))

    def eval_matrix(value):
        poly = polynomial(value)
        if pari.poldegree(poly) <= 0:
            return pari.polcoef(poly, 0)*identity
        return substitute(poly, generator)

    vectors = pari.matker(eval_matrix(factor))
    for q, value in fingerprint.items():
        matrix = pari.mshecke(symbols, q, new)-eval_matrix(value)
        vectors = vectors*(matrix*vectors).matker()
    basis = new[0]*vectors
    if len(basis) != len(factor)-1:
        raise ArithmeticError("The odd modular-symbol packet has the wrong dimension.")

    disc = int(pari.coredisc(-character))
    infinity = pari("oo")
    values = [0]*len(basis)
    for a in range(1,abs(disc)):
        weight = pari.kronecker(disc,a)
        if weight:
            path = [infinity, pari(a)/abs(disc)]
            for j, value in enumerate(pari.mseval(symbols,basis,path)):
                values[j] += weight*value
    if any(values):
        raise ArithmeticError(f"The central symbol did not vanish at level {conductor}: {values}.")
    return {"conductor": conductor, "quadratic_twist": disc, "symbol_basis_dimension": len(basis),
            "exact_twisted_symbol": list(map(int, values))}


def analyse(cases):
    targets = known_packets(cases)
    groups = defaultdict(list)
    for index, entry in cases:
        _, theta, primes = case_data(entry)
        groups[entry["parent"][0]].append((index, theta, primes))
    cases_out = {}
    symbols_out = {}
    for determinant, entries in groups.items():
        level = entries[0][1].level
        char = entries[0][1].conductor
        split = 7 if char == 2860 else 3
        bound = zero_bound(level)
        forms, labels, adjacency, theta = neighbour_operator(
            [h.forms for _, h, _ in entries], split, bound
        )
        for (index, _, primes), indices in zip(entries, labels):
            v = [0]*len(forms)
            for j, weight in zip(indices, (1,1,-1,-1)):
                v[j] += weight
            minimal, basis, coefficient_basis = cyclic_span(pari.Col(v), adjacency, theta)
            factors, multiplicities = pari.factor(minimal)
            if any(m != 1 for m in multiplicities):
                raise ArithmeticError("Expected semisimple Hecke factors.")
            blocks = []
            for factor in factors:
                key = tuple(map(int, pari.Vecrev(factor)))
                projector = project_block(minimal, factor)
                raw = basis*coefficient_vector(projector, len(basis))
                content = pari.content(raw)
                h = WeightedTheta(forms, raw/content)
                coefficients = h.through(bound)
                witness = next((n for n in range(1,bound+1)
                                if coefficients[n-1] and all(local_square(n,p) for p in primes)),None)
                item = {"factor": key}
                if (char, key) in targets:
                    if witness is None:
                        raise ArithmeticError("The global block lost its simultaneous local-square witness.")
                    conductor, fingerprint, common_level = targets[char,key]
                    verify_fingerprint(h,minimal,factor,projector/content,coefficient_basis,fingerprint,determinant)
                    if (char,key) not in symbols_out:
                        symbols_out[char,key] = central_symbol(conductor,split,key,fingerprint,char)
                        symbols_out[char,key].update({"character_discriminant": char, "factor": key,
                                                     "fingerprint": fingerprint,
                                                     "newform_search_level": common_level})
                    item.update({"central_packet": conductor, "simultaneous_local_square_norm": witness})
                else:
                    obstruction = h.obstruction(2)
                    if witness is not None or "norm" in obstruction:
                        raise ArithmeticError("An unclassified nonlocal block remains.")
                    item["local_square_exclusion"] = obstruction
                blocks.append(item)
            cases_out[index] = {"case": index, "parent_determinant": determinant, "theta_level": level,
                                "splitting_prime": split, "neighbour_classes": len(forms),
                                "theta_sturm_bound": bound,
                                "minimal_polynomial": tuple(map(int,pari.Vecrev(minimal))),
                                "blocks": blocks}
    return {"cases": [cases_out[i] for i in sorted(cases_out)],
            "central_packets": sorted(symbols_out.values(), key=lambda row: row["conductor"])}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    config = parser.parse_args()
    pari.allocatemem(64_000_000,768_000_000,silent=True)
    data = json.loads(Path(__file__).with_name("local-square.json").read_text())
    cases = [(i,r) for i,r in enumerate(data)
             if "simultaneously_local_square" in r and "waldspurger_packet" not in r]
    result = json.dumps(analyse(cases), indent=2)
    if config.verify and json.loads(result) != json.loads(config.verify.read_text()):
        raise ArithmeticError("The residual-packet replay differs from its receipt.")
    print(result)


if __name__ == "__main__":
    main()
