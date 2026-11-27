"""Exact good-prime neighbours of integral ternary lattices."""
from collections import defaultdict
from itertools import product

from cypari2 import Pari

pari = Pari()


def gram_matrix(gram):
    return pari.matrix(3, 3, [entry for row in gram for entry in row])


def quadratic(gram, v):
    return sum(v[i] * gram[i, j] * v[j] for i, j in product(range(3), repeat=2))


def neighbour_basis(gram, p, point):
    """Return the basis of a p-neighbour in the parent coordinates."""
    value = int(quadratic(gram, point))
    if value % p:
        raise ValueError("The point is not isotropic.")
    normal = [sum(int(gram[i, j]) * point[j] for j in range(3)) for i in range(3)]
    pivot = next(i for i, x in enumerate(normal) if x % p)
    inverse = pow(normal[pivot], -1, p)
    lift = list(point)
    lift[pivot] += p * ((-value // p * pow(2 * normal[pivot], -1, p)) % p)
    columns = [[0, 0, 0] for _ in range(3)]
    columns[pivot][pivot] = p * p
    for j in range(3):
        if j != pivot:
            columns[j][j] = p
            columns[j][pivot] = p * (-normal[j] * inverse % p)
    return pari.mathnf(pari.Mat([pari.Col(v) for v in columns + [lift]])) / p


def neighbours(gram, p):
    """Yield the p+1 neighbours, one for each isotropic projective line."""
    points = [(1, y, z) for y, z in product(range(p), repeat=2)]
    points += [(0, 1, z) for z in range(p)] + [(0, 0, 1)]
    for point in points:
        if quadratic(gram, point) % p:
            continue
        basis = neighbour_basis(gram, p, point)
        target = basis.mattranspose() * gram * basis
        if any(target[i, j].type() != "t_INT" for i, j in product(range(3), repeat=2)):
            raise ArithmeticError("The lifted neighbour is not integral.")
        yield target


class Closure:
    """List isometry classes and record each neighbour with multiplicity."""

    def __init__(self, seeds, primes):
        self.forms = []
        self.buckets = defaultdict(list)
        self.exact = {}
        for q in seeds:
            self.classify(gram_matrix(q))
        self.nseeds = len(self.forms)
        if self.nseeds != len(seeds):
            raise ArithmeticError("Seed list contains an isometry duplicate.")
        self.rows = {p: [] for p in primes}
        index = 0
        while index < len(self.forms):
            for p in primes:
                row = [self.classify(q) for q in neighbours(self.forms[index], p)]
                if len(row) != p + 1:
                    raise ArithmeticError("Unexpected valency.")
                self.rows[p].append(row)
            index += 1

    def classify(self, gram):
        transform = gram.qflllgram()
        reduced = transform.mattranspose() * gram * transform
        key = tuple(int(reduced[i, j]) for i in range(3) for j in range(i, 3))
        if key in self.exact:
            return self.exact[key]
        shell = tuple(map(int, reduced.qfrep(32)))
        for i in self.buckets[shell]:
            if self.forms[i].qfisom(reduced):
                self.exact[key] = i
                return i
        i = len(self.forms)
        self.forms.append(reduced)
        self.exact[key] = i
        self.buckets[shell].append(i)
        return i
