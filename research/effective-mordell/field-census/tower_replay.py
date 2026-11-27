"""Reproduce exact local, torsion and point checks for the fixed curve y^2=x^6-1."""
from itertools import combinations, product
from math import prod
from cypari2 import Pari

pari = Pari()


def require(statement, message):
    if not statement:
        raise ArithmeticError(message)


def gaussian_residues():
    """All elements of Z[i]/4. The ideal (4) is supported at one prime."""
    def multiply(a, b):
        return ((a[0]*b[0] - a[1]*b[1]) % 4, (a[0]*b[1] + a[1]*b[0]) % 4)

    def total(*values):
        return tuple(sum(coordinate) % 4 for coordinate in zip(*values))

    residues = list(product(range(4), repeat=2))
    squares = {multiply(t, t) for t in residues}
    for delta in (1, -1):
        for a, b in ((3*delta, -3*delta**2), (-6*delta, 3*delta**2)):
            for first, last in (((0, 1), (0, b)), ((0, b), (0, 1))):
                for t2 in squares:
                    value = total(multiply(first, multiply(t2, t2)), multiply((a, 0), t2), last)
                    require(value not in squares, 'dyadic isogeny obstruction')
    for c in ((0, 1), (0, 3)):
        for t2 in squares:
            value = multiply(c, total(multiply(t2, t2), t2, (1, 0)))
            require(value not in squares, 'quartic dyadic obstruction')


def finite_targets():
    x, w = pari('x'), pari('w')
    field = pari.nfinit(w**2 + 1)
    for parameter, size, expected in ((1, 8, {'0', 'w', '-w'}), (3, 4, {'0', '1', '-1'})):
        data = [0, parameter, 0, parameter**2, 0]
        rational = pari.ellinit(data)
        for target in (rational, rational.elltwist(-4)):
            require(list(target.ellrank()[:2]) == [0, 0], 'Gaussian cover target rank')
        curve = pari.ellinit(data, field)
        order, structure, gens = curve.elltors()
        require(int(order) == size, 'Gaussian cover target torsion order')
        points = [pari([0])]
        for n, g in zip(structure, gens):
            points = [curve.elladd(P, curve.ellmul(g, j)) for P in points for j in range(int(n))]
        require(len(set(map(str, points))) == size, 'missing torsion points')
        abscissas = {str(pari.lift(r)) for P in points if len(P) == 2
                     for r in pari.nfroots(field, x**2-P[0]/parameter)}
        require(abscissas == expected, 'Gaussian cover target fibres')
    require(pari.polresultant(x*x-3, pari('y')**2-9-6*x, x) ==
            pari('y')**4-18*pari('y')**2-27, '4-torsion norm')
    require(len(pari.nfroots(field, x*x+27)) == 0, 'norm accidentally a Gaussian square')
    for prime in (11, 23):
        require(pari.ellinit([0, -1], prime).ellcard() == prime+1, 'vanishing trace')


def restricted_twists():
    primes = (5, 17, 29, 41, 53, 89)
    count = 0
    for n in range(4):
        for subset in combinations(primes, n):
            for factor in (2, -6):
                d = factor*prod(subset)
                require(list(pari.ellinit([0, -d**3]).ellrank()[:2]) == [0, 0], 'even twist')
                count += 1
    return count


def coordinates():
    w = pari('w')
    phi8, phi12 = w**4+1, w**4-w**2+1
    for poly in (phi8, phi12):
        require(pari.polisirreducible(poly), 'cyclotomic defining polynomial')
    for poly, expected in ((phi8, 18), (phi12, 10)):
        nf = pari.nfinit(poly)
        x = pari('x')
        points = {(str(pari.lift(a)), str(pari.lift(b)))
                  for a in list(pari.nfroots(nf, x**6-1)) + list(pari.nfroots(nf, x*x+1))
                         + list(pari.nfroots(nf, 2*x*x+1)) + list(pari.nfroots(nf, x*x+2)) + [pari(0)]
                  for b in pari.nfroots(nf, x*x-a**6+1)}
        require(len(points)+2 == expected, 'cyclotomic point list')
    for field, val in (('w^4+10*w^2+1', 20),):
        nf = pari.nfinit(field)  # w = sqrt(-2)+sqrt(-3).
        xs = list(pari.nfroots(nf, pari('x')**6-1)) + list(pari.nfroots(nf, 8*pari('x')**6+1))
        pts = {(str(pari.lift(a)), str(pari.lift(b))) for a in xs
               for b in pari.nfroots(nf, pari('x')**2 - a**6 + 1)}
        require(len(pts)+2 == val, 'negative even tower point list')


if __name__ == '__main__':
    gaussian_residues()
    finite_targets()
    print(restricted_twists(), 'rational even-character twists verified.')
    coordinates()
    print('Dyadic certificates, rank-zero targets, and point lists of sizes 10, 18, 20 verified.')
