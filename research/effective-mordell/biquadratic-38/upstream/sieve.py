"""Squareclass sieve over L = Q(zeta12) for the system A^2 = s^4+3, 7 D^2 = s^4+27, s in Q(i).

Over L,  s^2 + sqrt(-3) = delta * alpha^2  and  s^2 + 3 sqrt(-3) = gamma * beta^2  with
delta, gamma in L(S,2), S = primes above 2, 3, 7.  A and D lie in Q(i), so
delta * tau(delta) is a square and gamma * tau(gamma) is 7 times a square, where tau fixes i.
Eliminations use only necessary residue conditions, so every elimination is a proof.
"""
from itertools import product
from cypari2 import Pari

pari = Pari(); pari.allocatemem(2*10**9)
pari('z')
POL = pari('w^4-w^2+1')
L = pari.bnfinit(POL, 1)
M = lambda a: pari.Mod(pari(a), POL)
ZETA, I, R3 = M('w'), M('w^3'), M('2*w^2-1')
assert I**2 == -1 and R3**2 == -3


def tau(a):                          # fixes i = w^3, sends sqrt(-3) to -sqrt(-3)
    return pari.Mod(pari.subst(pari.lift(M(a)), 'w', pari('w^5')), POL)


assert tau(I) == I and tau(R3) == -R3


def is_square(a):
    a = M(a)
    return a == 0 or len(pari.nfroots(L, pari('z^2') - pari.lift(a))) > 0


S = [P for q in (2, 3, 7) for P in pari.idealprimedec(L, q)]
SUNITS = [M(g) for g in pari.bnfsunit(L, S)[0]]
BASIS = [ZETA, M('w-1')] + SUNITS        # F2-basis of L(S,2): dimension 6


def element(bits):
    out = M(1)
    for b, g in zip(bits, BASIS):
        if b:
            out *= g
    return out


CLASSES = [element(b) for b in product((0, 1), repeat=6)]
# Sanity: the 64 classes are distinct mod squares.
assert all(not is_square(CLASSES[a]*CLASSES[b]) for a in range(64) for b in range(a+1, 64))


def coeffs(a):
    a = pari.lift(M(a))
    return tuple(int(pari.polcoef(a, k, 'w')) for k in range(4))


def denominator(a):
    return int(pari.denominator(pari.content(pari.lift(M(a))))) if pari.lift(M(a)).type() == 't_POL' else int(pari.denominator(pari.lift(M(a))))


# ---------------------------------------------------------------- residue rings
def mul(x, y, m):
    """Product in Z[w]/(m, w^4-w^2+1)."""
    c = [0]*7
    for i, a in enumerate(x):
        if a:
            for j, b in enumerate(y):
                c[i+j] += a*b
    for k in (6, 5, 4):              # w^k = w^(k-2) - w^(k-4)
        c[k-2] += c[k]; c[k-4] -= c[k]; c[k] = 0
    return tuple(v % m for v in c[:4])


def square_set(m):
    return {mul(x, x, m) for x in product(range(m), repeat=4)}


def gaussian(a, b):                  # a + b i as element of Z[w]; i = w^3
    return (a, 0, 0, b)


def add(x, y, m):
    return tuple((u+v) % m for u, v in zip(x, y))


def global_filter(pairs):
    out = []
    for d, g in pairs:
        if is_square(d*tau(d)) and is_square(7*g*tau(g)):
            out.append((d, g))
    return out


def ramified_test(pairs, q, m):
    """Necessary conditions modulo m at the ramified prime q (2 or 3)."""
    sq = square_set(m)
    r3 = coeffs(R3); r33 = tuple(3*c for c in r3)
    seven = (7, 0, 0, 0)
    reps = [(coeffs(d), coeffs(g)) for d, g in pairs]
    unit = (lambda a, b: (a+b) % 2 == 1) if q == 2 else (lambda a, b: a % 3 or b % 3)
    found = [False]*len(pairs)
    for ua, ub, wa, wb in product(range(m), repeat=4):
        if not (unit(ua, ub) or unit(wa, wb)):
            continue
        u, w = gaussian(ua, ub), gaussian(wa, wb)
        u2, w2 = mul(u, u, m), mul(w, w, m)
        u4, w4 = mul(u2, u2, m), mul(w2, w2, m)
        if add(u4, mul((3, 0, 0, 0), w4, m), m) not in sq:
            continue
        if mul(seven, add(u4, mul((27, 0, 0, 0), w4, m), m), m) not in sq:
            continue
        a = add(u2, mul(r3, w2, m), m)
        b = add(u2, mul(r33, w2, m), m)
        for k, (dc, gc) in enumerate(reps):
            if not found[k] and mul(dc, a, m) in sq and mul(gc, b, m) in sq:
                found[k] = True
    return [pr for pr, f in zip(pairs, found) if f]


def residue_test(pairs, p):
    """Necessary conditions in the residue fields of L above an odd prime p, one prime of Q(i) at a time.

    The value s lies in the projective line over the residue field of each prime q of Q(i)
    independently, so a split p is treated as two separate tests.
    """
    primes = pari.idealprimedec(L, p)
    issq = lambda x: x == 0 or bool(pari.issquare(x))
    if p % 4 == 3:
        groups = [([pari.nfmodprinit(L, P) for P in primes],
                   [(M(a + b*I), M(1)) for a in range(p) for b in range(p)] + [(M(1), M(0))])]
    else:
        r = next(c for c in range(p) if (c*c + 1) % p == 0)
        groups = []
        for root in (r, p - r):
            above = [P for P in primes if pari.nfmodpr(L, pari.lift(I - root), pari.nfmodprinit(L, P)) == 0]
            groups.append(([pari.nfmodprinit(L, P) for P in above],
                           [(M(c), M(1)) for c in range(p)] + [(M(1), M(0))]))
    red = lambda a, mp: pari.nfmodpr(L, pari.lift(M(a)), mp)
    survivors = pairs
    for modprs, pts in groups:
        found = [False]*len(survivors)
        for u, w in pts:
            s4, w4 = u**4, w**4
            if not all(issq(red(s4 + 3*w4, mp)) and issq(red(7*(s4 + 27*w4), mp)) for mp in modprs):
                continue
            a = u**2 + R3*w**2; b = u**2 + 3*R3*w**2
            for k, (d, g) in enumerate(survivors):
                if not found[k] and all(issq(red(d*a, mp)) and issq(red(g*b, mp)) for mp in modprs):
                    found[k] = True
        survivors = [pr for pr, f in zip(survivors, found) if f]
    return survivors


if __name__ == '__main__':
    pairs = global_filter([(d, g) for d in CLASSES for g in CLASSES])
    print('global norm conditions:', len(pairs), flush=True)
    for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37, 41, 43):
        pairs = residue_test(pairs, p)
        print('after residue fields above', p, len(pairs), flush=True)
    for q, m in ((3, 9), (2, 8), (2, 16)):
        pairs = ramified_test(pairs, q, m)
        print('after ramified', q, 'mod', m, len(pairs), flush=True)
    for d, g in pairs:
        print(pari.lift(d), '|', pari.lift(g))
