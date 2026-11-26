"""Inverse problem for quotient spectra in CM isogeny classes E x E, E = C/O_K.

For each prime level l and every admissible gluing psi (det psi = -1 in GL_2(F_l)) the refined Humbert
form q_psi is computed exactly (gluing.py).  Products (q represents 1: theta reducible) are discarded.
Forms are classified
  (i)   up to isometry (PARI qfisom: exact),
  (ii)  up to genus: the even lattice 2q has genus determined by its signature and discriminant
        quadratic form (Nikulin, Thm 1.9.4 / Cor 1.9.3); isomorphism of finite quadratic forms is
        decided exactly by backtracking on Smith generators,
  (iii) up to square-equivalence r_q(n^2) = r_q'(n^2) for all n.  For q, q' in one genus,
        theta_q - theta_q' is a weight-3/2 cusp form on Gamma_0(N) (N = level of 2q); its Shimura
        lift with t = 1 is a weight-2 modular form on Gamma_0(N/2) with coefficients
        A(n) = sum_{d | n} chi_1(d) (r_q - r_q')(n^2/d^2) (Shimura 1973, Niwa 1975), so equality of
        r(m^2) for all m <= B, B the Sturm bound (1/6) N prod_{p|N} (1 + 1/p), proves equality for
        every m.  Forms in different genera are separated by the first square where they differ.
"""
import itertools
import json
import math
import sys
import time

import cypari2

from gluing import ns_and_q, reduce_form, isometric

pari = cypari2.Pari()
pari.allocatemem(6 * 10 ** 9)

CM = {-3: "Q(sqrt-3)", -4: "Q(i)", -7: "Q(sqrt-7)", -8: "Q(sqrt-2)", -11: "Q(sqrt-11)", -19: "Q(sqrt-19)"}


def kind(d, l):
    s = pari.kronecker(d, l)
    return {1: "split", -1: "inert", 0: "ramified"}[int(s)]


def is_scalar(d, l, psi, c=None):
    """psi = multiplication by some lambda in O / l on the basis (1, w)?"""
    c = (d * d - d) // 4 if c is None else c
    for p, q in itertools.product(range(l), repeat=2):
        # lambda = p + q w : 1 -> (p, q), w -> lambda w = p w + q (d w - c) = (-q c, p + q d)
        M = [[p % l, (-q * c) % l], [q % l, (p + q * d) % l]]
        if M == [[x % l for x in r] for r in psi]:
            return True
    return False


def mat(q):
    return pari.matrix(3, 3, [x for r in q for x in r])


def level(q):
    A = mat(q) * 2
    Ai = pari.matsolve(A, pari.matid(3))
    N = 1
    while True:
        M = Ai * N
        if all(str(pari.type(M[i, j])) == "t_INT" for i in range(3) for j in range(3)) and \
                all(int(M[i, i]) % 2 == 0 for i in range(3)):
            return N
        N += 1


def disc_form(q):
    """(invariants, generators' Gram data) of the discriminant form of the even lattice 2q."""
    A = mat(q) * 2
    S = pari.matsnf(A, 1)
    U, V, D = S[0], S[1], S[2]
    # A^# / A : generators are columns of A^{-1} V... use: L^#/L ~ Z^3 / A Z^3 via x -> A^{-1} x
    Ai = pari.matsolve(A, pari.matid(3))
    invs, gens = [], []
    Ui = pari.matsolve(U, pari.matid(3))
    for i in range(3):
        e = int(D[i, i])
        if e > 1:
            invs.append(e)
            gens.append(Ai * Ui[:, i] if False else Ai * pari.matrix(3, 1, [Ui[k, i] for k in range(3)]))
    return invs, gens, A


def qval(A, x):
    from fractions import Fraction
    v = pari.mattranspose(x) * A * x
    return v[0, 0]


def disc_iso(q1, q2):
    i1, g1, A1 = disc_form(q1)
    i2, g2, A2 = disc_form(q2)
    if sorted(i1) != sorted(i2):
        return False
    # enumerate all elements of group 2 as rational vectors mod L
    def elements(invs, gens):
        out = []
        for co in itertools.product(*[range(e) for e in invs]):
            x = pari.matrix(3, 1, [0, 0, 0])
            for c, g in zip(co, gens):
                x = x + g * c
            out.append(x)
        return out
    E2 = elements(i2, g2)
    fq = lambda A, x: pari.frac(qval(A, x) / 2) * 2          # q(x) mod 2Z  (x^T A x mod 2)
    fb = lambda A, x, y: pari.frac((pari.mattranspose(x) * A * y)[0, 0])
    def order(x, A):
        k = 1
        while not all(str(pari.type(c)) == "t_INT" for c in (pari.matsolve(pari.matid(3), x * k)).list()):
            k += 1
        return k
    tgt = [[y for y in E2 if order(y, A2) == e and fq(A2, y) == fq(A1, g)] for e, g in zip(i1, g1)]
    n = len(g1)

    def bt(k, imgs):
        if k == n:
            # must be injective/bijective: group sizes equal, check kernel trivial by building images
            seen = set()
            for co in itertools.product(*[range(e) for e in i1]):
                y = pari.matrix(3, 1, [0, 0, 0])
                for c, im in zip(co, imgs):
                    y = y + im * c
                seen.add(tuple(str(pari.frac(c)) for c in y.list()))
            return len(seen) == math.prod(i1)
        for y in tgt[k]:
            if all(fb(A2, y, imgs[j]) == fb(A1, g1[k], g1[j]) for j in range(k)):
                if bt(k + 1, imgs + [y]):
                    return True
        return False
    return bt(0, [])


def square_counts(q, B):
    """r_q(m^2) for m = 1..B (number of vectors, both signs)."""
    r = pari.qfrep(mat(q), B * B)
    return [2 * int(r[m * m - 1]) for m in range(1, B + 1)]


def sturm(N):
    ps = [int(p) for p in pari.factor(N)[0]]
    idx = N
    for p in ps:
        idx = idx * (p + 1) // p
    return math.ceil(2 * idx / 12)


def run(ds, ls, Ncmp=120):
    out = {}
    for d in ds:
        classes = []                      # one representative per isometry class across all levels
        for l in ls:
            t = time.time()
            seen = []
            for a, b, c, e in itertools.product(range(l), repeat=4):
                if (a * e - b * c + 1) % l:
                    continue
                psi = [[a, b], [c, e]]
                r = ns_and_q(d, l, psi)
                if r is None:
                    continue
                q = reduce_form(r["q"])
                if q[0][0] == 1 or min(q[i][i] for i in range(3)) == 1:
                    continue
                if pari.qfminim(mat(q), 1, 1)[0] > 0:    # represents 1: product polarisation
                    continue
                sc = is_scalar(d, l, psi)
                for cl in seen:
                    if isometric(cl["q"], q):
                        cl["gluings"] += 1
                        cl["scalar"] = cl["scalar"] or sc
                        break
                else:
                    seen.append({"q": q, "level": l, "gluings": 1, "scalar": sc,
                                 "det": int(pari.matdet(mat(q)))})
            for cl in seen:
                cl["kind"] = kind(d, l)
            classes += seen
            print(d, l, kind(d, l), "classes", len(seen), [(c["q"], c["det"], c["gluings"], c["scalar"]) for c in seen],
                  round(time.time() - t, 1), flush=True)
        out[d] = classes
    return out


if __name__ == "__main__":
    ds = [int(x) for x in sys.argv[1].split(",")]
    ls = [int(x) for x in sys.argv[2].split(",")]
    res = run(ds, ls)
    json.dump({str(k): v for k, v in res.items()}, open(f"classes_{sys.argv[1]}.json", "w"), indent=1)
