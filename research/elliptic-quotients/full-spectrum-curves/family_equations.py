"""Explicit determinantal models for the opposite-secant family."""

import argparse
from fractions import Fraction
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path[:0] = [str(HERE.parent / name) for name in
                ("twin-curves", "inverse-spectrum", "opposite-secant-spectra/delay")]
from flint import acb, arb, ctx
from gluing import E_of, K, gram, isometric, q_from_ns
from local_square import delay_census
from polarised import ns_general, pari


def children(t):
    if t < 0 or t % 5 == 3:
        raise ValueError("Require t >= 0 and t != 3 modulo 5.")
    parent = ((5+36*t, 1-36*t, 2), (1-36*t, 5+36*t, 2), (2, 2, 8))
    return {name: delay_census.sublattice(parent, normal, 3)
            for name, normal in (("A", (1, 2, 2)), ("B", (1, 0, 1)))}


def product_humbert(t, a, b, n, m, k):
    """Use the full integral NS lattice of C/O x C/<a,b+sqrt(-N)>."""
    N = 90 * (1+18*t)
    if (not 0 <= b < a or (b*b+N) % a or gcd(a, 2*b, (b*b+N)//a) != 1
            or min(n, m, k) < 1 or n*m-a*k*k != 1):
        raise ValueError("Invalid primitive ideal or principal polarization parameters.")
    field = K(0, N)
    basis = [tuple(map(Fraction, row)) for row in
             ((1, 0, 0, 0), (0, 1, 0, 0), (0, 0, a, 0), (0, 0, b, 1))]
    H = (Fraction(n), Fraction(m, a), Fraction(k), Fraction(0))
    alternating = [[E_of(field, H, u, v) for v in basis] for u in basis]
    assert all(z.denominator == 1 for row in alternating for z in row)
    assert pari.matrix(4, 4, [int(z) for row in alternating for z in row]).matdet() == 1
    ns, intersection = ns_general(field, basis)
    vectors = pari.matrix(4, 4, [pari(str(ns[j][i])) for i, j in product(range(4), repeat=2)])
    theta = vectors.matsolve(pari.Col([pari(str(v)) for v in H]))
    assert all(v.denominator() == 1 for v in theta)
    return q_from_ns(intersection, list(map(int, theta)))


def certificate(t, name, parameters):
    a, b, n, m, k = (parameters[key] for key in ("a", "b", "n", "m", "k"))
    targets = children(t)
    q = product_humbert(t, a, b, n, m, k)
    target = targets[name]
    if not isometric(q, target) or isometric(q, targets["B" if name == "A" else "A"]):
        raise ArithmeticError("The polarization does not isolate its target Humbert form.")
    assert gram(target).matdet() == 1440*(1+18*t)
    return {"t": t, "curve": name, **parameters,
            "cm_discriminant": -360*(1+18*t), "humbert_gram": target,
            "section_matrix_columns": (n+k)*(m+a*k),
            "torsion_labels_in_basis_1_r": [[b, 1, a], [1, 0, a]]}


def isolate_torsion(record):
    """Give exact dyadic boxes for the normalized division-polynomial roots."""
    ctx.prec = 384
    a, b = record["a"], record["b"]
    N = -record["cm_discriminant"]//4
    tau = acb(0, arb(N).sqrt())
    g2, g3 = tau.elliptic_invariants()
    scale = 2*g2/(3*g3)
    bits = 80

    def coordinate(u, v):
        return scale*((u+v*tau)/a).elliptic_p(tau)

    def box(point):
        centres = [int((component.mid()*2**bits + arb(1)/2).floor().unique_fmpz())
                   for component in (point.real, point.imag)]
        enclosure = acb(*(arb(arb(c)/2**bits, arb(1)/2**bits) for c in centres))
        assert enclosure.contains(point)
        return centres, enclosure

    centres, boxes = zip(*(box(coordinate(u, v)) for u, v in ((b, 1), (1, 0))), strict=True)
    counts = [0, 0]
    for v in range(a//2+1):
        # The x-coordinate identifies a point with its negative.
        if v == 0:
            available = range(1, a//2+1)
        elif 2*v == a:
            available = range(a//2+1)
        else:
            available = range(a)
        for u in available:
            z = coordinate(u, v)
            for j, enclosure in enumerate(boxes):
                counts[j] += z.overlaps(enclosure)
    assert counts == [1, 1], "A root box is not isolated among the division points."
    return {"denominator_exponent": bits, "centres": list(centres)}


def search(t):
    """Find a product presentation by expanding the polarization degree bound."""
    targets = children(t)
    matrices = {name: gram(q) for name, q in targets.items()}
    found = {}
    last_bound = 0
    bound = 16
    N = 90*(1+18*t)
    while len(found) != 2:
        counts = matrices["A"].qfrep(bound**2)
        degrees = [d for d in range(1, bound+1) if counts[d*d-1]]
        for n, m in product(degrees, repeat=2):
            if max(n, m) <= last_bound:
                continue
            norm = n*m-1
            for k in range(1, isqrt(norm)+1):
                if norm % (k*k):
                    continue
                a = norm//(k*k)
                # Square roots of -N modulo a, using prime-power roots and CRT.
                roots = [0]
                modulus = 1
                primes, powers = pari.factor(a)
                for p, e in zip(primes, powers, strict=True):
                    q = int(p**e)
                    local = square_roots(-N, int(p), int(e))
                    roots = [(x+modulus*((y-x)*pow(modulus, -1, q) % q))
                             for x in roots for y in local]
                    modulus *= q
                for b in roots:
                    if gcd(a, 2*b, (b*b+N)//a) != 1:
                        continue
                    q = product_humbert(t, a, b, n, m, k)
                    for name in targets.keys() - found.keys():
                        if isometric(q, targets[name]):
                            found[name] = certificate(t, name, dict(a=a, b=b, n=n, m=m, k=k))
                    if len(found) == 2:
                        return found
        last_bound, bound = bound, 2*bound


def square_roots(value, p, e):
    x = pari("x")
    roots = [int(pari.lift(r)) for r in pari.polrootsmod(x*x-value, p)]
    power = p
    for _ in range(1, e):
        roots = [r+j*power for r in roots for j in range(p)
                 if ((r+j*power)**2-value) % (p*power) == 0]
        power *= p
    return roots



def velu_graph(first, p, q, degree):
    """Build the dual graph isogeny from the two labelled division points."""
    x, y = pari("x"), pari("y")
    source = pari.ellinit(first)
    factors = [int(prime) for prime in pari.factor(degree)[0]]

    def check_order(curve, point):
        if pari.ellmul(curve, point, degree) != pari([0]) or any(
            pari.ellmul(curve, point, degree//prime) == pari([0]) for prime in factors
        ):
            raise ArithmeticError("A labelled point has the wrong order.")

    check_order(source, p)
    check_order(source, q)
    quotient, forward = pari.ellisogeny(source, p, 0, x, y)
    image = pari.ellisogenyapply(forward, q)
    check_order(pari.ellinit(quotient), image)
    back, (U, VY, D) = pari.ellisogeny(pari.ellinit(quotient), image, 0, x, y)
    if any(quotient[:3]) or any(back[:3]) or back[3] != first[0]*degree**4 or back[4] != first[1]*degree**6:
        raise ArithmeticError("The second Vélu map has the wrong scaled target.")
    return tuple(quotient[3:5]), (U/degree**2, pari.polcoef(VY, 1, y)/degree**3, D)

def elliptic_basis(bound):
    """All (power of x, power of y) with pole order at most bound at O."""
    return [(i, w) for w in range(2) for i in range((bound-3*w)//2+1)]


def section_matrix(first, second, h, n, m, k, degree):
    """Coefficient matrix for vanishing to order k along P = h(Q).

    Inputs are short Weierstrass [A,B] for each factor and [U,V,D] from
    Vélu, where h(x,y)=(U/D², V*y/D³). The columns are products of
    Weierstrass monomials on the factors, in the returned basis order.
    """
    x = pari("x")
    A, B = first
    AA, BB = second
    f, g = x**3+A*x+B, x**3+AA*x+BB
    U, V, D = h
    X, Y = U/D**2, V/D**3
    if Y**2*g != X**3+A*X+B:
        raise ArithmeticError("The supplied graph map fails the curve equation.")
    b1, b2 = elliptic_basis(n+k), elliptic_basis(m+degree*k)
    transformed = [(x**i, 0) if w == 0 else (0, x**i) for i, w in b1]
    columns = [[] for _ in range(len(b1)*len(b2))]
    for _ in range(k):
        restrictions = [(pari.subst(u, x, X), Y*pari.subst(v, x, X)) for u, v in transformed]
        terms = []
        for (u, v), (j, w) in product(restrictions, b2):
            terms.append((u*x**j, v*x**j) if w == 0 else (v*g*x**j, u*x**j))
        denominator = 1
        for term in terms:
            for z in term:
                denominator = pari.lcm(denominator, pari.denominator(z, x))
        polynomials = [tuple(pari.simplify(denominator*z) for z in term) for term in terms]
        degrees = [max(int(pari.poldegree(family[i], x)) for family in polynomials if family[i])
                   if any(family[i] for family in polynomials) else -1 for i in (0, 1)]
        for column, polynomials_by_parity in zip(columns, polynomials, strict=True):
            for f0, d in zip(polynomials_by_parity, degrees, strict=True):
                column.extend(pari.polcoef(f0, j, x) for j in range(d+1))
        transformed = [(2*f*pari.deriv(v, x)+pari.deriv(f, x)*v, 2*pari.deriv(u, x))
                       for u, v in transformed]
    return pari.Mat([pari.Col(column) for column in columns]), list(product(b1, b2))


def check_control():
    """Verify both one jet and multiple jets over exact coefficient fields."""
    x, y = pari("x"), pari("y")
    E = pari.ellinit([pari.Mod(1, 101), pari.Mod(13, 101)])
    P = pari.ellmul(E, [pari.Mod(37, 101), pari.Mod(100, 101)], 9)
    E2 = pari.ellinit([pari(1), pari(0)])
    for source, point, n, m, k, degree in ((E, P, 3, 4, 1, 11), (E2, [0, 0], 3, 3, 2, 2)):
        assert pari.ellorder(source, point) == degree
        target, (U, VY, D) = pari.ellisogeny(source, point)
        assert not any(target[:3])
        V = pari.polcoef(VY, 1, y)
        matrix, basis = section_matrix(tuple(target[3:5]), tuple(source[3:5]), (U, V, D), n, m, k, degree)
        kernel = matrix.matker()
        assert len(kernel) == 1 and len(kernel[0]) == len(basis) and not any(matrix*kernel[i] for i in range(len(kernel)))
        print(f"degree {degree}, order {k}: {len(basis)} columns and a one-dimensional exact kernel")
    # The two-step graph construction can also be checked over a split 2-division field.
    field = lambda value: pari.Mod(value, 101)
    first = (field(-1), field(0))
    second, h = velu_graph(first, [field(0), field(0)], [field(1), field(0)], 2)
    graph, _basis = section_matrix(first, second, h, 3, 3, 2, 2)
    assert len(graph.matker()) == 1
    print("The normalized dual Vélu graph gives the same section dimension.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--parameter", type=int)
    parser.add_argument("--isolate-torsion", action="store_true")
    args = parser.parse_args()
    if args.parameter is not None:
        print(json.dumps(search(args.parameter), indent=2))
        return
    data = json.loads((HERE / "family-equations.json").read_text())
    for record in data:
        params = {k: record[k] for k in ("a", "b", "n", "m", "k")}
        assert record == certificate(record["t"], record["curve"], params)
        frequencies = gram(record["humbert_gram"]).qfrep(81)
        assert not any(frequencies[d*d-1] for d in range(1, 9)) and frequencies[80] == 1
    print(f"{len(data)} specialized polarized products and all four degree-nine minima verified.")
    if args.isolate_torsion:
        labels = [{"t": record["t"], "curve": record["curve"], **isolate_torsion(record)} for record in data]
        archive = json.loads((HERE / "family-torsion-roots.json").read_text())
        assert labels == archive, "The labelled division roots differ from the receipt."
        print(f"{2*len(labels)} exact root boxes verified against all division points.")
    check_control()


if __name__ == "__main__":
    main()
