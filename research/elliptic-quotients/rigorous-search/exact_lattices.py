"""Recover exact homology actions of known algebraic maps and saturate their span."""
from fractions import Fraction
from itertools import product
from math import lcm

import cypari2
import sympy as sp
from flint import acb, acb_mat, arb
from sympy.matrices.normalforms import hermite_normal_form

from enumeration import short_vectors
from search import PAIRS, integer, normalize, real_frame


pari = cypari2.Pari()
I = sp.eye(4)


def rational(value):
    value = Fraction(value)
    return arb(value.numerator) / value.denominator


def gaussian(value):
    return acb(*(rational(t) for t in value))


def exact_entries(matrix):
    return sp.Matrix(matrix.nrows(), matrix.ncols(), lambda i, j: integer(matrix[i, j]))


def homology_action(data, tangent):
    P = acb_mat(data["Pi"])
    image = tangent * P
    return exact_entries(real_frame(data["Pi"]).inv() * real_frame(
        [[image[i, j] for j in range(4)] for i in range(2)]))


def norm_endomorphism(data, G, a, b, degree):
    """Use the polarization to recover f^* f_* of an exact degree-d map."""
    F = real_frame([[a * p + b * q for p, q in zip(*data["Pi"])]])
    dual = G.inv() * F.transpose()
    return exact_entries(degree * dual * (F * dual).inv() * F)


def columns(matrix):
    return [matrix[:, j] for j in range(matrix.cols)]


def flat(matrix):
    return matrix.reshape(16, 1)


def integer_kernel(matrix):
    """Full integral kernel via PARI's HNF-based exact matkerint."""
    A = pari.matrix(matrix.rows, matrix.cols, list(map(int, matrix)))
    K = pari.matkerint(A)
    B = sp.Matrix(matrix.cols, len(K), lambda i, j: int(K[i, j]))
    if matrix * B != sp.zeros(matrix.rows, B.cols) or B.cols != matrix.cols - matrix.rank():
        raise ArithmeticError("Invalid integral kernel.")
    return B


def generated_order(generators, rank):
    """The products are honest endomorphisms, so reaching rank fixes the rational span."""
    basis = hermite_normal_form(sp.Matrix.hstack(*(flat(R) for R in generators)))
    while basis.cols < rank:
        mats = [v.reshape(4, 4) for v in columns(basis)]
        larger = hermite_normal_form(sp.Matrix.hstack(basis, *(flat(a*b) for a, b in product(mats, repeat=2))))
        if larger.cols == basis.cols:
            raise ArithmeticError("The algebraic anchors did not generate the required rank.")
        basis = larger
    if basis.cols != rank:
        raise ArithmeticError("The algebraic anchors exceeded the proven rank bound.")
    return basis


def saturate(basis):
    return hermite_normal_form(integer_kernel(integer_kernel(basis.T).T))


def intersection(A, B):
    return (sp.trace(A) * sp.trace(B)) / 4 - sp.trace(A * B) / 2


def lattices(E, span):
    end = saturate(span)
    matrices = [v.reshape(4, 4) for v in columns(end)]
    pairing = sp.Matrix(E)
    adjoints = [pairing.inv() * R.T * pairing for R in matrices]
    change = sp.Matrix.hstack(*(flat(R - S) for R, S in zip(adjoints, matrices)))
    ns = [v.reshape(4, 4) for v in columns(end * integer_kernel(change))]
    gram = sp.Matrix([[intersection(A, B) for B in ns] for A in ns])
    theta, _ = sp.Matrix.hstack(*(flat(A) for A in ns)).gauss_jordan_solve(flat(I))
    trace_gram = [[int(sp.trace(A * B)) for B in matrices] for A in adjoints]
    if any(not x.is_Integer for x in [*gram, *theta]):
        raise ArithmeticError("The Néron–Severi pairing is not integral.")
    return end, ns, gram, theta, trace_gram


def plane_set(E, ns, classes):
    planes = {}
    for degree, coeffs in classes:
        R = sp.zeros(4)
        for x, A in zip(coeffs, ns):
            R += x * A
        # The class of an elliptic subvariety has norm zero on its tangent.
        kernel = integer_kernel(R)
        if kernel.cols != 2:
            raise ArithmeticError("An isotropic class did not give a rank-two plane.")
        pi = normalize(tuple(int(kernel[i, 0] * kernel[j, 1] - kernel[j, 0] * kernel[i, 1])
                             for i, j in PAIRS))
        if abs(sum(E[i][j] * c for (i, j), c in zip(PAIRS, pi))) != degree or pi in planes:
            raise ArithmeticError("The elliptic plane disagrees with its primitive class.")
        planes[pi] = degree
    return planes


def rosati_count(gram, bound):
    """Counts nonzero endomorphisms up to sign, including imprimitive ones."""
    return sum(1 for _ in short_vectors(gram, bound))


def check_qm(end, gens, gram):
    A, B = gens
    words = (I, A, B, A * B, B * A, A * A, B * B)
    if hermite_normal_form(sp.Matrix.hstack(*(flat(R) for R in words))) != end:
        raise ArithmeticError("The exact QM generators did not generate the maximal order.")
    matrices = [v.reshape(4, 4) for v in columns(end)]
    reduced_trace = sp.Matrix([[sp.trace(R * S) / 2 for S in matrices] for R in matrices])
    if reduced_trace.det() != -36:
        raise ArithmeticError("The reduced trace discriminant differs from -36.")
    for modulus in (3, 9, 27, 81):
        if not any(
            any(x % 3 for x in v)
            and sum(v[i] * gram[i, j] * v[j] for i in range(3) for j in range(3)) % modulus == 0
            for v in product(range(modulus), repeat=3)
        ):
            return modulus
    raise ArithmeticError("The QM form failed the local anisotropy check.")


def check_end_shortlist(G, J, end, trace_gram, bound, penalty=2 ** 20):
    """Compare a bounded ambient search to the exact rank-and-saturation result."""
    from flint import arb_mat
    from search import dyadic_bounds, endomorphism_form, integral_lower_bound

    Q = endomorphism_form(G, J, penalty)
    scale = 2 ** 24
    A = sp.Matrix(integral_lower_bound(Q, scale))
    proposal = pari.qflllgram(pari.matrix(16, 16, list(map(int, A))))
    U = sp.Matrix(16, 16, lambda i, j: int(proposal[i, j]))
    if abs(U.det()) != 1:
        raise ArithmeticError("The enumeration basis change is not unimodular.")
    reduced = (U.T * A * U).tolist()
    exact_change = U.inv() * end
    equations = integer_kernel(exact_change.T).T
    conditions = [[(i, int(value)) for i, value in enumerate(row) if value]
                  for row in equations.tolist()]
    pivot = list(exact_change.T.rref()[1])
    coordinates = exact_change[pivot, :].inv()
    gram = coordinates.T * sp.Matrix(trace_gram) * coordinates
    denominator = lcm(*(int(value.q) for value in gram))
    terms = [(i, j, int(denominator * gram[a, b] * (2 if a != b else 1)))
             for a, i in enumerate(pivot) for b, j in enumerate(pivot[:a + 1]) if gram[a, b]]


    resolution = 2 ** 50
    rulings = []
    for v in columns(U):
        R = arb_mat([list(map(int, row)) for row in v.reshape(4, 4).tolist()])
        C = J * R - R * J
        rulings.append([dyadic_bounds(C[i, j], resolution) for i in range(4) for j in range(4)])

    accepted, rejected = 0, 0
    for u in short_vectors(reduced, scale * bound):
        if any(sum(u[i] * value for i, value in condition) for condition in conditions):
            parts = [(u[i], row) for i, row in enumerate(rulings) if u[i]]
            for j in range(16):
                lo = sum(x * row[j][0 if x > 0 else 1] for x, row in parts)
                hi = sum(x * row[j][1 if x > 0 else 0] for x, row in parts)
                if lo > 0 or hi < 0:
                    break
            else:
                raise ArithmeticError("The ambient shortlist encountered an unresolved commutator.")
            rejected += 1
        elif sum(value * u[i] * u[j] for i, j, value in terms) <= denominator * bound:
            accepted += 1
    return accepted, rejected


def check_stored_ns(anchors, maps, ns, stored):
    """Realize the torsion-certified generators in the independently saturated lattice."""
    classes = [I, *(d * I - R for (_, _, d), R in zip(maps, anchors[1:]))]
    G = sp.Matrix([[intersection(A, B) for B in classes] for A in classes])
    if G != sp.Matrix(stored["gram"]):
        raise ArithmeticError("The certified map classes give a different Gram matrix.")
    generators = sp.Matrix.hstack(*(flat(A) for A in classes))
    stored_basis = generators * sp.Matrix(stored["NS_basis_in_generators"]).applyfunc(sp.Rational).T
    ns_basis = sp.Matrix.hstack(*(flat(R) for R in ns))
    transform, _ = ns_basis.gauss_jordan_solve(stored_basis)
    if any(not x.is_Integer for x in transform) or abs(transform.det()) != 1:
        raise ArithmeticError("The torsion-certified basis differs from the exact lattice.")
