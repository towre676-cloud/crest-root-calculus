"""Certify a prime-power collision and an unshared square-null Hecke eigenline."""
from fractions import Fraction
from itertools import product
import json
from math import gcd, isqrt, lcm

from neighbours import Closure, gram_matrix, pari
from observe import degree_counts, form_gram, residue_field, scalar_orbits, square_shells


def genera(D, ell, p):
    field = residue_field(D, ell)
    splits = {}
    for sign in (1, -1):
        centers = [a for a in scalar_orbits(field) if field.legendre(1 + field.norm(a)) == sign]
        closure = Closure([form_gram(field, a) for a in centers], (p,))
        splits[sign] = centers, closure
    return splits


def prime_collision():
    D, ell, p = -3, 5, 17
    field = residue_field(D, ell)
    examples = []
    for sign, (centers, closure) in genera(D, ell, p).items():
        if len(centers) != 2:
            raise ArithmeticError("Expected two gluing forms.")
        rows = [[r.count(j) for j in range(len(closure.forms))] for r in closure.rows[p]]
        norm_one = [2 * int(q.qfrep(1)[0]) for q in closure.forms]
        if norm_one[0] != norm_one[1] or rows[0][2:] != rows[1][2:] or sum(rows[0][:2]) != sum(rows[1][:2]):
            raise ArithmeticError("Neighbour matrix does not preserve equality on the pair.")
        shells = square_shells(field, 2)
        ordinary_witness = [degree_counts(field, a, [shells[1]])[0] for a in centers]
        if len(set(ordinary_witness)) != 2:
            raise ArithmeticError("Expected a degree-two primitive distinction.")
        examples.append({"sign": sign, "centers": centers, "neighbour_counts": rows,
                         "norm_one": norm_one, "degree_two_counts": ordinary_witness})
    return {"D": D, "ell": ell, "p": p, "examples": examples}


def forms_of_determinant(det):
    """Enumerate a superset of size-reduced integral LLL Gram representatives."""
    a = 1
    while a**3 <= 8 * det:
        for b in range(1, isqrt(2 * det // a) + a//4 + 2):
            for u in range(-(a//2), a//2 + 1):
                h = a*b - u*u
                if h <= 0:
                    continue
                for v in range(-(a//2), a//2 + 1):
                    lower = -((h - 2*u*v)//(2*a))
                    upper = (h + 2*u*v)//(2*a)
                    for w in range(lower, upper + 1):
                        c, remainder = divmod(det + a*w*w + b*v*v - 2*u*v*w, h)
                        if remainder or any(s % 4 not in (0, 1) for s in (a,b,c)):
                            continue
                        if all((a*x*x + b*y*y + c*z*z + 2*(u*x*y + v*x*z + w*y*z)) % 4 < 2
                               for x, y, z in product(range(4), repeat=3)):
                            yield [[a,u,v],[u,b,w],[v,w,c]]
        a += 1


def primary_discriminant_form(gram, p):
    """Give the coordinate orders and Gram pairing of the p-primary part."""
    even_gram = 2 * gram
    U, _, diagonal = even_gram.matsnf(1)
    dual_basis = even_gram**(-1) * U**(-1)
    generators, orders = [], []
    for i in range(3):
        order = int(diagonal[i,i])
        power = 1
        while order % (power*p) == 0:
            power *= p
        if power > 1:
            orders.append(power)
            generators.append(order//power * dual_basis[i])
    pairing = [[Fraction(str((v.mattranspose() * even_gram * w)[0])) for w in generators]
               for v in generators]
    return orders, pairing


def isometric_primary(source, target):
    orders, A = source
    target_orders, B = target
    if orders != target_orders:
        return False

    def pair(v, w):
        return sum(x * B[i][j] * y for i, x in enumerate(v) for j, y in enumerate(w))

    elements = list(product(*(range(n) for n in orders)))
    candidates = [
        [v for v in elements
         if lcm(*(n//gcd(x,n) for x,n in zip(v,orders))) == orders[i] and (pair(v,v)-A[i][i])%2 == 0]
        for i in range(len(orders))
    ]
    for images in product(*candidates):
        if any((pair(images[i],images[j])-A[i][j])%1 for i in range(len(orders)) for j in range(i)):
            continue
        image = {tuple(sum(c*v[j] for c,v in zip(coords,images)) % n for j,n in enumerate(orders))
                 for coords in elements}
        if len(image) == len(elements):
            return True
    return False


def certify_complete_genera(closures, det, bad_primes):
    """Exhaust the determinant, using the primary discriminant-form genus criterion."""
    targets = [[primary_discriminant_form(c.forms[0], p) for p in bad_primes] for c in closures]
    matched = [set(),set()]
    candidates = 0
    for gram in forms_of_determinant(det):
        candidates += 1
        q = gram_matrix(gram)
        data = [primary_discriminant_form(q, p) for p in bad_primes]
        for j, target in enumerate(targets):
            if not all(isometric_primary(s,t) for s,t in zip(data,target)):
                continue
            i = next((i for i,m in enumerate(closures[j].forms) if m.qfisom(q)),None)
            if i is None:
                raise ArithmeticError("The neighbour closure omits a genus class.")
            matched[j].add(i)
    if any(len(m) != len(c.forms) for m,c in zip(matched,closures)):
        raise ArithmeticError("The exhaustive enumeration did not recover every neighbour class.")
    return {"size_reduced_gram_candidates_after_mod_4_filter": candidates,
            "genus_class_counts": list(map(len,matched))}


def shimura_coefficients(theta_coefficients, D, t, level):
    bound = len(theta_coefficients) - 1
    lift = [0] * (bound + 1)
    for d in range(1, bound + 1):
        if gcd(d, level) == 1:
            character = int(pari.kronecker(D*t, d))
            for m in range(1, bound//d + 1):
                lift[d*m] += character * theta_coefficients[m]
    return lift


def eigenline():
    """The (-7,5) genus minus has a nonzero theta eigenspace with T(3^2)=1 only there."""
    D, ell, p = -7, 5, 3
    split = genera(D, ell, p)
    closures = [split[s][1] for s in (1,-1)]
    completeness = certify_complete_genera(closures, 700, (2, 5, 7))
    matrices = []
    for c in closures:
        n = len(c.forms)
        matrices.append(pari.matrix(n,n,[r.count(j) for r in c.rows[p] for j in range(n)]).mattranspose())
    eigenspaces = [(transposed - pari.matid(len(transposed))).matker() for transposed in matrices]
    if list(map(len,eigenspaces)) != [0,1]:
        raise ArithmeticError("The exact eigenspace multiplicities changed.")
    vector = list(map(int,eigenspaces[1][0]))
    if sum(vector):
        raise ArithmeticError("The signed combination is not cuspidal.")
    # M=700, index 1440, so the weight-two Sturm bound at M is 240.
    cutoff = 240
    coefficients = [q.qfrep(5*cutoff**2) for q in closures[1].forms]
    def theta(n):
        return sum(2*c*int(r[n-1]) for c,r in zip(vector,coefficients))
    if theta(5) == 0 or any(theta(m*m) for m in range(1,cutoff+1)):
        raise ArithmeticError("The unshared eigenline must have nonzero theta and zero square counts.")
    ainvs = [0, 1, 1, -131, -650]
    curve = pari.ellinit(ainvs)
    if curve.ellglobalred()[0] != 35:
        raise ArithmeticError("Wrong conductor for the comparison form.")
    a = [0] + list(map(int, curve.ellan(cutoff)))
    lift = shimura_coefficients([0] + [theta(5*m*m) for m in range(1,cutoff+1)], D, 5, 700)
    expected = [theta(5) * (a[m] + (a[m//2] if m % 2 == 0 else 0)) for m in range(cutoff+1)]
    if lift != expected:
        raise ArithmeticError("The level-35 oldform identity failed.")
    twist = pari.ellinit(curve.elltwist(D))
    symbols, symbol = curve.msfromell(1)
    central_ratio = str(pari.mseval(symbols, symbol, [pari('oo'), 0]))
    if twist.ellrootno() != -1 or central_ratio != "1":
        raise ArithmeticError("Unexpected central-value data for the level-35 form.")
    return {
        "D": D, "ell": ell, "p": p, "lambda": 1, "multiplicities_by_genus": [0,1],
        "genus_completeness": completeness,
        "neighbour_characteristic_polynomials": [str(m.charpoly()) for m in matrices],
        "genus_minus_forms": [[[int(q[i,j]) for j in range(3)] for i in range(3)] for q in closures[1].forms],
        "genus_minus_vector": vector, "r_5": theta(5), "square_sturm_bound": cutoff,
        "weight_two_curve_conductor": 35, "weight_two_curve_ainvs": ainvs,
        "shimura_t5_oldform_weights_at_1_and_2": [theta(5)] * 2,
        "untwisted_L_over_real_period": central_ratio,
        "twist_conductor": int(twist.ellglobalred()[0]), "twist_root_number": int(twist.ellrootno()),
    }


def shared_packet():
    """Resolve all three copies of the level-11 eigensystem at (-11,7)."""
    D, ell, p = -11, 7, 3
    split = genera(D, ell, p)
    closures = [split[s][1] for s in (1, -1)]
    completeness = certify_complete_genera(closures, 2156, (2, 7, 11))
    # M=2156, index 4032, so the weight-two Sturm bound at M is 672.
    level, bound = 2156, 672
    ainvs = [0, -1, 1, -10, -20]
    curve = pari.ellinit(ainvs)
    if curve.ellglobalred()[0] != 11:
        raise ArithmeticError("Wrong conductor for the comparison form.")
    a = [0] + list(map(int, curve.ellan(bound)))
    oldform = [a[m] + (a[m//2] if m % 2 == 0 else 0) for m in range(bound+1)]
    weight_rows, dimensions = [], []
    for c in closures:
        n = len(c.forms)
        T = pari.matrix(n,n,[r.count(j) for r in c.rows[p] for j in range(n)]).mattranspose()
        vectors = (T + pari.matid(n)).matker()
        dimensions.append(len(vectors))
        counts = []
        for q in c.forms:
            r = q.qfrep(bound**2)
            counts.append([int(r[m*m-1]) for m in range(1,bound+1)])
        for raw in vectors:
            v = list(map(int, raw))
            if sum(v):
                raise ArithmeticError("Expected a cuspidal theta combination.")
            squares = [0] + [sum(2*b*x for b,x in zip(v,row)) for row in zip(*counts)]
            lift = shimura_coefficients(squares, D, 1, level)
            weights = [lift[1], lift[7] - lift[1]*oldform[7], 0]
            weights[2] = lift[49] - weights[0]*oldform[49] - weights[1]*oldform[7]
            expected = [sum(w*oldform[m//d] for w,d in zip(weights,(1,7,49)) if m % d == 0)
                        for m in range(bound+1)]
            if lift != expected:
                raise ArithmeticError("The level-11 oldform identity failed.")
            weight_rows.append(weights)
    determinant = int(gram_matrix(weight_rows).matdet())
    if dimensions != [1,2] or determinant == 0:
        raise ArithmeticError("Three independent square lifts were expected.")
    return {"D": D, "ell": ell, "p": p, "lambda": -1, "multiplicities_by_genus": dimensions,
            "genus_completeness": completeness, "shimura_t1_sturm_bound": bound,
            "weight_two_curve_conductor": 11, "weight_two_curve_ainvs": ainvs,
            "oldform_weights_at_1_7_and_49": weight_rows, "oldform_weight_determinant": determinant}


if __name__ == "__main__":
    pari.allocatemem(32_000_000, 256_000_000, silent=True)
    print(json.dumps({"prime_power_collision": prime_collision(), "unshared_eigenline": eigenline(),
                      "shared_packet": shared_packet()}, indent=2))
