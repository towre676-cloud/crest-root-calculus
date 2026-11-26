"""Audit the imported inert-level census with exact lattice and degree checks."""
from itertools import combinations, product
import json
from math import gcd, isqrt
from pathlib import Path
import sys

from census import residue_field, scalar_orbits
from gluing import ns_and_q, pari, q_from_ns
from verify import degree_counts, form_gram, square_shells

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "rigorous-search"))
from enumeration import short_vectors  # noqa: E402


def pari_matrix(rows):
    return pari.matrix(len(rows), len(rows[0]), [x for row in rows for x in row])


def check_isometry(first, second):
    A, B = map(pari_matrix, (first, second))
    U = pari.qfisom(A, B)
    if U == 0 or abs(pari.matdet(U)) != 1:
        raise ArithmeticError("The two integral Humbert forms are not isometric.")
    if pari.mattranspose(U)*A*U != B and pari.mattranspose(U)*B*U != A:
        raise ArithmeticError("The returned isometry fails its exact Gram identity.")


def scalar_parts(field, discriminant, psi):
    """Read the scalar and conjugate-linear parts in the basis used by census.py."""
    p = field.p
    shift = (discriminant-field.trace_omega)//2

    def to_standard(x, y):
        return field.reduce((x+shift*y, y))

    image_one = to_standard(*[row[0] for row in psi])
    image_tau = to_standard(*[2*row[1]-discriminant*row[0] for row in psi])
    tau = (-field.trace_omega, 2)
    plus = field.mul(image_tau, field.inv(tau))
    half = (pow(2, -1, p), 0)
    return (field.mul(half, field.add(image_one, plus)),
            field.mul(half, field.add(image_one, field.neg(plus))))


def exact_primitive_counts(gram, bound):
    out = [0]*bound
    for v in short_vectors(gram, bound*bound):
        if gcd(*v) != 1:
            continue
        length = sum(x*y*g for x, row in zip(v, gram) for y, g in zip(v, row))
        d = isqrt(length)
        if d*d == length:
            out[d-1] += 2
    return out


def audit_level(row):
    D, p = row["D"], row["ell"]
    field = residue_field(D, p)
    orbits = scalar_orbits(field)
    orbit_lookup = {a: rep for rep, candidates in orbits.items() for a in candidates}
    if row["products"] or row["nonscalar_gluings"] != (p+1)*len(orbit_lookup):
        raise ArithmeticError("The presentation total differs from the structural count.")
    seen, genera = set(), set()
    max_degree = max((w[0] for g in row["genera"] for w in g["witnesses"]), default=1)
    shells = square_shells(field, max_degree)
    for genus in row["genera"]:
        if genus["level"] != -4*D*p*p:
            raise ArithmeticError("The theta level disagrees with the local formula.")
        signs, spectra = set(), []
        for form in genus["forms"]:
            psi, q = form["psi"], form["q"]
            a, b = scalar_parts(field, D, psi)
            if (1+field.norm(a)-field.norm(b)) % p or b == (0, 0):
                raise ArithmeticError("An archived graph is not a nonscalar anti-isometry.")
            rep = orbit_lookup[a]
            if rep in seen or form["gluings"] != (p+1)*len(orbits[rep]):
                raise ArithmeticError("A Humbert orbit was repeated or had the wrong graph count.")
            seen.add(rep)
            signs.add(field.legendre(1+field.norm(a)))
            result = ns_and_q(D, p, psi)
            if result is None:
                raise ArithmeticError("A graph representative failed the polarization check.")
            check_isometry(q, result["q"])
            check_isometry(q, form_gram(field, a))
            counts = exact_primitive_counts(q, max_degree)
            if counts != degree_counts(field, a, shells):
                raise ArithmeticError("Direct lattice counts disagree with the scalar-part form.")
            spectra.append(counts)
        if len(signs) != 1 or genera.intersection(signs):
            raise ArithmeticError("The archived genus disagrees with the structural genus.")
        genera.update(signs)
        compared = []
        for left, right in combinations(spectra, 2):
            d = next((i+1 for i, (x,y) in enumerate(zip(left, right)) if x != y), None)
            if d is None:
                raise ArithmeticError("The archived pair lacks a separation witness.")
            # The uploaded qfrep receipts counted one representative per vector sign.
            compared.append([d, left[d-1]//2, right[d-1]//2])
        if compared != genus["witnesses"]:
            raise ArithmeticError("An archived first witness differs from exact enumeration.")
    if seen != set(orbits) or genera != {1, -1}:
        raise ArithmeticError("The archived census omits a form or a genus.")
    return len(seen), sum(len(g["witnesses"]) for g in row["genera"])


def audit_level_two():
    candidates = []
    for a, b, c, d in product(range(2), repeat=4):
        if (a*d-b*c) % 2 == 0:
            continue
        data = ns_and_q(-4, 2, [[a,b],[c,d]])
        if data is None:
            raise ArithmeticError("An admissible level-two Gaussian graph failed.")
        if not exact_primitive_counts(data["q"], 1)[0]:
            candidates.append(data["q"])
    if len(candidates) != 4:
        raise ArithmeticError("Unexpected number of nonproduct level-two presentations.")
    for suffix in (102, 120):
        receipt = json.loads((HERE.parent/"ns-certificates"/f"cert_cm_{suffix}.json").read_text())
        q = q_from_ns(receipt["NS_gram"], receipt["theta"])
        for candidate in candidates:
            check_isometry(q, candidate)
        counts = exact_primitive_counts(q, receipt["X"])
        for d, count in enumerate(counts, 1):
            if count != receipt["spectrum_through_X"].get(str(d), 0):
                raise ArithmeticError("The geometric degree count disagrees with its NS receipt.")


def audit_counterexample():
    field = residue_field(-163, 5)
    presentations = [[ [2,2], [1,3] ], [ [2,3], [4,3] ]]
    centers = [(0,0), (4,2)]   # zero and sqrt(-163) in the basis (1, (1+sqrt(-163))/2)
    spectra = []
    shells = square_shells(field, 42)
    for psi, expected in zip(presentations, centers):
        a, _ = scalar_parts(field, -163, psi)
        if a != expected or field.legendre(1+field.norm(a)) != 1:
            raise ArithmeticError("The counterexample has the wrong scalar part or genus.")
        data = ns_and_q(-163, 5, psi)
        check_isometry(data["q"], form_gram(field, a))
        counts = exact_primitive_counts(data["q"], 42)
        if counts != degree_counts(field, a, shells):
            raise ArithmeticError("The two exact counterexample calculations disagree.")
        spectra.append(counts)
    if spectra[0][:-1] != spectra[1][:-1] or [v[-1] for v in spectra] != [8, 0]:
        raise ArithmeticError("The degree-42 counterexample did not replay.")


def audit_nonmaximal_orders():
    """Independently reconstruct a conductor-dependent family over Q(i)."""
    l = 7
    for f in (1, 8, 15, 36, 71):
        D = -4*f*f
        # The lattice extractor uses w=(D+sqrt(D))/2 rather than f*i.
        shift = (D//2) % l
        basis = pari_matrix([[1, shift], [0, 1]])
        basis_inverse = pari_matrix([[1, -shift], [0, 1]])
        counts = []
        for linear in ([[1, 0], [0, -1]], [[1, 0], [2, -1]]):
            psi = basis_inverse*pari_matrix(linear)*basis
            psi = [[int(psi[i, j]) % l for j in range(2)] for i in range(2)]
            data = ns_and_q(D, l, psi)
            if data is None:
                raise ArithmeticError("A conductor-family graph failed the polarization check.")
            counts.append(exact_primitive_counts(data["q"], 2*f))
        if (counts[0][:-1] != counts[1][:-1]
                or counts[0][-1]-counts[1][-1] != 2
                or any(c[0] for c in counts)):
            raise ArithmeticError("The conductor-family separation did not replay.")


def main():
    receipt = json.loads((HERE/"imported-census.json").read_text())
    audit_level_two()
    print("Gaussian level 2 agrees with both exact NS certificates.", flush=True)
    total_forms = total_witnesses = 0
    for row in receipt["levels"]:
        forms, witnesses = audit_level(row)
        total_forms += forms
        total_witnesses += witnesses
        print(f'D={row["D"]}, ell={row["ell"]}: {forms} forms and {witnesses} first witnesses verified.',
              flush=True)
    audit_counterexample()
    audit_nonmaximal_orders()
    print(f"{total_forms} forms and {total_witnesses} witnesses verified across the imported census.")
    print("D=-163, ell=5: first difference at degree 42, with 8 versus 0 geometric classes.")
    print("Gaussian level 7: the first difference is 2f at conductors 1, 8, 15, 36 and 71.")


if __name__ == "__main__":
    main()
