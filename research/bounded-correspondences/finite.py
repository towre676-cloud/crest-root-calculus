"""Exact rational points of a zero-dimensional polynomial incidence over Q.

The trace construction follows rational univariate representation (Rouillier).
The local stage uses a derived numerator/denominator bound and exact replay.
"""
from collections import deque
from math import prod
import os
import shutil
import subprocess

import sympy as sp


T = sp.Symbol("T")


def coefficients(poly):
    """Serialize a univariate rational polynomial in increasing powers."""
    return [str(poly.nth(i)) for i in range(poly.degree() + 1)] if poly else ["0"]


def polynomial(values, variable=T):
    return sp.Poly(sum(sp.Rational(a) * variable**i for i, a in enumerate(values)),
                   variable, domain=sp.QQ)


def staircase(basis):
    count = len(basis.gens)
    leading = [p.LM(order=basis.order).exponents for p in basis.polys]
    seen = set()
    pending = deque([(0,) * count])
    while pending:
        monomial = pending.popleft()
        if monomial in seen or any(all(a >= b for a, b in zip(monomial, lm))
                                   for lm in leading):
            continue
        seen.add(monomial)
        for j in range(count):
            next_monomial = list(monomial)
            next_monomial[j] += 1
            pending.append(tuple(next_monomial))
    return sorted(seen, key=lambda exponents: (sum(exponents), exponents))


def multiplication_algebra(equations, variables, method="f5b", supplied_basis=None):
    proof = None
    executable = os.environ.get("CREST_SINGULAR") or shutil.which("Singular")
    if supplied_basis is not None:
        basis = sp.groebner(supplied_basis, *variables, order="grevlex", domain=sp.QQ, method=method)
    elif executable:
        algorithm = "slimgb" if method == "f5b" else "std"
        names = ",".join(map(str, variables))
        rendered = []
        for equation in equations:
            terms = []
            for exponents, coefficient in sp.Poly(equation, *variables, domain=sp.QQ).terms():
                term = "(" + str(coefficient) + ")"
                for variable, exponent in zip(variables, exponents):
                    if exponent:
                        term += "*" + str(variable) + "^" + str(exponent)
                terms.append(term)
            rendered.append("+".join(terms))
        expressions = ",".join(rendered)
        script = (f"ring R=0,({names}),dp;\nshort=0; option(redSB);\n"
                  f"ideal I={expressions};\nideal G=interred({algorithm}(I));\n"
                  'print("BEGIN_BASIS"); int i;\n'
                  'for(i=1;i<=size(G);i++){string(G[i]);}\n'
                  'print("END_BASIS");\n')
        if method == "f5b":
            script += ('matrix L=lift(I,G); print("BEGIN_COFACTORS"); int j;\n'
                       f'for(i=1;i<=size(G);i++){{for(j=1;j<={len(equations)};j++){{string(L[j,i]);}}}}\n'
                       'print("END_COFACTORS");\n')
        script += "quit;\n"
        process = subprocess.run([executable, "-q"], input=script, text=True,
                                 capture_output=True, check=True)
        if "?" in process.stdout or "BEGIN_BASIS" not in process.stdout or "END_BASIS" not in process.stdout:
            raise ArithmeticError("Singular did not return a complete exact basis: " + process.stdout[-600:])
        lines = process.stdout.split("BEGIN_BASIS", 1)[1].split("END_BASIS", 1)[0].splitlines()
        context = dict(zip(map(str, variables), variables))
        polynomials = [sp.sympify(line.replace("^", "**"), locals=context)
                       for line in lines if line.strip()]
        basis = sp.groebner(polynomials, *variables, order="grevlex", domain=sp.QQ,
                           method=method)
        if method == "f5b":
            if "BEGIN_COFACTORS" not in process.stdout or "END_COFACTORS" not in process.stdout:
                raise ArithmeticError("The original-ideal cofactors were not returned.")
            lines = process.stdout.split("BEGIN_COFACTORS", 1)[1].split("END_COFACTORS", 1)[0].splitlines()
            values = [sp.sympify(line.replace("^", "**"), locals=context)
                      for line in lines if line.strip()]
            if len(values) != len(equations) * len(polynomials):
                raise ArithmeticError("The original-ideal cofactor matrix is incomplete.")
            original = [sp.Poly(q, *variables, domain=sp.QQ) for q in polynomials]
            proof = []
            for q in basis.polys:
                index = next(i for i, native in enumerate(original)
                             if native.mul_ground(1 / native.LC(order="grevlex")) == q)
                start = index * len(equations)
                proof.append([str(value / original[index].LC(order="grevlex"))
                              for value in values[start:start + len(equations)]])
    else:
        basis = sp.groebner(equations, *variables, order="grevlex", domain=sp.QQ,
                           method=method)
    if list(basis) == [sp.Integer(1)]:
        return basis, [], [], proof
    if not basis.is_zero_dimensional:
        raise ValueError("The supplied incidence has positive dimension.")
    monomials = staircase(basis)
    expressions = [prod(x**e for x, e in zip(variables, monomial))
                   for monomial in monomials]
    indices = {monomial: j for j, monomial in enumerate(monomials)}
    matrices = []
    for x in variables:
        matrix = sp.zeros(len(monomials))
        for j, expression in enumerate(expressions):
            remainder = sp.Poly(basis.reduce(x * expression)[1], *variables,
                                domain=sp.QQ)
            for exponents, value in remainder.terms():
                matrix[indices[exponents], j] = value
        matrices.append(matrix)
    return basis, monomials, matrices, proof


def evaluate_matrix(poly, matrix):
    result = sp.zeros(matrix.rows)
    identity = sp.eye(matrix.rows)
    for value in poly.all_coeffs():
        result = result * matrix + value * identity
    return result


def separating_presentation(matrices, monomials):
    dimension = matrices[0].rows
    identity = sp.eye(dimension)
    basis_matrices = [prod((matrix**e for matrix, e in zip(matrices, exponents)),
                           start=identity) for exponents in monomials]
    trace_pairing = sp.Matrix([[sp.trace(a * b) for b in basis_matrices]
                              for a in basis_matrices])
    support_size = trace_pairing.rank()
    maximum_trials = (len(matrices) - 1) * support_size * (support_size - 1) // 2 + 1
    for k in range(maximum_trials):
        linear_coefficients = [k**i for i in range(len(matrices))]
        linear_matrix = sum((a * matrix for a, matrix in zip(linear_coefficients, matrices)),
                            sp.zeros(dimension))
        characteristic = sp.Poly(linear_matrix.charpoly(T).as_expr(), T, domain=sp.QQ)
        squarefree = characteristic.sqf_part().monic()
        if squarefree.degree() == support_size:
            break
    else:
        raise ArithmeticError("The finite separating-form bound failed.")
    powers = [identity]
    for _ in range(2 * support_size - 2):
        powers.append(powers[-1] * linear_matrix)
    moments = [sp.trace(matrix) for matrix in powers]
    hankel = sp.Matrix(support_size, support_size, lambda i, j: moments[i + j])
    coordinate_polynomials = []
    for matrix in matrices:
        right = sp.Matrix([sp.trace(matrix * powers[j]) for j in range(support_size)])
        values = hankel.inv() * right
        coordinate_polynomials.append(sp.Poly(sum(values[j] * T**j
                                                  for j in range(support_size)), T))
    _, primitive = squarefree.clear_denoms(convert=True)
    primitive = primitive.primitive()[1]
    if primitive.LC() < 0:
        primitive = -primitive
    return linear_coefficients, linear_matrix, primitive, coordinate_polynomials, support_size


def rational_reconstruction(residue, modulus, numerator_bound, denominator_bound):
    r0, r1 = modulus, residue % modulus
    s0, s1 = 0, 1
    while r1 > numerator_bound:
        quotient, remainder = divmod(r0, r1)
        r0, r1 = r1, remainder
        s0, s1 = s1, s0 - quotient * s1
    if s1 == 0:
        return None
    numerator, denominator = (r1, s1) if s1 > 0 else (-r1, -s1)
    if denominator > denominator_bound or sp.gcd(numerator, denominator) != 1:
        return None
    if (numerator - residue * denominator) % modulus:
        return None
    return sp.Rational(numerator, denominator)


def local_rational_roots(poly, coordinate_polynomials, prime=None):
    if (poly.degree() < 1 or any(value.q != 1 for value in poly.all_coeffs())
            or poly.LC() <= 0 or poly.primitive()[0] != 1):
        raise ValueError("Local recovery requires a primitive integral polynomial with positive leading coefficient.")
    leading = abs(int(poly.LC()))
    numerator_bound = leading + max(map(abs, map(int, poly.all_coeffs()[1:])))
    denominator_bound = leading
    excluded = leading * abs(int(sp.discriminant(poly)))
    for coordinate in coordinate_polynomials:
        excluded *= sp.ilcm(*[int(a.q) for a in coordinate.all_coeffs()], 1)
    if prime is None:
        prime = 2
        while excluded % prime == 0:
            prime = int(sp.nextprime(prime))
    if not sp.isprime(prime) or excluded % prime == 0:
        raise ValueError("The supplied prime divides a leading coefficient, discriminant, or denominator.")
    roots = [r for r in range(prime) if int(poly.eval(r)) % prime == 0]
    modulus = prime
    lifted = list(roots)
    derivative = poly.diff()
    while modulus <= 2 * numerator_bound * denominator_bound:
        modulus *= modulus
        lifted = [(r - int(poly.eval(r)) * pow(int(derivative.eval(r)), -1, modulus)) % modulus
                  for r in lifted]
    records = []
    rational_roots = []
    for initial, residue in zip(roots, lifted):
        candidate = rational_reconstruction(residue, modulus, numerator_bound, denominator_bound)
        valid = candidate is not None and poly.eval(candidate) == 0
        records.append({"initial_root": initial, "lifted_root": str(residue),
                        "candidate": str(candidate) if candidate is not None else None,
                        "exact_root": bool(valid)})
        if valid:
            rational_roots.append(candidate)
    return rational_roots, {"prime": prime, "modulus": str(modulus),
                            "numerator_bound": str(numerator_bound),
                            "denominator_bound": str(denominator_bound),
                            "branches": records}


def rational_points(equations, variables, prime=None):
    basis, monomials, matrices, proof = multiplication_algebra(equations, variables)
    result = {"variables": list(map(str, variables)),
              "equations": list(map(str, equations)),
              "groebner_basis": list(map(str, basis)),
              "monomials": monomials,
              "algebra_dimension": len(monomials)}
    if proof is not None:
        result["original_cofactors"] = proof
    if not monomials:
        return {**result, "geometric_points": 0, "rational_points": []}
    linear, linear_matrix, poly, coordinates, support_size = separating_presentation(matrices, monomials)
    roots, lifting = local_rational_roots(poly, coordinates, prime)
    points = [[str(coordinate.eval(root)) for coordinate in coordinates] for root in roots]
    for point in points:
        substitution = dict(zip(variables, map(sp.Rational, point)))
        if any(sp.expand(equation).subs(substitution) != 0 for equation in equations):
            raise ArithmeticError("A reconstructed point failed its original equation.")
    return {**result, "geometric_points": support_size,
            "multiplication_matrices": [[[str(v) for v in row] for row in matrix.tolist()]
                                         for matrix in matrices],
            "linear_form": linear, "minimal_polynomial": coefficients(poly),
            "coordinates": [coefficients(coordinate) for coordinate in coordinates],
            "lifting": lifting, "rational_points": points}


def verify(result):
    """Prove ideal equality by cofactors or independent elimination, then replay."""
    variables = sp.symbols(result["variables"])
    context = dict(zip(result["variables"], variables))
    equations = [sp.sympify(s, locals=context) for s in result["equations"]]
    supplied_basis = None
    if "original_cofactors" in result:
        supplied_basis = [sp.sympify(q, locals=context) for q in result["groebner_basis"]]
        proof = result["original_cofactors"]
        if len(proof) != len(supplied_basis) or any(len(row) != len(equations) for row in proof):
            raise ArithmeticError("The original-ideal cofactor matrix is incomplete.")
        original = [sp.Poly(equation, *variables, domain=sp.QQ) for equation in equations]
        for target, row in zip(supplied_basis, proof):
            combination = sum((sp.Poly(sp.sympify(value, locals=context), *variables, domain=sp.QQ) * equation
                               for value, equation in zip(row, original)),
                              sp.Poly(0, *variables, domain=sp.QQ))
            if combination != sp.Poly(target, *variables, domain=sp.QQ):
                raise ArithmeticError("An original-ideal cofactor identity failed.")
    basis, monomials, matrices, _ = multiplication_algebra(equations, variables, method="buchberger",
                                                          supplied_basis=supplied_basis)
    if any(basis.reduce(equation)[1] != 0 for equation in equations):
        raise ArithmeticError("An original equation is absent from the presented ideal.")
    if list(map(str, basis)) != result["groebner_basis"]:
        raise ArithmeticError("The supplied basis differs from the original incidence ideal.")
    if result["algebra_dimension"] != len(monomials):
        raise ArithmeticError("Incorrect quotient-algebra dimension.")
    if monomials != list(map(tuple, result["monomials"])):
        raise ArithmeticError("The standard-monomial basis differs.")
    if not matrices:
        if result["geometric_points"] != 0 or result["rational_points"] != []:
            raise ArithmeticError("Incorrect empty incidence.")
        return
    if [[[str(v) for v in row] for row in matrix.tolist()] for matrix in matrices] != result["multiplication_matrices"]:
        raise ArithmeticError("Incorrect multiplication matrices.")
    dimension = len(monomials)
    if len(result["linear_form"]) != len(matrices) or len(result["coordinates"]) != len(matrices):
        raise ArithmeticError("The presentation omits or adds a coordinate.")
    linear_matrix = sum((a * matrix for a, matrix in zip(result["linear_form"], matrices)),
                        sp.zeros(dimension))
    poly = polynomial(result["minimal_polynomial"])
    coordinates = [polynomial(values) for values in result["coordinates"]]
    expected = sp.Poly(linear_matrix.charpoly(T).as_expr(), T).sqf_part()
    expected = expected.clear_denoms(convert=True)[1].primitive()[1]
    if expected.LC() < 0:
        expected = -expected
    if poly != expected.set_domain(sp.QQ) or poly.degree() != result["geometric_points"]:
        raise ArithmeticError("Incorrect reduced separating polynomial.")
    for matrix, coordinate in zip(matrices, coordinates):
        if (matrix - evaluate_matrix(coordinate, linear_matrix))**dimension != sp.zeros(dimension):
            raise ArithmeticError("A coordinate fails on the geometric support.")
    if sp.rem(sum(a * q.as_expr() for a, q in zip(result["linear_form"], coordinates)) - T,
              poly.as_expr(), T) != 0:
        raise ArithmeticError("The presentation does not recover its separating coordinate.")
    substitution = dict(zip(variables, [q.as_expr() for q in coordinates]))
    for equation in equations:
        if sp.rem(equation.subs(substitution, simultaneous=True), poly.as_expr(), T) != 0:
            raise ArithmeticError("A presentation coordinate fails an original equation.")
    roots, lifting = local_rational_roots(poly, coordinates, result["lifting"]["prime"])
    points = [[str(coordinate.eval(root)) for coordinate in coordinates] for root in roots]
    if lifting != result["lifting"] or points != result["rational_points"]:
        raise ArithmeticError("The local recovery certificate differs.")
