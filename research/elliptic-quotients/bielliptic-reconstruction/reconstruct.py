"""Exact geometric elliptic quotient spectra from bielliptic genus-two equations."""

import argparse
from collections import Counter
from functools import reduce
from itertools import product
import json
from math import gcd, isqrt, pi, sqrt
import subprocess

import sympy as sp

x, u, z, a, b = sp.symbols("x u z a b")

# Class-number-one imaginary quadratic orders give the 13 rational CM j-values.
RATIONAL_CM = {0, 1728, -3375, 8000, -32768, 54000, 287496, -884736,
               -12288000, 16581375, -884736000, -147197952000, -262537412640768000}


def algebraic_involutions(basis, variables, C, count):
    """A primitive-element presentation of every nonrational reduced involution."""
    T = sp.Symbol("T")
    A = a if C else 1
    # Any pair of distinct solutions excludes at most one choice of slope.
    for c in range(count*(count-1)//2+1):
        extension = sp.groebner([p.as_expr() for p in basis.polys]+[T-A-c*b],
                                *variables, T, domain=sp.QQ)
        minimal = sp.Poly(extension.polys[-1].as_expr(), T, domain=sp.QQ)
        if minimal.degree() == count:
            break
    else:
        raise ArithmeticError("Could not separate the reduced involutions.")
    coords = [sp.Poly(extension.reduce(sp.sympify(v))[1], T, domain=sp.QQ)
              for v in (A, b, C)]
    return [
        {"minimal_polynomial": str(factor.as_expr()),
         "coefficients": [str(coord.rem(factor).as_expr()) for coord in coords]}
        for factor, exponent in minimal.factor_list()[1] if factor.degree() > 1
    ]


def involutions(coefficients):
    """Return the rational involutions and exact number-field presentations."""
    f = sp.Poly.from_list(list(reversed(coefficients)), x, domain=sp.QQ)
    if f.degree() not in (5, 6) or not f.is_sqf:
        raise ValueError("Use a squarefree rational polynomial of degree five or six.")
    solutions, algebraic, geometric_count = [], [], 0
    for A, B, C, variables in ((a, b, 1, (z, b, a)), (1, b, 0, (b,))):
        delta = A*A+B*C
        equations = sp.Poly(sp.expand(sum(
            f.nth(i)*(A*x+B)**i*(C*x-A)**(6-i) for i in range(7)
        )) - delta**3*f.as_expr(), x).all_coeffs()
        if C:
            equations.append(z*delta-1)
        basis = sp.groebner(equations, *variables, domain=sp.QQ)
        if list(basis) == [1]:
            continue
        # The stabilizer injects into permutations of the six branch points.
        # The saturated scheme is reduced, as explained in the proof.
        leads = [p.LM(order=basis.order).exponents for p in basis.polys]
        bounds = [min(m[i] for m in leads if m[i] and sum(m) == m[i])
                  for i in range(len(variables))]
        chart_count = sum(
            all(not all(e >= d for e, d in zip(m, lead)) for lead in leads)
            for m in product(*(range(bound) for bound in bounds))
        )

        def rational_zeros(index, substitution):
            var = variables[index]
            polys = [sp.Poly(p.as_expr().subs(substitution), var, domain=sp.QQ)
                     for p in basis.polys if p.as_expr().subs(substitution).free_symbols <= {var}]
            univariate = reduce(sp.gcd, polys)
            for root in univariate.ground_roots():
                values = {**substitution, var: root}
                if index == 0:
                    solutions.append((values[a], values[b], 1) if C else (1, values[b], 0))
                else:
                    rational_zeros(index-1, values)

        previous = len(solutions)
        rational_zeros(len(variables)-1, {})
        geometric_count += chart_count
        if len(solutions)-previous < chart_count:
            algebraic.extend(algebraic_involutions(basis, variables, C, chart_count))
    return f, solutions, algebraic, geometric_count


def quotient_data(f, involution, field=sp.QQ):
    A, B, C = involution
    num, den = A*x+B, C*x-A
    t = 0
    while (f.eval(t) == 0 or not field.convert(den.subs(x, t))
           or not field.convert((num-t*den).subs(x, t))):
        t += 1
    # The invariant coordinate is 1/(x-t) + 1/(sigma(x)-t).
    N = num-t*den+den*(x-t)
    D = (x-t)*(num-t*den)
    parameter_domain = field.poly_ring(u)
    fibre = sp.Poly(u*D-N, x, domain=parameter_domain)
    orbit = f.set_domain(parameter_domain).resultant(fibre)
    if f.degree() == 5:
        orbit *= u - (0 if not C else 1/(A/C-t))
    orbit = sp.Poly(orbit, u, domain=field).sqf_part().monic()
    fixed = sp.Poly(C*x*x-2*A*x-B, x, domain=parameter_domain).resultant(fibre)
    if not C:
        fixed *= u
    fixed = sp.Poly(fixed, u, domain=field).monic()
    if orbit.degree() != 3 or fixed.degree() != 2 or not fixed.is_sqf:
        raise ArithmeticError("The branch orbit and fixed-point polynomials are inconsistent.")
    return orbit, fixed


def cubic_model(orbit, point, field):
    """The roots of the monic cubic are 1/(label - point)."""
    orbit = orbit.set_domain(field)
    X = sp.Poly(x, x, domain=field)
    linear = sp.Poly.from_list([field.convert(point), field.one], x, domain=field)
    model = sum((X**(3-i)*linear**i).mul_ground(orbit.rep.nth(i)) for i in range(4))
    return model.monic()


def short_model(cubic):
    _, A, B, C = cubic.rep.to_list()
    field = cubic.domain
    three = field.convert(3)
    p = B-A*A/three
    q = C-A*B/three+2*A**3/three**3
    return field.to_sympy(A/three), sp.Poly.from_list([field.one, field.zero, p, q], x, domain=field)


def invariant(cubic):
    p, q = cubic.rep.nth(1), cubic.rep.nth(0)
    return cubic.domain.to_sympy(1728*4*p**3/(4*p**3+27*q**2))


def nonintegral(value, field):
    if field == sp.QQ:
        return not value.is_Integer
    t = sp.Symbol("t")
    modulus = sp.Poly.from_list(field.mod.to_list(), t, domain=sp.QQ)
    polynomial = sp.Poly.from_list(field.convert(value).to_list(), t, domain=sp.QQ)
    characteristic = sp.resultant(modulus.as_expr(), z-polynomial.as_expr(), t)
    minimal = sp.Poly(characteristic, z, domain=sp.QQ).sqf_part().monic()
    return any(not coefficient.is_Integer for coefficient in minimal.all_coeffs())


def non_cm_witness(value, field):
    if value.is_Rational:
        return None if value in RATIONAL_CM else {"test": "rational singular moduli"}
    if nonintegral(value, field):
        return {"test": "nonintegral j"}

    t = sp.Symbol("t")
    modulus = sp.Poly.from_list(field.mod.to_list(), t, domain=sp.QQ)
    coordinates = sp.Poly.from_list(field.convert(value).to_list(), t, domain=sp.QQ)
    witnesses = {}
    for p in map(int, sp.primerange(5, 150)):
        if any(c.q % p == 0 for poly in (modulus, coordinates) for c in poly.all_coeffs()):
            continue

        def mod_p(poly):
            values = [int(c.p)*pow(int(c.q), -1, p) % p for c in poly.all_coeffs()]
            return sp.Poly.from_list(values, t, modulus=p)

        over_fp = mod_p(modulus)
        if not over_fp.is_sqf:
            continue
        image = mod_p(coordinates)
        squares = {a*a % p for a in range(1, p)}
        for root in sorted(over_fp.ground_roots().keys()):
            j = int(image.eval(root)) % p
            if j in (0, 1728 % p):
                continue
            unit = j*pow(1728-j, -1, p) % p
            values = [(x**3+3*unit*x+2*unit) % p for x in range(p)]
            a = -sum(1 if v in squares else -1 for v in values if v)
            if not a:
                continue
            discriminant = -int(sp.prod(prime for prime, e in sp.factorint(4*p-a*a).items() if e % 2))
            witnesses.setdefault(discriminant, {"p": p, "root": int(root) % p,
                                               "trace": a, "quadratic_field": discriminant})
            if len(witnesses) == 2:
                return {"test": "distinct ordinary endomorphism fields", "reductions": list(witnesses.values())}
    return None


def horizontal_scale(source, target):
    field = source.domain if source.domain != sp.QQ else target.domain
    source, target = source.set_domain(field), target.set_domain(field)
    a, b = source.rep.nth(1), source.rep.nth(0)
    c, d = target.rep.nth(1), target.rep.nth(0)
    return field.to_sympy(d*a/(c*b))


def pari_isogeny(j1, j2, field, gp):
    """Find the minimum-degree kernel over the common field of the invariants."""
    t = sp.Symbol("t")
    if field == sp.QQ:
        field_setup, field_arg = "", ""
    else:
        modulus = str(sp.Poly.from_list(field.mod.to_list(), t).as_expr()).replace("**", "^")
        field_setup = f"K=nfinit({modulus});"
        field_arg = ",K"

    def gp_expression(value):
        if field == sp.QQ:
            return str(value)
        expression = sp.Poly.from_list(value.to_list(), t, domain=sp.QQ).as_expr()
        return "Mod(" + str(expression).replace("**", "^") + ",K.pol)"

    j1, j2 = field.convert(j1), field.convert(j2)
    c = 3*j1/(1728-j1)
    d = 2*j1/(1728-j1)
    script = f"""
pack(v) = Str(Vecrev(lift(v)));
{{
iferr(
    {field_setup}
    E=ellinit([{gp_expression(c)},{gp_expression(d)}]{field_arg}); S=ellisomat(E);
    for(i=1,#S[1],R=S[1][i]; J=ellinit(R[1]{field_arg}).j;
        if(J=={gp_expression(j2)},U=R[2];
            print(S[2][1,i],"|",pack(R[1][1]),"|",pack(R[1][2]),"|",
                apply(pack,Vecrev(U[1])),"|",apply(pack,Vecrev(U[3]^2))))),
    err, print("ERROR ",err); quit(1)
);
quit(0)
}}
"""
    proc = subprocess.run([gp, "-q", "-f", "-s", "128M"], input=script, text=True,
                          capture_output=True, timeout=90)
    if proc.returncode or proc.stderr.strip():
        raise ArithmeticError(proc.stderr.strip() or proc.stdout.strip() or "PARI did not certify the isogeny.")

    def coefficient(value):
        interior = value.strip()[1:-1]
        values = list(map(sp.Rational, interior.split(","))) if interior else []
        if field == sp.QQ:
            return field.convert(values[0]) if values else field.zero
        result = field.zero
        for c in reversed(values):
            result = result*field.unit+field.convert(c)
        return result

    def cubic(coefficient_x, constant):
        return sp.Poly.from_list([field.one, field.zero, coefficient_x, constant], x, domain=field)

    def polynomial(vector):
        # GP's coefficient lists are quoted inside the outer vector.
        vectors = json.loads(vector)
        cs = [coefficient(value) for value in vectors]
        return sp.Poly.from_list(cs[::-1], x, domain=field)

    rows = []
    for line in proc.stdout.splitlines():
        k, p, q, numerator, denominator = line.split("|")
        target = cubic(coefficient(p), coefficient(q))
        rows.append((int(k), target, polynomial(numerator), polynomial(denominator)))
    if not rows:
        return None
    k, target, N, D = min(rows, key=lambda row: row[0])
    source = cubic(c, d)
    common = N.gcd(D)
    N, D = N.exquo(common), D.exquo(common)
    if N.degree() != k or D.degree() != k-1 or N.gcd(D).degree() != 0:
        raise ArithmeticError("The isogeny x-map has the wrong degree or poles.")
    left = source*(N.diff()*D-N*D.diff())**2
    right = D*(N**3+(N*D**2).mul_ground(target.rep.nth(1))
               +(D**3).mul_ground(target.rep.nth(0)))
    if left.mul_ground(right.rep.LC()) != right.mul_ground(left.rep.LC()):
        raise ArithmeticError("The isogeny x-map fails the exact pullback identity.")
    return k, source, target, N, D

def split_fixed_images(orbit, fixed):
    """Find both fixed labels over one explicitly presented number field."""
    field = fixed.domain
    points = list(fixed.ground_roots())
    if len(points) == 2:
        return orbit, field, sorted(points, key=sp.default_sort_key)
    if field == sp.QQ:
        points = fixed.all_roots()
        return orbit, field.algebraic_field(points[0]), points

    T, W, Z = sp.symbols("T W Z")
    modulus = sp.Poly.from_list(field.mod.to_list(), T, domain=sp.QQ)
    size = 2*modulus.degree()

    def old_coordinates(value):
        coefficients = value.to_list()
        return sp.Poly.from_list(coefficients, T, domain=sp.QQ)

    _, linear, constant = fixed.rep.to_list()
    relative = W**2+old_coordinates(linear).as_expr()*W+old_coordinates(constant).as_expr()
    for slope in range(size*(size-1)//2+1):
        basis = sp.groebner([modulus.as_expr(), relative, Z-W-slope*T], T, W, Z, domain=sp.QQ)
        primitive = sp.Poly(basis.polys[-1].as_expr(), Z, domain=sp.QQ)
        if primitive.degree() == size:
            break
    else:
        raise ArithmeticError("Could not separate the fixed labels over the involution field.")
    if not primitive.is_irreducible:
        raise ArithmeticError("The fixed-label extension should be a field.")
    root = sp.CRootOf(primitive, 0)
    extended = sp.QQ.algebraic_field(root)
    def evaluate(poly, point):
        value = extended.zero
        for coefficient in poly.rep.to_list():
            value = value*point + extended.convert(coefficient, poly.domain)
        return value

    old_root, p1 = (evaluate(sp.Poly(basis.reduce(var)[1], Z, domain=sp.QQ), extended.unit)
                    for var in (T, W))

    def transfer(poly):
        return sp.Poly.from_list([evaluate(old_coordinates(c), old_root) for c in poly.rep.to_list()],
                                 u, domain=extended)

    orbit, fixed = transfer(orbit), transfer(fixed)
    p2 = -fixed.rep.nth(1)-p1
    if any(evaluate(fixed, point) for point in (p1, p2)):
        raise ArithmeticError("A proposed fixed label fails its exact defining equation.")
    return orbit, extended, tuple(extended.to_sympy(point) for point in (p1, p2))


def torsion_action(orbit, fixed, gp):
    orbit, field, (p1, p2) = split_fixed_images(orbit, fixed)
    shift1, source = short_model(cubic_model(orbit, p1, field))
    shift2, target = short_model(cubic_model(orbit, p2, field))
    j1, j2 = invariant(source), invariant(target)
    witness = non_cm_witness(j1, field)
    if witness is None:
        witness = non_cm_witness(j2, field)
        if witness is None:
            return {"j": [str(j1), str(j2)], "status": "non-CM not certified"}
        p1, p2 = p2, p1
        j1, j2 = j2, j1
        source, target = target, source
        shift1, shift2 = shift2, shift1
    info = {"j": [str(j1), str(j2)], "non_cm_witness": witness}
    if field != sp.QQ:
        info["label_field"] = [str(c) for c in reversed(field.mod.to_list())]
    if j1 == j2:
        k = 1
        first, second = source, source
        scale1 = 1
        scale2 = horizontal_scale(source, target)
        N, D = sp.Poly(x, x), sp.Poly(1, x)
    else:
        isogeny_field = sp.QQ if j1.is_Rational and j2.is_Rational else field
        result = pari_isogeny(j1, j2, isogeny_field, gp)
        if result is None:
            return {**info, "status": "geometric elliptic factors are nonisogenous"}
        k, first, second, N, D = result
        scale1, scale2 = horizontal_scale(source, first), horizontal_scale(second, target)

    label = sp.Poly(u, u, domain=field)

    def constant(value, value_field=None):
        return sp.Poly.from_list([field.convert(value, value_field)], u, domain=field)

    source_x = ((sp.invert(label-constant(p1), orbit)+constant(shift1)).mul_ground(scale1)).rem(orbit)

    def evaluate(poly):
        value = constant(0)
        for term in poly.rep.to_list():
            value = (value*source_x+constant(term, poly.domain)).rem(orbit)
        return value

    numerator, denominator = evaluate(N), evaluate(D)
    image_numerator = (numerator.mul_ground(scale2)-denominator.mul_ground(shift2)).rem(orbit)
    # The label map is p2 + denominator / image_numerator. At a kernel
    # point it gives p2, which is not a root of the orbit cubic.
    image = (constant(p2)+denominator*sp.invert(image_numerator, orbit)).rem(orbit)
    # Check the image of every non-kernel label in the target branch polynomial.
    labels = orbit.exquo(orbit.gcd(denominator))
    if not orbit.compose(image).rem(labels).is_zero:
        raise ArithmeticError("A non-kernel label did not map to a target label.")
    fixed_count = orbit.gcd(image-label).degree()
    info.update({"k": k, "fixed_nonzero_two_torsion": fixed_count,
                 "fixed_images_in_order": [str(p1), str(p2)],
                 "quotient_x_shifts": [str(shift1), str(shift2)],
                 "isogeny_x_scales": [str(scale1), str(scale2)],
                 "isogeny_models": [str(cubic.as_expr()) for cubic in (first, second)],
                 "torsion_label_map": str(image.as_expr()),
                 "isogeny_x_numerator": str(N.as_expr()),
                 "isogeny_x_denominator": str(D.as_expr())})
    if k % 2:
        info["pi_type"] = {0: "3-cycle", 1: "transposition", 3: "identity"}[fixed_count]
    else:
        kernel_label = orbit.gcd(denominator)
        if kernel_label.degree() != 1 or fixed_count not in (0, 1):
            raise ArithmeticError("A cyclic even isogeny must have a single nonzero two-torsion kernel point.")
        info["kernel_label"] = str(-kernel_label.nth(0))
        info["gluing_sends_kernel_to_image"] = fixed_count == 0
    return info


def kernel_order(a, b, k, fixed_count):
    """Kernel of the primitive slope under the product-to-Jacobian map."""
    if gcd(a, b) != 1:
        raise ValueError("Use a primitive integer direction.")
    m = (k & -k).bit_length()-1
    g = gcd(a, k >> m)
    if not m:
        return g*(1+fixed_count if a % 2 and b % 2 else 1)
    if a % 2:
        return g*(1+fixed_count if b % 2 else 1)
    d = (a & -a).bit_length()-1 if a else m+1
    L = min(d, m)
    extra = (d < m and not fixed_count) or (d == m and fixed_count)
    return g << (L+bool(extra))


def degree(a, b, k, fixed_count):
    numerator = 2*(a*a+k*b*b)
    order = kernel_order(a, b, k, fixed_count)
    quotient, remainder = divmod(numerator, order)
    if remainder:
        raise ArithmeticError("The gluing order does not divide the pulled-back degree.")
    return quotient


def spectrum(k, fixed_count, bound):
    # The kernel of the product isogeny has order four, and the kernel of
    # the slope into the product has order gcd(a,k).
    counts = Counter({2: 1}) if bound >= 2 else Counter()
    for a in range(1, isqrt(2*k*bound)+1):
        window = isqrt((2*k*bound-a*a)//k)
        for b in range(-window, window+1):
            if gcd(a, b) == 1:
                d = degree(a, b, k, fixed_count)
                if d <= bound:
                    counts[d] += 1
    return dict(sorted(counts.items()))


def leading_constant(k, fixed_count):
    """Return the leading coefficient as a rational multiple of 1/(pi sqrt(k))."""
    m = (k & -k).bit_length()-1
    two = (3+fixed_count if not m else (m+5 if fixed_count else 2*m+2))
    return sp.Rational(two, 2)*sp.prod(
        1+sp.Rational(e*(p-1), p+1) for p, e in sp.factorint(k >> m).items()
    )


def reconstruct_quotient(f, involution, field, bound, gp):
    orbit, fixed = quotient_data(f, involution, field)
    record = torsion_action(orbit, fixed, gp)
    record.update({"involution": list(map(str, involution)), "orbit_cubic": str(orbit.as_expr()),
                   "fixed_image_poly": str(fixed.as_expr())})
    if "k" in record:
        k, fixed_count = record["k"], record["fixed_nonzero_two_torsion"]
        c = leading_constant(k, fixed_count)
        record.update({"spectrum": spectrum(k, fixed_count, bound),
                       "leading_numerator": str(c),
                       "leading_constant": float(c)/(pi*sqrt(k)), "status": "certified"})
    elif record["status"] == "geometric elliptic factors are nonisogenous":
        record["spectrum"] = {2: 2} if bound >= 2 else {}
    return record


def reconstruct(coefficients, bound=200, gp="gp"):
    if not isinstance(coefficients, list) or any(type(c) is not int for c in coefficients):
        raise ValueError("Use a list of integer coefficients.")
    if bound < 0:
        raise ValueError("Choose a nonnegative degree bound.")
    f, invols, algebraic, geometric_count = involutions(coefficients)
    results = [reconstruct_quotient(f, involution, sp.QQ, bound, gp) for involution in invols]
    number_field_results = []
    T = sp.Symbol("T")
    for presentation in algebraic:
        defining_polynomial = sp.Poly(presentation["minimal_polynomial"], T, domain=sp.QQ)
        root = sp.CRootOf(defining_polynomial, 0)
        field = sp.QQ.algebraic_field(root)
        involution = tuple(sp.Poly(c, T, domain=sp.QQ).eval(root) for c in presentation["coefficients"])
        record = reconstruct_quotient(f, involution, field, bound, gp)
        record["involution_field"] = presentation["minimal_polynomial"]
        record["conjugate_involutions"] = defining_polynomial.degree()
        number_field_results.append(record)
    return {"f": coefficients, "geometric_reduced_involutions": geometric_count,
            "rational_reduced_involutions": len(invols), "algebraic_reduced_involutions": algebraic,
            "quotients": results, "number_field_quotients": number_field_results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coefficients", type=json.loads, help="JSON array in ascending coefficient order")
    parser.add_argument("--bound", type=int, default=200)
    parser.add_argument("--gp", default="gp", help="PARI/GP executable")
    args = parser.parse_args()
    print(json.dumps(reconstruct(args.coefficients, args.bound, args.gp), indent=2))
