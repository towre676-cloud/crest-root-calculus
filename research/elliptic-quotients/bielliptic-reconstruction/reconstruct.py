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


def quotient_data(f, involution):
    A, B, C = involution
    num, den = A*x+B, C*x-A
    t = 0
    while (f.eval(t) == 0 or den.subs(x, t) == 0
           or num.subs(x, t) == t*den.subs(x, t)):
        t += 1
    N, D = sp.fraction(sp.cancel(1/(x-t)+den/(num-t*den)))
    fibre = u*D-N
    orbit = sp.resultant(f.as_expr(), fibre, x)
    if f.degree() == 5:
        orbit *= u - (0 if not C else 1/(sp.Rational(A, C)-t))
    orbit = sp.Poly(orbit, u, domain=sp.QQ).sqf_part().monic()
    fixed = sp.resultant(C*x*x-2*A*x-B, fibre, x)
    if not C:
        fixed *= u
    fixed = sp.Poly(fixed, u, domain=sp.QQ).monic()
    if orbit.degree() != 3 or fixed.degree() != 2 or not fixed.is_sqf:
        raise ArithmeticError("The branch orbit and fixed-point polynomials are inconsistent.")
    return orbit, fixed


def cubic_model(orbit, point, field):
    """The roots of the monic cubic are 1/(label - point)."""
    expression = sum(orbit.nth(i)*x**(3-i)*(point*x+1)**i for i in range(4))
    return sp.Poly(expression, x, domain=field).monic()


def short_model(cubic):
    _, A, B, C = cubic.all_coeffs()
    p = B-A*A/3
    q = C-A*B/3+2*A**3/27
    return A/3, sp.Poly(x**3+p*x+q, x, domain=cubic.domain)


def invariant(cubic):
    p, q = cubic.nth(1), cubic.nth(0)
    field = cubic.domain
    p, q = field.convert(p), field.convert(q)
    return field.to_sympy(1728*4*p**3/(4*p**3+27*q**2))


def nonintegral(value):
    minimal = sp.Poly(sp.minpoly(value, z), z, domain=sp.QQ).monic()
    return any(not coefficient.is_Integer for coefficient in minimal.all_coeffs())


def horizontal_scale(source, target):
    return (target.nth(0)*source.nth(1))/(target.nth(1)*source.nth(0))


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
        expression = sp.Poly.from_list(field.convert(value).to_list(), t).as_expr()
        return f"Mod({expression},K.pol)"

    c = 3*j1/(1728-j1)
    d = 2*j1/(1728-j1)
    script = f"""
pack(v) = Str(polcoef(lift(v),0),":",polcoef(lift(v),1));
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
                          capture_output=True, check=True, timeout=90)
    if proc.stderr.strip():
        raise ArithmeticError(proc.stderr.strip())

    def coefficient(value):
        constant, linear = map(sp.Rational, value.strip().strip('"').split(":"))
        return constant if field == sp.QQ else constant+linear*field.ext.as_expr()

    def polynomial(vector):
        cs = [coefficient(value) for value in vector.strip()[1:-1].split(",")]
        return sp.Poly.from_list(cs[::-1], x, domain=field)

    rows = []
    for line in proc.stdout.splitlines():
        k, p, q, numerator, denominator = line.split("|")
        target = sp.Poly(x**3+coefficient(p)*x+coefficient(q), x, domain=field)
        rows.append((int(k), target, polynomial(numerator), polynomial(denominator)))
    if not rows:
        return None
    k, target, N, D = min(rows, key=lambda row: row[0])
    source = sp.Poly(x**3+c*x+d, x, domain=field)
    common = N.gcd(D)
    N, D = N.exquo(common), D.exquo(common)
    if N.degree() != k or D.degree() != k-1 or N.gcd(D).degree() != 0:
        raise ArithmeticError("The isogeny x-map has the wrong degree or poles.")
    left = source*(N.diff()*D-N*D.diff())**2
    right = D*(N**3+target.nth(1)*N*D**2+target.nth(0)*D**3)
    if left.mul_ground(right.LC()) != right.mul_ground(left.LC()):
        raise ArithmeticError("The isogeny x-map fails the exact pullback identity.")
    return k, source, target, N, D

def torsion_action(orbit, fixed, gp):
    points = list(fixed.ground_roots())
    if len(points) == 2:
        p1, p2 = sorted(points)
        field = sp.QQ
    else:
        p1, p2 = fixed.all_roots()
        field = sp.QQ.algebraic_field(p1)
    shift1, source = short_model(cubic_model(orbit, p1, field))
    shift2, target = short_model(cubic_model(orbit, p2, field))
    j1, j2 = invariant(source), invariant(target)
    info = {"j": [str(j1), str(j2)]}
    if not (nonintegral(j1) or nonintegral(j2)):
        return {**info, "status": "non-CM not certified by nonintegral j"}
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

    orbit = orbit.set_domain(field)
    def reduce(expression):
        return sp.Poly(expression, u, domain=field).rem(orbit)

    source_x = reduce(scale1*(sp.invert(sp.Poly(u-p1, u, domain=field), orbit).as_expr()+shift1))

    def evaluate(poly):
        value = reduce(0)
        for term in poly.all_coeffs():
            value = (value*source_x+reduce(term)).rem(orbit)
        return value

    numerator, denominator = evaluate(N), evaluate(D)
    image_numerator = (reduce(scale2)*numerator-reduce(shift2)*denominator).rem(orbit)
    # The label map is p2 + denominator / image_numerator. At a kernel
    # point it gives p2, which is not a root of the orbit cubic.
    image = (reduce(p2)+denominator*sp.invert(image_numerator, orbit)).rem(orbit)
    # Check the image of every non-kernel label in the target branch polynomial.
    labels = orbit.exquo(orbit.gcd(denominator))
    if not orbit.compose(image).rem(labels).is_zero:
        raise ArithmeticError("A non-kernel label did not map to a target label.")
    fixed_count = orbit.gcd(image-reduce(u)).degree()
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


def reconstruct(coefficients, bound=200, gp="gp"):
    if not isinstance(coefficients, list) or any(type(c) is not int for c in coefficients):
        raise ValueError("Use a list of integer coefficients.")
    if bound < 0:
        raise ValueError("Choose a nonnegative degree bound.")
    f, invols, algebraic, geometric_count = involutions(coefficients)
    results = []
    for invol in invols:
        orbit, fixed = quotient_data(f, invol)
        record = torsion_action(orbit, fixed, gp)
        record.update({"involution": list(map(str, invol)), "orbit_cubic": str(orbit.as_expr()),
                       "fixed_image_poly": str(fixed.as_expr())})
        if "k" in record:
            k, fixed_count = record["k"], record["fixed_nonzero_two_torsion"]
            c = leading_constant(k, fixed_count)
            record.update({"spectrum": spectrum(k, fixed_count, bound),
                           "leading_numerator": str(c),
                           "leading_constant": float(c)/(pi*sqrt(k)), "status": "certified"})
        results.append(record)
    return {"f": coefficients, "geometric_reduced_involutions": geometric_count,
            "rational_reduced_involutions": len(invols),
            "algebraic_reduced_involutions": algebraic, "quotients": results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("coefficients", type=json.loads, help="JSON array in ascending coefficient order")
    parser.add_argument("--bound", type=int, default=200)
    parser.add_argument("--gp", default="gp", help="PARI/GP executable")
    args = parser.parse_args()
    print(json.dumps(reconstruct(args.coefficients, args.bound, args.gp), indent=2))
