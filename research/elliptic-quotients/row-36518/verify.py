"""Classify elliptic map degrees for Sutherland row 36518 and its family."""
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "row-23640"))
from verify_involution import SmallField, require


def count_mod_five(n):
    field = SmallField(5, n)
    add, scale, mul = field.add, field.scale, field.mul
    roots = [[] for _ in range(field.order)]
    for t in range(field.order):
        roots[mul(t, t)].append(t)

    def evaluate(coefficients, x):
        value = 0
        for c in reversed(coefficients):
            value = add(mul(value, x), c % 5)
        return value

    c_count, e_count = 0, 0
    for u in range(field.order + 1):
        if u == field.order:  # [u:v] = [1:0]
            a, b = 4, (-2) % 5
        else:
            a = evaluate([2, -8, 4], u)
            b = evaluate([-3, 8, -4, 2, -2], u)
        d = add(mul(a, a), scale(-4, b))
        for w in roots[d]:
            e_count += 1
            square = scale(3, add(w, scale(-1, a)))  # 1/2 = 3 mod 5
            c_count += len(roots[square])
    return c_count, e_count


def rational_sieve():
    primes = [5, 7, 11, 13, 17, 19, 37, 43, 67, 163]
    excluded, witnesses = set(), []
    for p in (5, 7, 11, 17, 41):
        squares = [0]*p
        for y in range(p):
            squares[(y*y) % p] += 1
        count = 1 + sum(squares[(x**3+16) % p] for x in range(p))
        trace = p+1-count
        delta = trace**2-4*p
        disallowed = [ell for ell in primes
                      if ell != p and pow(delta, (ell-1)//2, ell) == ell-1]
        excluded.update(disallowed)
        witnesses.append({"p": p, "trace": trace, "nonsquare_mod": disallowed})
    require(excluded == set(primes), "The prime-isogeny sieve has a gap.")
    x = sp.symbols("x")
    require(sp.Poly(x**3+16, x).is_irreducible, "A rational 2-isogeny survives.")
    require(sp.expand(x**3+16-4**2) == x**3, "The flex at (0,4) failed.")
    return witnesses


def main():
    x, y, z, u, v, t, w, k, T, U = sp.symbols("x y z u v t w k T U")
    terms = json.loads((HERE / "quartic.json").read_text(encoding="utf-8"))["terms"]
    original = sum(c*x**i*y**j*z**l for i, j, l, c in terms)
    a = 2*v**2-8*u*v+4*u**2
    b = -3*v**4+8*u*v**3-4*u**2*v**2+2*u**3*v-2*u**4
    r = t**2+5*k*u*v
    f = r**2+a*r+b
    change = original.subs({x:u-2*v, y:u, z:v+t}, simultaneous=True)
    require(sp.expand(change-f.subs(k, 0)) == 0, "The row coordinate change failed.")
    require(f.subs({u: 0, v: 1, t: 1}) == 0, "Missing rational source point.")
    difference = sp.Poly(f-f.subs(k, 0), t, u, v, k)
    require(all(c % 5 == 0 for c in difference.coeffs()), "Family reduction varies at 5.")
    d = sp.expand(a**2-4*b)
    require(sp.expand((2*r+a)**2-d-4*f) == 0, "The fixed elliptic quotient identity failed.")
    # In homogeneous coordinates X=2u/(u-v), Y=w/(u-v)^2.
    require(sp.expand(d - (2*u)**3*(u-v) - 16*(u-v)**4) == 0,
            "The rational Mordell model failed.")
    coords = (u, v, t)
    gradients = [f.subs(k, 0).diff(q) for q in coords]
    for q in coords:
        gb = sp.groebner([g.subs(q, 1) for g in gradients], *(h for h in coords if h != q), modulus=5)
        require(list(gb) == [1], "The plane quartic is singular modulo 5.")

    counts = [count_mod_five(n) for n in (1, 2, 3)]
    c_counts, e_counts = map(list,zip(*counts))
    require(c_counts == [8, 40, 140] and e_counts == [6, 36, 126], "Point counts disagree.")
    s1, s2, s3 = [5**n+1-N for n, N in enumerate(c_counts, 1)]
    e2 = sp.Rational(s1*s1-s2, 2)
    e3 = sp.Rational(s1**3-3*s1*s2+2*s3, 6)
    frob = T**6-s1*T**5+e2*T**4-e3*T**3+5*e2*T**2-25*s1*T+125
    h = T**4+2*T**3+4*T**2+10*T+25
    require(sp.expand(frob-(T**2+5)*h) == 0, "Failed to identify the complementary surface.")
    q1 = 5*U**4+6*U**3+6*U**2+6*U+5
    q2 = 25*U**4+20*U**3+26*U**2+20*U+25
    resultant = sp.resultant(h, h.subs(T,U*T), T)
    require(sp.expand(resultant-625*(U-1)**4*q1**2*q2) == 0, "Root-ratio resultant failed.")
    for expr, var in ((h,T),(q1,U),(q2,U)):
        fz = sp.Poly(expr, var)
        fm = sp.Poly(expr, var, modulus=11)
        require(fz.content()==1 and fm.degree()==4 and fm.is_irreducible,
                "Modulo-11 irreducibility witness failed.")
    require(all(sp.Poly(q, U).LC()>1 for q in (q1, q2)), "Root-ratio factor is monic.")
    print(json.dumps({"row": 36518, "p": 5, "C_counts": c_counts, "E_counts": e_counts,
                      "elliptic_curve": "Y^2=X^3+16", "j": 0,
                      "prym_frobenius": str(h), "nontrivial_root_ratio_factors": [str(q1), str(q2)],
                      "coefficient_prime": 11, "rational_isogeny_sieve": rational_sieve(),
                      "geometric_degrees_all_elliptic_targets": "2*n",
                      "rational_degrees_all_elliptic_targets": "2*n^2 or 6*n^2",
                      "family": "r=t^2+5*k*u*v, k in Z"}, indent=2))


if __name__ == "__main__":
    main()
