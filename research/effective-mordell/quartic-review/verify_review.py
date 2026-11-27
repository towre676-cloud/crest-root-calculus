"""Independent exact checks of six archived records and three further genus-two inputs.

Run python verify_review.py with SymPy and cypari2 installed. Coordinates in
review-inputs.json are those of Claude's paired-fibre chart, as labelled there.
"""
import json
from math import isqrt
from pathlib import Path

import sympy as sp
from cypari2 import Pari


pari = Pari()
pari.allocatemem(128_000_000)
x, t = pari("x"), pari("t")
w = pari('w')
z = sp.Symbol("x")
weights = (1, 2, 3, 5, 7, 11)


def require(claim, message):
    if not claim:
        raise ValueError(message)


def poly_in_x(coeffs):
    return sum(int(c)*x**k for k, c in enumerate(coeffs))


def check_orbits(record, f):
    seen = []
    for orbit in record["pairing_orbits"]:
        field = orbit["field"]
        unit = 1 if field == "Q" else pari.Mod(1, pari(field))
        sigma = pari(orbit["involution"])*unit
        at_sigma = pari.Mod(sigma, f)
        require(pari.subst(f, x, at_sigma) == 0, "branch permutation")
        require(pari.subst(sigma, x, at_sigma) == pari.Mod(x, f), "involution")
        require(pari.poldegree(pari.gcd(x-sigma, f), x) == 0, "fixed branch")
        key = sum(c*pari.polcoef(sigma, i, x) for i, c in enumerate(weights))
        mp = pari.minpoly(key, t)
        require(pari.poldegree(mp) == orbit["orbit_size"], "orbit degree")
        require(all(pari.poldegree(pari.gcd(mp, q)) == 0 for q in seen),
                "pairing orbits coincide")
        seen.append(mp)
    require(sum(int(pari.poldegree(q)) for q in seen) == 15, "missing pairing")


def check_maps(record, f):
    for m in record["maps"]:
        field = m["point_field"]
        unit = 1 if field == "Q" else pari.Mod(1, pari(field))
        n, s, h, g, omega = [pari(m[k])*unit for k in
                             ("numerator", "half_poles", "lift", "cubic", "differential")]
        a, b = pari.polcoef(g, 1, t), pari.polcoef(g, 0, t)
        require(g == t**3+a*t+b and 4*a**3+27*b**2 != 0, "target smoothness")
        require(pari.poldegree(pari.gcd(n,s),x) == 0 and
                max(pari.poldegree(n,x),2*pari.poldegree(s,x)) == 4, "cover degree")
        require(n**3+a*n*s**4+b*s**6 == f*h**2, "cover identity")
        require(pari.deriv(n,x)*s-2*n*pari.deriv(s,x) == omega*h and
                omega != 0 and pari.poldegree(omega,x)<=1, "differential")
        require(1728*4*a**3/(4*a**3+27*b**2) == pari(m["j"]), "j-invariant")


def rational_square(v):
    q = sp.Rational(v)
    return q >= 0 and isqrt(q.p)**2 == q.p and isqrt(q.q)**2 == q.q


def local_square(n, p):
    if n == 0:
        return True
    exponent = 0
    while n % p == 0:
        n //= p
        exponent += 1
    unit_square = n % 8 == 1 if p == 2 else pow(n % p, (p-1)//2, p) == 1
    return exponent % 2 == 0 and unit_square


def check_twins(rec):
    m = next(m for m in rec["maps"] if m["cubic"] == "t^3 - 105903*t - 11165202")
    E = pari.ellinit([-105903, -11165202])
    require(list(E.ellrank())[:2] == [0, 0] and E.elltors()[0] == 4, "target group")
    n, s = [sp.sympify(m[key].replace("^","**")) for key in ("numerator","half_poles")]
    f = sum(c*z**i for i,c in enumerate(rec["coefficients"]))
    require(not rational_square(rec["coefficients"][-1]), "point at infinity")
    require(not sp.ground_roots(s,z), "rational pole")
    roots = [v for v in (-246, -123, 369) if v**3-105903*v-11165202 == 0]
    require(len(roots) == 3, "three remaining target points")
    candidates = {q for v in roots for q in sp.ground_roots(n-v*s**2,z)}
    require(all(not rational_square(f.subs(z,q)) for q in candidates), "rational lift")
    # Good primes >= 17 have rational smooth reduction points by Hasse–Weil.
    # Positive leading coefficient gives a real point. Test all other places.
    bad = set(sp.factorint(abs(sp.discriminant(f,z)))) | {2}
    require(bad == {2,3,5,41}, "bad-prime candidates")
    places = (2,3,5,7,11,13,41)
    witnesses = (0,11,2,5,0,3,3) if rec["name"]=="T1" else (3,1,0,2,1,1,0)
    for p, witness in zip(places,witnesses):
        require(local_square(int(f.subs(z,witness)), p), "local witness failed")
    require(sp.LC(sp.Poly(f,z)) > 0, "real witness")
    print(rec["name"], "has local points everywhere and no rational point.",
          "Finite candidate abscissas:", sorted(candidates))


def check_stoll(rec):
    E = pari.ellinit([-108, 0])
    require(list(E.ellrank())[:2] == [0, 0] and E.elltors()[0] == 2, "Stoll target")
    m = next(m for m in rec["maps"] if m["cubic"] == "t^3 - 108*t")
    n,s,h = [sp.sympify(m[k].replace("^","**")) for k in ("numerator","half_poles","lift")]
    u = 1/(z-1)
    horiz = sp.factor((n/s**2).subs(z,u))
    vert = sp.factor((h/s**3).subs(z,u)*u**3)
    require(sp.cancel(horiz+32*z/(z**2+1)**2) == 0 and
            sp.cancel(vert-8*(3*z**2-1)/(z**2+1)**3) == 0,
            "back-substitution in original Stoll coordinates")
    f = 6*z**5+28*z**3+54*z
    require(sp.cancel(horiz**3-108*horiz-f*vert**2)==0, "Stoll original identity")
    # Denominator has no rational zero, and the only finite rational target is (0,0).
    # Thus a finite rational source point needs x=0; the odd model has one rational infinity.
    print("Stoll 5.3(4) has exactly (0,0) and the point at infinity.")


def check_nonstable_pairing(rec):
    f = x**6+1
    K = pari("w^2-w+1")
    unit = pari.Mod(1,K)
    # The pairing polynomial has degree < 6, so unequal conjugate
    # coefficient vectors give different permutations on the six roots.
    sigma = ((w-1)*x**5-2*w*x**3-2*x)*unit/3
    sigma_bar = pari.subst(pari.lift(sigma),w,1-w)*unit
    require(sigma != sigma_bar, "pairing is rational")
    selected = [m for m in rec["maps"] if m["pairing_field"] != "Q"]
    require(len(selected)==4 and rec["geometric_maps"]==12, "nonstable orbit counts")
    poles = set()
    for m in rec["maps"]:
        s = pari(m["half_poles"])*unit
        poles.update((str(s), str(pari.subst(pari.lift(s),w,1-w)*unit)))
    require(len(poles)==12, "conjugate half-poles coincide")
    for m in selected:
        s,n = [pari(m[key])*unit for key in ("half_poles","numerator")]
        r = pari.Mod(n,f)/pari.Mod(s**2,f)
        arg = pari.Mod(sigma,f)
        require(pari.subst(n,x,arg)/pari.subst(s,x,arg)**2 == r, "incorrect pairing")
        arg_bar = pari.Mod(sigma_bar,f)
        require(pari.subst(n,x,arg_bar)/pari.subst(s,x,arg_bar)**2 != r,
                "conjugation preserved the pairing")
        branch = pari.polresultant(f,t*s**2-n,x)
        sqfree = branch/pari.gcd(branch,pari.deriv(branch,t))
        require(pari.poldegree(sqfree,t)==3 and
                branch/pari.pollead(branch,t) == (sqfree/pari.pollead(sqfree,t))**2,
                "three different branch values with two Weierstrass points each")
    print("Four quadratic orbits on y^2=x^6+1 have a genuinely nonstable pairing.")


def check_twists():
    d = sp.Symbol("d")
    f = 6*z**5+28*z**3+54*z
    maps = [
        (-108,-32*z/(z**2+1)**2,8*(3*z**2-1)/(z**2+1)**3),
        (-12,-32*z**3/(z**2+9)**2,-8*z*(z**2-27)/(z**2+9)**3),
    ]
    for a,horizontal,vertical in maps:
        X = d*horizontal
        require(sp.cancel(X**3+a*d**2*X-d*f*(d*vertical)**2)==0, "twist identity")
    open_parameters = []
    for value in range(1,31):
        if any(value%(i*i)==0 for i in range(2,isqrt(value)+1)):
            continue
        ranks = []
        for a,_,_ in maps:
            curve = pari.ellinit([a*value**2,0])
            low,high = list(curve.ellrank())[:2]
            require(low == high and curve.elltors()[0]==2, "rank or torsion changed")
            ranks.append(int(low))
        if min(ranks)>0:
            open_parameters.append(value)
        print("d =",value,"ranks =",ranks)
    require(open_parameters == [14,22], "exceptional parameters")
    # Either rank-zero target has only (0,0) in its affine chart.
    # Both horizontal functions have only the finite rational zero z=0
    # and no rational poles. Each odd-degree source has one point at infinity.
    print("17 of 19 twist parameters have exactly (0,0) and infinity by these maps.")


if __name__ == "__main__":
    records = json.loads(Path(__file__).with_name("review-inputs.json").read_text())
    for rec in records:
        f = poly_in_x(rec["coefficients"])
        center = rec["reciprocal_chart_center"]
        if center is not None:
            f = sum(c*(center*x+1)**i*x**(6-i) for i,c in enumerate(rec["coefficients"]))
        check_orbits(rec, f)
        check_maps(rec, f)
        if rec["name"] in ("T1","T2"):
            check_twins(rec)
        if rec["name"] == "Stoll-5.3(4)":
            check_stoll(rec)
        if rec["name"] == "Stoll-d26":
            targets = {str(pari(m["cubic"])) for m in rec["maps"]}
            require(targets == {"t^3 - 4563*t","t^3 - 507*t"}, "cold d=26 targets")
            for a,r in ((-4563,0),(-507,2)):
                curve = pari.ellinit([a,0])
                require(list(curve.ellrank())[:2]==[r,r] and curve.elltors()[0]==2,
                        "cold d=26 arithmetic")
        if rec["name"] == "hexagonal":
            check_nonstable_pairing(rec)
    for b, expect in ((-63269478,0),(63269478,1)):
        require(list(pari.ellinit([196677,b]).ellrank())[:2] == [expect,expect],
                "quadratic-field rank decomposition")
    print("Exact orbit disjointness, map identities and differentials passed for all nine inputs.")
    check_twists()
