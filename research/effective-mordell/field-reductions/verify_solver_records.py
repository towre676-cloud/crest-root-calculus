"""Independent cover and pairing checks for the 36-input replay."""

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
