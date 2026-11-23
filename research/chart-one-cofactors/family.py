"""A family of transverse cubic elliptic quotients with complete degree spectra."""
from pathlib import Path
import sys
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"elliptic-quotients/polarized-reduction"))
from screen import curve_profile

t,z = sp.symbols("t z")
MAZUR = (2,3,5,7,11,13,17,19,37,43,67,163)


def source(k):
    square = t**3+t**2+t*z+1+z
    K = z**3+k*z**2+2*z+1
    G = sp.Poly(square**2-K,t,z)
    rows = []
    for (r,n),a in G.terms():
        j = r+n-2
        i = 4-j-n
        if min(i,j,n)<0:
            raise ArithmeticError("The specified plane model does not exist.")
        rows.append([i,j,n,(-1)**(i+j)*int(a)])
    return {"terms":rows}


def invariants(k):
    # Standard Weierstrass invariants for y²=x³+kx²+2x+1.
    cubic = z**3+k*z**2+2*z+1
    disc = 16*sp.discriminant(cubic,z)
    c4 = 16*(k*k-6)
    return disc,sp.cancel(c4**3/disc)


def verify_family():
    data = source(6)
    simplified = {tuple(k[:3]):k[3] for k in data["terms"]}
    profile = curve_profile(simplified,5)
    expected = [1,-1,12,-10,60,-25,125]
    if profile["point_counts"]!=[5,49,131] or profile["weil"]!=expected:
        raise ArithmeticError("The smooth reduction/point-count receipt changed.")
    blocks = profile["elliptic_blocks"]
    if len(blocks)!=1 or blocks[0]["geometric_multiplicity"]!=1 or blocks[0]["trace"]!=0:
        raise ArithmeticError("The special fibre did not bound elliptic directions by one.")
    k = sp.Symbol("k")
    disc,j = invariants(k)
    if not sp.diff(j,k) or disc.subs(k,6)%5==0:
        raise ArithmeticError("The family must vary geometrically and have good reduction.")
    forbidden = set()
    profiles = []
    for prime in (5,7,11,13):
        if disc.subs(k,6)%prime==0:
            raise ArithmeticError("The elliptic sieve used a bad prime.")
        trace = -sum(sp.legendre_symbol(a**3+6*a*a+2*a+1,prime) for a in range(prime))
        excluded = [ell for ell in MAZUR if ell!=prime and
                    (trace%2==1 if ell==2 else sp.legendre_symbol(trace**2-4*prime,ell)==-1)]
        forbidden.update(excluded)
        profiles.append([prime,int(trace),excluded])
    if set(MAZUR)!=forbidden:
        raise ArithmeticError("Some rational prime-isogeny degrees survived.")
    return {"p5_counts":profile["point_counts"],"p5_weil":expected,
            "geometric_degrees":"3n for k = 1 mod 5","rational_degrees":"3n^2 for k = 6 mod 5005",
            "isogeny_sieve":profiles}
