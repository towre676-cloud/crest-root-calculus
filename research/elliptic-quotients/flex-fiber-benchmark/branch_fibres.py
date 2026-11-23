"""Recover a degree-four elliptic map from trigonal branch fibres and certify exhaustion."""
import argparse
from itertools import combinations
import json
from pathlib import Path
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "polarized-reduction"))
from screen import curve_profile, read_terms
sys.path.insert(0, str(HERE.parent / "involution-atlas"))
from classify import classify

X,Y,Z = sp.symbols("X Y Z")
t,s,l,m = sp.symbols("t s l m")
Q = sp.Rational


def projections(form):
    """Take rational points on the line at infinity as centres of trigonal pencils."""
    infinity = form.subs(Z,0)
    for p in sp.Poly(infinity.subs(Y,1),X).ground_roots():
        yield (X-p*Y,Y), form.subs({X:t+p*s,Y:s,Z:1}, simultaneous=True)
    if infinity.subs({X:1,Y:0}) == 0:
        yield (Y,X), form.subs({X:s,Y:t,Z:1}, simultaneous=True)


def solve_fibre(a, b, H):
    """Solve the two normalized linear forms in A-B=H, allowing a constant form."""
    h = H.as_expr()
    left = H.eval(b)/(b-a)*(t-a)*(1+l*(t-b))**3
    right = H.eval(a)/(b-a)*(t-b)*(1+m*(t-a))**3
    polynomial = sp.Poly(sp.cancel((left-right-h)/((t-a)*(t-b))), t)
    gb = sp.groebner(polynomial.all_coeffs(),l,m,order="lex",domain=sp.QQ)
    if not gb.is_zero_dimensional:
        return
    for lv,mv in sp.solve_poly_system(list(gb),l,m):
        if lv.is_Rational and mv.is_Rational:
            g = sp.factor(sp.cancel((left/right).subs({l:lv,m:mv})))
            if max(sp.degree(e,t) for e in sp.fraction(g)) == 4:
                yield g, list(map(str,gb))


def recover(trigonal, stats):
    poly = sp.Poly(trigonal,s)
    if poly.degree() != 3:
        return
    branch = sp.Poly(sp.discriminant(poly.as_expr(),s),t)
    factors = branch.factor_list()[1]
    double = [f for f,e in factors if e==2 and f.degree()==1]
    simple = [f.monic() for f,e in factors if e==1]
    if len(double) != 2 or sum(f.degree() for f in simple)!=6 or branch.degree()!=10:
        return
    a,b = sorted(-f.nth(0)/f.nth(1) for f in double)
    lead, c2, c1, c0 = poly.all_coeffs()
    c2,c1,c0 = (sp.cancel(z/lead) for z in (c2,c1,c0))
    dep_a = sp.factor(c1-c2**2/3)
    dep_b = sp.factor(2*c2**3/27-c1*c2/3+c0)
    if dep_a==0 or dep_b==0:
        return
    for n in range(1,len(simple)+1):
        for group in combinations(simple,n):
            if sum(f.degree() for f in group)!=4:
                continue
            H=sp.prod(group)
            other=sp.prod(f for f in simple if f not in group)
            stats["two_variable_systems"] += 1
            for g, basis in solve_fibre(a,b,H):
                num,den = map(lambda v:sp.Poly(v,t), sp.fraction(g))
                if den.gcd(other).degree()!=0:
                    continue
                beta = (num*den.invert(other)).rem(other)
                if beta.degree()!=0 or beta.nth(0) in (0,1):
                    continue
                value = beta.nth(0)
                shift = sp.sqrt(value)
                if not shift.is_Rational:
                    continue
                for B in (-shift,shift):
                    D = 27*(1+value+2*B)/4
                    if D==0:
                        continue
                    scale = sp.factor((g+B)*dep_a/dep_b)
                    if sp.cancel(scale**2*dep_a+D*g) or sp.cancel(scale**3*dep_b+D*g*(g+B)):
                        continue
                    x = sp.factor(scale*(s+c2/3))
                    # Independent substitution in the original cubic function field.
                    residual = sp.together(x**3-D*g*x-D*g*(g+B)).as_numer_denom()[0]
                    if sp.Poly(residual,s,domain=sp.QQ.frac_field(t)).rem(
                            sp.Poly(trigonal,s,domain=sp.QQ.frac_field(t))):
                        raise ArithmeticError("Failed the target identity")
                    # In the function field, dz/(3x^2-Dz) is a scalar multiple of
                    # the canonical differential represented by a line in the source plane.
                    K = sp.QQ.frac_field(t)
                    modulus = sp.Poly(trigonal,s,domain=K)
                    tangent_num = sp.Poly(trigonal.diff(s)*g.diff(t),s,domain=K)
                    tangent_den = sp.Poly(3*x**2-D*g,s,domain=K)
                    line = (tangent_num*tangent_den.invert(modulus)).rem(modulus).as_expr().expand()
                    if sp.Poly(line,t,s).total_degree() > 1:
                        raise ArithmeticError("Pullback failed the canonical line check")
                    # Distinct 0, infinity, 1, beta and nonzero B,D ensure a smooth target.
                    yield {"quartic_fibre":str(H.as_expr()), "other_branch_pair":str(other.as_expr()),
                           "triple_branch_positions":[str(a),str(b)], "elimination_basis":basis,
                           "D":str(D),"B":str(B),"x":str(x),"z":str(g),"degree":4,
                           "canonical_line_in_pencil_coordinates":str(line)}


def discover(path, prime):
    terms = read_terms(path)
    form = sum(c*X**i*Y**j*Z**k for (i,j,k),c in terms.items())
    profile = curve_profile(terms,prime)
    stats = {"rational_projections":0,"two_variable_systems":0}
    maps = []
    for coords, curve in projections(form):
        stats["rational_projections"] += 1
        for recovered in recover(curve,stats):
            maps.append({"t":str(coords[0]),"s":str(coords[1]),**recovered})
    atlas = classify({"terms":[[*k,c] for k,c in terms.items()]})
    blocks = profile["elliptic_blocks"]
    unique = len(blocks)==1 and blocks[0]["geometric_multiplicity"]==1
    exhausted = bool(maps) and unique and atlas["geometric_involutions"]==0
    return {"status":"unique primitive elliptic quotient in all degrees" if exhausted else "unresolved",
            "prime":prime,"point_counts":profile["point_counts"],"weil":profile["weil"],
            "elliptic_blocks":blocks,"involutions":atlas["geometric_involutions"],
            "search":stats,"presentations":maps}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quartic",type=Path)
    parser.add_argument("--prime",type=int,default=5)
    args=parser.parse_args()
    print(json.dumps(discover(args.quartic,args.prime),indent=2))
