"""Replay the blind recovery, the symbolic family and its rational isogeny sieve."""
import json
import sympy as sp

from branch_fibres import HERE, X, Y, Z, discover, read_terms


def family(r,t,s):
    edge = (4*t-5)*(t-2)
    a1, a3 = 3*r**2, -1-2*r**3
    A = a1*t*edge
    B = edge*(t**3*(4*t-5)+a3*(t-2))
    b = (t-1)*s-r
    F = sp.cancel((b**3-A*b-B)/(t-1)**2)
    return F, a1, a3


def main():
    path = HERE / "quartic.json"
    result = discover(path,5)
    if result != json.loads((HERE/"receipt.json").read_text()):
        raise ArithmeticError("Equation-only reconstruction differs from receipt")
    if result["status"] != "unique primitive elliptic quotient in all degrees":
        raise ArithmeticError("The global upper and lower bounds must meet")

    single, = result["presentations"]
    tt, ss = sp.symbols("t s")
    aff_line = sp.sympify(single["canonical_line_in_pencil_coordinates"]).subs(
        {tt:sp.sympify(single["t"]), ss:sp.sympify(single["s"])}, simultaneous=True)
    homogeneous_line = sp.Poly(aff_line,X,Y).homogenize(Z).as_expr()
    if sp.expand(homogeneous_line-12*(-2*X-2*Y+3*Z)) != 0:
        raise ArithmeticError("The solution of the unrestricted positive chart differs")
    r,t,s = sp.symbols("r t s")
    F,a1,a3 = family(r,t,s)
    if sp.Poly(F,t,s).total_degree()!=4 or sp.Poly(F,r,t,s).get_domain()!=sp.ZZ:
        raise ArithmeticError("The family is not integral of plane degree four")
    g = t**3*(4*t-5)/(t-2)
    x = t*((t-1)*s-r)/(t-2)
    if sp.cancel(x**3-a1*g*x-g*(g+a3)-t**3*(t-1)**2*F/(t-2)**3):
        raise ArithmeticError("Global family map failed its elliptic equation")
    base = sp.Poly(F.subs(r,1),t,s).homogenize(Z).as_expr()
    transformed = base.subs({t:X+Y+Z,s:X-Y+2*Z}, simultaneous=True).expand()
    if sp.Poly(transformed,X,Y,Z).as_dict()!=read_terms(path):
        raise ArithmeticError("Mixed-coordinate input and integral family disagree")
    if sp.factor(g.diff(t)-3*t**2*(t-1)*(4*t-10)/(t-2)**2)!=0:
        raise ArithmeticError("Unexpected critical points on the rational degree-four cover")
    if g.subs(t,sp.Rational(5,2)) != sp.Rational(625,4):
        raise ArithmeticError("Unexpected final branch value")
    z=sp.Symbol("z")
    x0 = sp.Symbol("x0")
    cubic = x0**3-a1*z*x0-z*(z+a3)
    if sp.expand(sp.discriminant(cubic,x0)+27*z**2*(z-1)*(z-a3**2)):
        raise ArithmeticError("Unexpected branch locus for the elliptic projection")

    candidates = {2,3,5,7,11,13,17,19,37,43,67,163}
    traces = {}
    for p in (5,7,13,17):
        trace = -sum(int(sp.legendre_symbol((4*i**3+(3*i-3)**2)%p,p)) for i in range(p))
        traces[p] = trace
        candidates = {ell for ell in candidates if ell==p or
                      any((v*v-trace*v+p)%ell==0 for v in range(ell))}
    if candidates!={3} or traces!={5:-3,7:-4,13:-1,17:-3}:
        raise ArithmeticError("Incomplete exclusion of other prime isogenies")
    if any((v*v-traces[7]*v+7)%9==0 for v in range(9)) or (8-traces[7])%9==0:
        raise ArithmeticError("The 9-isogeny and second 3-subgroup controls failed")
    print("Verified the unique geometric quotient, the infinite family and the rational degree spectrum.")


if __name__ == "__main__":
    main()
