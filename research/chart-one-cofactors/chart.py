"""Eliminate the conic parameters of a marked transverse chart over QQ."""
from dataclasses import dataclass
from itertools import permutations
from math import comb
import sympy as sp
from sympy.polys.rings import ring


X = sp.symbols("X Y Z")
R, c, s = ring("c,s",sp.QQ)
C, S = sp.symbols("c s")
# Order of the four remaining generators (coefficient of t^i z^j).
POSITIONS = ((2,0),(2,1),(1,0),(1,1))


def adapted_frame(data, point):
    terms = data["terms"]
    if any(len(row)!=4 or any(type(n)!=int or n<0 for n in row[:3])
           or sum(row[:3])!=4 for row in terms):
        raise ValueError("Expected ternary quartic monomials.")
    F = sp.Poly.from_dict({},X,domain=sp.QQ)
    for i,j,k,a in terms:
        F += sp.Poly(sp.Rational(a)*X[0]**i*X[1]**j*X[2]**k,*X,domain=sp.QQ)
    pt = tuple(sp.Rational(a) for a in point)
    if len(pt)!=3 or not any(pt):
        raise ValueError("Expected three rational projective point coordinates.")
    evaluate = lambda poly: poly.as_expr().subs(dict(zip(X,pt)))
    jac = [evaluate(F.diff(x)) for x in X]
    if evaluate(F) or not any(jac):
        raise ValueError("The specified mark must be a smooth point of the curve.")
    # Singularities elsewhere would invalidate the elliptic interpretation.
    for n in range(3):
        aff = [x for j,x in enumerate(X) if j!=n]
        partials = [F.diff(x).as_expr().subs(X[n],1) for x in X]
        if list(sp.groebner(partials,*aff,domain=sp.QQ)) != [1]:
            raise ValueError("The projective source must be smooth.")
    i,j,n = next(f for f in permutations(range(3)) if pt[f[2]] and jac[f[1]])
    u,v,d = pt[i]/pt[n], pt[j]/pt[n], -jac[i]/jac[j]
    L,H,N = sp.symbols("L H N")
    T = sp.Poly(F.as_expr().xreplace({X[i]:L+u*N,X[j]:H+d*L+v*N,X[n]:N}),L,H,N)
    # Primitive integral coefficients decrease the height of every following identity.
    _, T = T.clear_denoms(convert=True)
    _, T = T.primitive()
    if T.LC()<0:
        T = -T
    return T, {"axes":[i,j,n],"u":str(u),"v":str(v),"d":str(d)}


def hasse_table(T, polynomial_ring):
    """[u^i v^j] T(s-u,-1,c+v), with entries in QQ[c,s]."""
    cc,ss = polynomial_ring.gens
    table = {(i,j):polynomial_ring.zero for i in range(5) for j in range(4)}
    for (i,j,k),a in T.rep.to_dict().items():
        if k>3:
            raise ArithmeticError("Mark normalization failed.")
        for h in range(i+1):
            for l in range(k+1):
                coeff = polynomial_ring.domain.convert(a)*(-1)**(j+h)*comb(i,h)*comb(k,l)
                table[h,l] += polynomial_ring.ground_new(coeff)*ss**(i-h)*cc**(k-l)
    return table


def square_generators(g, mul):
    """Denominator-cleared coefficients N2 and N1, including their z² terms."""
    h,a = g(6,0),g(5,0)
    h2 = mul(h,h)
    h3,h4 = mul(h,h2),mul(h2,h2)
    b0 = 4*mul(g(4,0),h)-mul(a,a)
    b1 = g(4,1)
    c0 = 8*mul(g(3,0),h2)-mul(a,b0)
    c1 = 2*mul(g(3,1),h)-mul(a,b1)
    u,v = 4*mul(b1,h),4*mul(c1,h)
    return {
        (2,0):64*mul(g(2,0),h3)-mul(b0,b0)-4*mul(a,c0),
        (2,1):64*mul(g(2,1),h3)-2*mul(b0,u)-4*mul(a,v),
        (2,2):64*mul(g(2,2),h3)-mul(u,u),
        (1,0):64*mul(g(1,0),h4)-mul(b0,c0),
        (1,1):64*mul(g(1,1),h4)-mul(b0,v)-mul(u,c0),
        (1,2):64*mul(g(1,2),h4)-mul(u,v),
    }


@dataclass
class Chart:
    frame: dict
    table: dict
    q: object
    h: object
    denominator: object
    numerator: object
    rows: list


def build(data, point):
    T, frame = adapted_frame(data,point)
    return build_from_frame(T,frame,R)


def build_from_frame(T, frame, polynomial_ring):
    a = hasse_table(T,polynomial_ring)
    h,u,v,w = (a[0,i] for i in range(4))
    q = (4*h*v-u*u).monic()
    if q.degree(0)!=4 or q.LC!=polynomial_ring.domain.one or not w.is_ground or not w:
        raise ArithmeticError("Expected a monic osculation quartic and nonzero cubic coefficient.")
    reduce = lambda p:p.rem([q])
    mul = lambda x,y:reduce(x*y)
    d = 12*h*h*w
    n = -(4*h*h*a[1,2]-2*h*u*a[1,1]+u*u*a[1,0])

    def powers(value):
        result = [polynomial_ring.one]
        for _ in range(6):
            result.append(mul(result[-1],value))
        return result

    nums,dens = powers(reduce(n)),powers(reduce(d))

    def g(i,l):
        # For b=n/d this is d^(6-i-2l) times [t^i z^l]G(t,z) modulo q.
        weight = 6-i-2*l
        parts = (a[j,l+weight-j]*comb(l+weight-j,l)*mul(nums[weight-j],dens[j])
                 for j in range(weight+1) if j<=4 and l+weight-j<=3)
        return reduce(sum(parts,polynomial_ring.zero))

    residuals = square_generators(g,mul)
    if residuals[2,2] or residuals[1,2]:
        raise ArithmeticError("The two osculation equations failed identically.")
    return Chart(frame,a,q,h,d,n,[residuals[k] for k in POSITIONS])
