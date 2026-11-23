"""Local witnesses, descent and distinct moduli on an infinite Richelot family."""
from itertools import permutations,product
from math import comb
import sympy as sp

from frobenius import X, compare
from partitions import pairings


def curves(a):
    f=[a**3+2*a**2+a+1,-3*a-1,3*a*a+3*a+1,-2,3*a+1,0,1]
    g=[1-a**3-a*a,-2*(3*a+2),3*a*a+6*a+4,-2,-3*a-5,0,1]
    return [sp.Poly.from_list(c[::-1],X,domain=sp.QQ) for c in (f,g)]


def invariant(f):
    # The order-six transvectant up to scalar; total coefficient degree two.
    I=sum((-1)**i*f.nth(i)*f.nth(6-i)/comb(6,i) for i in range(7))
    return I**5/f.discriminant()


def packet():
    a,t,z=sp.symbols('a t z')
    quart=X**2+t*X+t*t+a
    m=t**3-t-1
    trace=sp.Matrix([[3,0,2],[0,2,3],[2,3,2]])
    basis=sp.Matrix([[1,0,a],[0,1,0],[0,0,1]])
    dual=(basis.T*trace).inv()*sp.Matrix([1,-2*z,z*z])
    ell=sp.Matrix([1,t,t*t]).dot(dual)
    f=sp.Poly(sp.resultant(m,quart,t),X)
    g=sp.Poly(23*sp.resultant(m,ell,t),z)
    F,G=curves(1)
    # The only possible collisions between distinct branch quadratics occur at a=-1.
    # Distinct cubic roots t,u satisfy t^2+t*u+u^2=1. This identity works integrally.
    u=sp.symbols('u')
    relation=t*t+t*u+u*u-1
    resultant=sp.resultant(quart,quart.subs(t,u),X)
    collision=sp.rem(resultant-(t-u)**2*(a+1),relation,u)
    disc_factor=64*a**3+96*a**2+36*a+27
    if (sp.expand(collision) or sp.expand(f.discriminant()+529*(a+1)**6*disc_factor)
        or sp.expand(f.as_expr().subs(a,-1)-(X**3-X-1)**2)):
        raise ArithmeticError('The symbolic branch-collision or discriminant certificate changed')
    # Direct fibre tests are separate from the symbolic argument covering all primes.
    reductions=singular_reductions=0
    for p in sp.primerange(3,98):
        if p==23:
            continue
        for residue in range(p):
            h=sp.Poly(curves(residue)[0].as_expr(),X,modulus=p)
            repeated=h.gcd(h.diff())
            if repeated.gcd(h.diff().diff()).degree()>0:
                raise ArithmeticError('Triple branch root in a claimed semistable fibre')
            singular_reductions+=int(repeated.degree()>0)
            reductions+=1
    if (f.eval(0)-(a**3+2*a**2+a+1)).expand() or trace.det()!=-23:
        raise ArithmeticError('Incorrect integral family or trace pairing')
    for j in range(7):
        ff=[a**3+2*a**2+a+1,-3*a-1,3*a*a+3*a+1,-2,3*a+1,0,1]
        gg=[1-a**3-a*a,-2*(3*a+2),3*a*a+6*a+4,-2,-3*a-5,0,1]
        if sp.expand(f.nth(j)-ff[j]) or sp.expand(g.nth(j)-gg[j]):
            raise ArithmeticError('Symbolic family norm changed')
    for h in (F,G):
        if not sp.Poly(h.as_expr(),X,modulus=3).is_irreducible:
            raise ArithmeticError('Missing common full-cycle Frobenius')
    fac=sp.Poly(F.as_expr(),X,modulus=59).factor_list()[1]
    if sorted(k.degree() for k,e in fac for _ in range(e))!=[1,1,1,1,2]:
        raise ArithmeticError('Missing one-pair flip')
    if invariant(F)==invariant(curves(87262)[0]):
        raise ArithmeticError('The moduli invariant does not separate the examples')
    # Enumerate C2 wr S3 on the six branches independently of the degree-15 resolvent.
    perms=[]
    for block in permutations(range(3)):
        for flips in product((0,1),repeat=3):
            perms.append([2*block[i//2]+(i%2 ^ flips[i//2]) for i in range(6)])
    canon=lambda pairs:tuple(sorted(tuple(sorted(q)) for q in pairs))
    unvisited={canon(m) for m in pairings(list(range(6)))}
    sizes=[]
    while unvisited:
        m=min(unvisited)
        orbit={canon((p[i],p[j]) for i,j in m) for p in perms}
        sizes.append(len(orbit));unvisited-=orbit
    return {'endomorphisms_Z_modulus':17*29,'unique_rational_two_power_class_modulus':3*17*29,
            'full_signed_permutation_group_modulus':17*29*59,
            'simultaneous_family_modulus':3*17*29*59, 'geometric_kernel_orbits':sorted(sizes),
            'source_discriminant_descending_cubic_factor':[64,96,36,27],
            'semistability_exceptions':[2,23],'checked_reduction_fibres':reductions,
            'checked_singular_nodal_fibres':singular_reductions,
            'field_witnesses':compare(F,G,[17,29])['source_endomorphisms_Z_witnesses']}
