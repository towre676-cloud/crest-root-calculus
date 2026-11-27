"""Check the Eisenstein obstruction, all three Chabauty lifts and the uniform prime ray."""
from fractions import Fraction
from itertools import product
from eisenstein_quadratic import (
    CHABAUTY, RAY_MODULUS, MODULI, local_filter, norm, pari, positive_pair,
    ray_generator, representatives, sieve, solve_residual_pair, square_rational, squarefree,
)

CASES={19:(1,19),31:(3,31),43:(1,43),67:(1,67),79:(3,79),95:(1,57),
       103:(3,103),127:(3,127),129:(1,129),139:(1,139),151:(3,151),
       163:(1,163),182:(7,13),193:(1,193),199:(3,199)}
def rational_survivors(p):
    pairs=set(product((1,3,p,3*p), repeat=2))
    for prime,m in ((2,64),(3,81)):
        squares={i*i%m for i in range(m)}
        lines=[(j,1) for j in range(m)]+[(1,j) for j in range(0,m,prime)]
        pairs={(a,b) for a,b in pairs if any(
            a*(u*u+u*v+v*v)%m in squares and b*(u*u-u*v+v*v)%m in squares
            and p*a*b*(u*u-v*v)%m in squares for u,v in lines)}
    return pairs

def verify_ray(p):
    pi=ray_generator(p)
    assert pi and norm(pi)==p and p%RAY_MODULUS==43
    assert rational_survivors(p)=={(1,3),(3,1),(1,p),(p,1),(p,3*p),(3*p,p)}
    options=set(product(representatives(1,1),representatives(p,1,{p:pi})))
    counts=[len(options)]
    for prime,m in MODULI[:3]:
        options=local_filter(options,1,p,1,prime,m)
        counts.append(len(options))
    assert counts==[32,2,1,0]
    assert not options
    print('prime ray:',p,pi,counts)

def verify():
    E=pari.ellinit([0,3,0,9,0])
    assert list(E.ellrank()[:2])==[0,0]
    n,orders,gens=E.elltors()
    points=[pari([0])]
    for order,T in zip(orders,gens):
        points=[E.elladd(P,E.ellmul(T,j)) for P in points for j in range(int(order))]
    assert len(points)==int(n)==4
    fibres=[Fraction(str(P[0]/3)) for P in points if len(P)>1]
    assert {x for x in fibres if square_rational(x)}=={0,1}
    for d,(a,b) in CASES.items():
        k=squarefree(d*a*b)
        rem,counts=sieve(a,b,k)
        target=CHABAUTY.get((d,a,b))
        assert rem==([target[0]] if target else [])
        print(d,a,b,k,counts)
        xs=solve_residual_pair(d,a,b)
        expected={19:{Fraction(5,3)},95:{Fraction(8,7)},182:{Fraction(-3)}}.get(d,set())
        assert xs==expected
        assert solve_residual_pair(-d,a,b)=={1/x for x in expected}
        assert solve_residual_pair(d,b,a)=={-x for x in expected}
    for p in (43,1317109291,1317290731):
        verify_ray(p)
    print('All fifteen residual twists, reciprocal lifts, and prime-ray checks passed.')
if __name__=='__main__':verify()
