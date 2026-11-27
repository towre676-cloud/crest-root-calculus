"""Compute every point of y²=x⁶−1 over Q(zeta_24), with rank two on both targets."""
from fractions import Fraction
from itertools import product
import json
from pathlib import Path

from cypari2 import Pari

p = Pari()
x,w = p('x'),p('w')
f = w**8-w**4+1
nf=p.nfinit(f)
z=p.Mod(w,f)
i=z**6
sqrt2=z**3+z**-3
sqrt3=z**2+z**-2
E=p.ellinit([0,-1],nf)
origin=p([0])

def numeric_bound_witness(d,coordinates):
    """Certify the full rational generator and a positive lower canonical height."""
    F=p.ellinit([0,-d**3])
    ranks=F.ellrank()
    if list(ranks[:2])!=[1,1] or int(F.elltors()[0])!=2:
        raise ArithmeticError('Unexpected rational rank or torsion.')
    multiples={str(F.elladd(T,F.ellmul(coordinates,n))) for T in (origin,p([d,0])) for n in range(-4,5)}
    small={Fraction(a,b) for a in range(-10,11) for b in range(1,11)}
    for q in small:
        X=p(str(q))*d
        for Y in F.ellordinate(X):
            if str(p([X,Y])) not in multiples:
                raise ArithmeticError('Generator saturation failed on the height window.')
    scaled=Fraction(str(F.ellmul(coordinates,4)[0]/d))
    height=max(abs(scaled.numerator),scaled.denominator)
    # h(4G)-A > 4 means H(G)>1/4, using exp(4)<55 and 132^(1/3)<6.
    if height<=330:raise ArithmeticError('Insufficient canonical height lower bound.')


def halves(P):
    if len(P)==1:raise ArithmeticError('Attempted to halve a zero free direction.')
    U=p('x');b=-1
    g=U**4-8*b*U-4*P[0]*(U**3+b)
    for X in p.nfroots(nf,g):
        for Y in p.nfroots(nf,U**2-X**3-b):
            R=p([X,Y])
            if E.ellmul(R,2)==P:return R
    return None


def all_torsion():
    order,orders,gens=E.elltors()
    pts=[origin]
    for n,T in zip(orders,gens):
        pts=[E.elladd(A,E.ellmul(T,j)) for A in pts for j in range(int(n))]
    if int(order)!=12 or len(set(map(str,pts)))!=12:raise ArithmeticError('Wrong cyclotomic torsion.')
    return pts


def saturate(B,T):
    # Trace projections over the elementary abelian group of order eight place the full free lattice in (1/8) <B>.
    coords=[(Fraction(1),Fraction(0)),(Fraction(0),Fraction(1))]
    index=1
    while True:
        updated=False
        for a,b in ((1,0),(0,1),(1,1)):
            P=E.elladd(E.ellmul(B[0],a),E.ellmul(B[1],b))
            for torsion in T:
                R=halves(E.elladd(P,torsion))
                if R is None:continue
                j=0 if a else 1
                coords[j]=tuple((a*u+b*v)/2 for u,v in zip(*coords))
                B[j]=R;index*=2;updated=True
                if index>64:raise ArithmeticError('Trace-index bound exceeded.')
                break
            if updated:break
        if not updated:return B,coords,index


def classify():
    numeric_bound_witness(-2,[1,3])
    numeric_bound_witness(6,[10,28])
    if not all(list(p.ellinit([0,-d**3]).ellrank()[:2])==[0,0] for d in (1,-1,2,3,-3,-6)):
        raise ArithmeticError('Rank-zero character check failed.')
    original=[p([-p(1)/2,3/(-2*i*sqrt2)]),p([p(5)/3,28/(6*sqrt2*sqrt3)])]
    T=all_torsion()
    B,C,index=saturate(original.copy(),T)
    det=C[0][0]*C[1][1]-C[0][1]*C[1][0]
    if abs(det)!=Fraction(1,index):raise ArithmeticError('Lattice determinant.')
    good=set()
    tested=0
    def save(a,b):
        if b*b!=a**6-1:raise ArithmeticError('Point failed source equation.')
        good.add((str(p.lift(a)),str(p.lift(b))))
    # Each orthogonal coefficient of the original basis has absolute value < 6.
    # Invert the change-of-basis rectangle to bound saturated lattice coordinates.
    from math import ceil
    limits=[ceil(6*(abs(C[1][1])+abs(C[1][0]))/abs(det)),
            ceil(6*(abs(C[0][0])+abs(C[0][1]))/abs(det))]
    for a,b in product(*(range(-n,n+1) for n in limits)):
        uv=[a*u+b*v for u,v in zip(*C)]
        if any(abs(t)>=6 for t in uv):continue
        Q=E.elladd(E.ellmul(B[0],a),E.ellmul(B[1],b))
        for t in T:
            P=E.elladd(Q,t)
            if len(P)==1:continue
            tested+=1
            for r in p.nfroots(nf,x*x-P[0]):
                save(r,P[1])
    expected=set()
    def record(a,b):
        expected.add((str(p.lift(a*z**0)),str(p.lift(b*z**0))))
    for h in range(6):
        unit=z**(4*h)
        record(unit,0)
        for sig in (-1,1):
            record(z**2*unit,sig*i*sqrt2)
            record(i*sqrt2*unit,sig*3*i)
            record(i*sqrt2*unit/2,sig*3*i*sqrt2/4)
    for sig in (-1,1):record(0,sig*i)
    if len(good)!=44 or good!=expected:
        raise ArithmeticError('Complete affine enumeration disagrees with the cyclotomic formulas.')
    return dict(field=str(f),status='complete',projective_points=46,infinities=2,index=index,basis_change=[[str(v) for v in row] for row in C],
                finite_elliptic_points_tested=tested,affine_points=sorted(good))

if __name__=='__main__':
    print(json.dumps(classify(),indent=2))
