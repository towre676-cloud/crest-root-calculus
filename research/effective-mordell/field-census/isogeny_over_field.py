"""Finite 2-isogeny image bounds and point fibres over class-number-one fields."""
from argparse import ArgumentParser
from functools import cached_property
from itertools import product
import json
from math import prod

from cypari2 import Pari

pari=Pari()


class LocalRing:
    """Arithmetic in the integral basis modulo a power of one prime ideal."""

    def __init__(self, field, ideal, basis_products):
        hnf=pari.idealhnf(field,ideal)
        self.cols=[tuple(map(int,col)) for col in hnf]
        self.degree=len(self.cols)
        self.diagonal=[col[i] for i,col in enumerate(self.cols)]
        self.basis_products=basis_products
        self.field=field

    def reduce(self, coordinates):
        result=list(coordinates)
        for j in reversed(range(self.degree)):
            col=self.cols[j]
            quotient=result[j]//col[j]
            for i in range(j+1):
                result[i]-=quotient*col[i]
        return tuple(result)

    def element(self, algebraic):
        coords=pari.nfalgtobasis(self.field,algebraic)
        if any(c.denominator()!=1 for c in coords):
            raise ArithmeticError('A coefficient of the integral local model was nonintegral.')
        return self.reduce(tuple(map(int,coords)))

    def mul(self,a,b):
        out=[0]*self.degree
        for i,ai in enumerate(a):
            if ai==0:continue
            for j,bj in enumerate(b):
                if bj==0:continue
                for t,c in enumerate(self.basis_products[i][j]):
                    out[t]+=ai*bj*c
        return self.reduce(out)

    @cached_property
    def squares(self):
        elements=product(*(range(n) for n in self.diagonal))
        return {self.mul(t,t) for t in elements}

    @cached_property
    def powers(self):
        return [(t2,self.mul(t2,t2)) for t2 in self.squares]

    def possible(self,e,a,b_over_e):
        # Both primitive projective charts. Unit scaling changes the quartic by a fourth power.
        e,a,last=map(self.element,(e,a,b_over_e))
        for t2,t4 in self.powers:
            mid=self.mul(a,t2)
            for first,constant in ((e,last),(last,e)):
                top=self.mul(first,t4)
                value=self.reduce(tuple(x+y+z for x,y,z in zip(top,mid,constant)))
                if value in self.squares:return True
        return False


def local_rings(field,limit):
    basis=list(field.nf_get_zk())
    table=[[tuple(map(int,pari.nfalgtobasis(field,a*b))) for b in basis] for a in basis]
    for prime in (3,2):
        for prime_ideal in pari.idealprimedec(field,prime):
            norm=int(pari.idealnorm(field,prime_ideal))
            ideal=prime_ideal
            exponent=1
            while norm**exponent<=limit:
                yield LocalRing(field,ideal,table)
                exponent+=1
                ideal=pari.idealmul(field,ideal,prime_ideal)


def isogeny_class_reps(field):
    # The unit, fundamental units, and generators of primes over 3 are independent mod squares.
    gens=[field.bnf_get_tu()[1],*list(field.bnf_get_fu())]
    for pr in pari.idealprimedec(field,3):
        exponents,vector=pari.bnfisprincipal(field,pr)
        if any(exponents):raise ArithmeticError('Prime ideal unexpectedly nonprincipal.')
        gen=pari.nfbasistoalg(field,vector)
        if pari.idealhnf(field,gen)!=pari.idealhnf(field,pr):
            raise ArithmeticError('Prime generator verification failed.')
        gens.append(gen)
    reps=[pari.Mod(1,field.nf_get_pol())]
    for gen in gens:
        reps += [r*gen for r in reps]
    return reps


def lift_torsion(field,sign):
    E=pari.ellinit([0,sign],field)
    order,orders,generators=E.elltors()
    points=[pari([0])]
    for n,T in zip(orders,generators):
        points=[E.elladd(P,E.ellmul(T,j)) for P in points for j in range(int(n))]
    if len({str(P) for P in points})!=int(order):raise ArithmeticError('Torsion enumeration.')
    X=pari('x')
    lifts=set()
    def add(x,y):
        if y**2!=x**6-1:raise ArithmeticError('Wrong source fibre.')
        lifts.add((str(pari.lift(x)),str(pari.lift(y))))
    for P in points:
        if len(P)==1:continue
        if sign==1 and P[0]==0:continue
        square=P[0] if sign==-1 else -1/P[0]
        for x in pari.nfroots(field,X**2-square):
            y=P[1] if sign==-1 else P[1]*x**3
            add(x,y)
    # The reciprocal map sends x=0 to infinity; restore these points independently.
    for y in pari.nfroots(field,X**2+1):add(pari(0),y)
    return sorted(lifts)


def classify(polynomial,max_residues=256,max_classes=256):
    if not isinstance(polynomial,list) or len(polynomial)<2 or any(type(c)!=int for c in polynomial) or polynomial[-1]!=1:
        raise ValueError('Supply ascending integer coefficients of a monic polynomial.')
    f=pari.Polrev(polynomial,pari('w'))
    if not pari.polisirreducible(f):raise ValueError('The defining polynomial must be irreducible.')
    field=pari.bnfinit(f,1)
    ans=dict(field=str(f),status='partial')
    if int(field.bnf_get_no())!=1:
        return dict(ans,reason='The current squareclass generator routine requires class number one.')
    if pari.bnfcertify(field)!=1:raise ArithmeticError('Unit and ideal class certificate failed.')
    count=2**(1+len(field.bnf_get_fu())+len(pari.idealprimedec(field,3)))
    if count>max_classes:return dict(ans,reason='Squareclass enumeration exceeded the chosen budget.')
    reps=isogeny_class_reps(field)
    positive=[r for r in reps if all(int(v)>0 for v in pari.nfeltsign(field,r))]
    tests={sign:[positive.copy(),reps.copy()] for sign in (-1,1)}
    for ring in local_rings(field,max_residues):
        for sign,images in tests.items():
            # Translate the rational root -sign to zero.
            a=-3*sign
            for j,(aj,bj) in enumerate(((a,3),(-2*a,-3))):
                images[j]=[e for e in images[j] if ring.possible(e,aj,bj/e)]
            bounds=list(map(len,images))
            if prod(bounds)<4:raise ArithmeticError('Local filters contradict the isogeny degree.')
            if prod(bounds)==4:
                return dict(ans,status='complete',sign=sign,isogeny_image_upper_bounds=bounds,
                            verified_affine_points=lift_torsion(field,sign),rational_infinities=2)
    return dict(ans,isogeny_image_upper_bounds={str(d):list(map(len,v)) for d,v in tests.items()},
                reason='Local image bounds leave positive rank possible.')

if __name__=='__main__':
    cli=ArgumentParser(description=__doc__)
    cli.add_argument('polynomial',help='Ascending integer coefficients as a JSON list.')
    cli.add_argument('--max-residues',type=int,default=256)
    cli.add_argument('--max-classes',type=int,default=256)
    args=cli.parse_args()
    print(json.dumps(classify(json.loads(args.polynomial),args.max_residues,args.max_classes),indent=2))
