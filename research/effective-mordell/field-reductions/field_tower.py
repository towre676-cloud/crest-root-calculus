"""Finite checks for the uniform rank and point proof in the accompanying note."""
from itertools import combinations
from math import prod
from cypari2 import Pari

p=Pari()

def check(c,msg):
    if not c: raise ArithmeticError(msg)

def verify():
    primes=(5,17,29,41,53,89,101)
    count=0
    for k in range(4):
        for subset in combinations(primes,k):
            for sign in (-3,-1,1,3):
                d=sign*prod(subset)
                E=p.ellinit([0,-d**3])
                check(list(E.ellrank()[:2])==[0,0],'sampled twist has other rank')
                check(int(E.elltors()[0])==(6 if d==-1 else 2),'rational torsion')
                count+=1
    x,d=p('x'),p('d')
    check((x+d)**3-d**3==x*(x**2+3*d*x+3*d**2),'2-isogeny model')
    check((3*d)**2-4*(3*d**2)==-3*d*d,'first discriminant')
    check((-6*d)**2-4*(-3*d*d)==48*d*d,'dual discriminant')
    X=p('x')
    K=p.nfinit('w^4-w^2+1')
    check(int(p.ellinit([0,-1],K).elltors()[0])==12,'torsion over the twelfth cyclotomic field')
    check(p.polresultant(X**2-3,p('y')**2-9-6*X,X)==p('y')**4-18*p('y')**2-27,
          'mixed-signature 4-torsion ordinate')
    for b in [-1,1]:
        check(list(p.ellinit([0,b]).ellrank()[:2])==[0,0],'base ranks')
        check(list(p.ellinit([0,b*7**3]).ellrank()[:2])==[1,1],'positive ranks at seven')
    for q in (11,23):
        for d0 in range(1,q):
            check(int(p.ellinit([0,-d0**3],q).ellcard())==q+1,'torsion reduction')
    print(count,'signed twist checks passed; both standard quotients have rank one over Q(sqrt(7)).')
if __name__=='__main__':verify()
