"""Check absolute simplicity for the full congruence family using its common reduction."""
from itertools import product
from math import isqrt
import argparse
import flint
import sympy as sp
from cypari2 import Pari

pari=Pari()
x,t=sp.symbols('x t')
g=x**4+1
f=g*(1+x*g)
def check(c,msg):
    if not c: raise ArithmeticError(msg)

def finite_count(n):
    k=flint.fq_default_ctx(5,n)
    F=flint.fq_default_poly_ctx(k)([k(int(a)) for a in sp.Poly(f,x).all_coeffs()[::-1]])
    ans=1
    for a in product(range(5),repeat=n):
        r=F(k(list(a)))
        ans+=1 if r.is_zero() else 2*int(r.is_square())
    return ans

def classify(k):
    m=1+5*k
    r=isqrt(abs(m))
    return (['infinity']+([(0,str(r)),(0,str(-r))] if m>0 and r*r==m else
                           [(0,str(r)+'i'),(0,str(-r)+'i')] if m<0 and r*r==-m else []))

def verify():
    check(sp.Poly(f,x,modulus=5).is_sqf,'smooth reduction')
    counts=[finite_count(n) for n in range(1,5)]
    powers=[5**n+1-v for n,v in enumerate(counts,1)]
    coeff=[1]
    for n in range(1,5):
        numerator=-sum(coeff[n-i]*powers[i-1] for i in range(1,n+1))
        check(numerator%n==0,'Newton identity');coeff.append(numerator//n)
    P=sp.Poly.from_list(coeff+[5**i*coeff[4-i] for i in range(1,5)],x)
    check(P.all_coeffs()==[1,2,4,-2,6,-10,100,250,625] and P.is_irreducible,'Frobenius')
    R=sp.Poly(sp.resultant(P.as_expr(),P.as_expr().subs(x,t*x),x),t)
    A=R.exquo(sp.Poly((t-1)**8,t))
    groups=A.factor_list()[1]
    check(all(not a.is_cyclotomic for a,n in groups),'root-of-unity quotient')
    check(sorted((a.degree(),n) for a,n in groups)==[(8,1),(8,2),(16,2)],'factor degrees')
    for b in (1,-1):
        E=pari.ellinit([b,0]);check(list(E.ellrank()[:2])==[0,0],'rank over Q')
        K=pari.nfinit('w^2+1');G=pari.ellinit([b,0],K).elltors()
        check(int(G[0])=={1:4,-1:8}[b],'torsion over Q(i)')
    print('Finite-field counts:',counts,'Frobenius:',P.as_expr())
    print('Noncyclotomic factor degrees and multiplicities:',[(a.degree(),n) for a,n in groups])
    print('Parameters -1, 0, 1:',[classify(k) for k in (-1,0,1)])
if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__);cli.add_argument('k',nargs='*',type=int)
    args=cli.parse_args()
    if not args.k:verify()
    for k in args.k:print(k,classify(k))
