"""Elliptic-factor obstruction at good primes for integral plane quartics.

Usage: python frobenius_obstruction.py quartic.json 2 3 5
Input "terms" lists [i,j,k,c] for c*x**i*y**j*z**k. Requires SymPy, NumPy.
"""
from itertools import product
from math import lcm
from pathlib import Path
import argparse
import json
import numpy as np
import sympy as sp

x,y,z,T,U = sp.symbols("x y z T U")


def smooth(terms,p):
    form = sum(c*x**i*y**j*z**k for i,j,k,c in terms)
    equations = [form]+[sp.diff(form, var) for var in (x,y,z)]
    for fixed in (z,y,x):
        others = [v for v in (x,y,z) if v != fixed]
        if list(sp.groebner([f.subs(fixed,1) for f in equations],*others,modulus=p)) != [1]:
            return False
    return True


def irreducibles(p, degree):
    """Coefficient lists in ascending order. Only degrees 1, 2 and 3 are used."""
    for coeff in product(range(p), repeat=degree):
        if degree == 1 or not any(sum(c*a**i for i,c in enumerate((*coeff,1)))%p == 0
                                  for a in range(p)):
            yield (*coeff,1)


class Field:
    """Polynomial-basis finite field, encoded in base p; operations on scalars."""
    def __init__(self,p,modulus):
        self.p=p
        self.m=tuple(modulus)
        self.n=len(modulus)-1
        self.q=p**self.n
        self.digits=[tuple((a//p**i)%p for i in range(self.n)) for a in range(self.q)]

    def encode(self,coeffs):
        return sum((c%self.p)*self.p**i for i,c in enumerate(coeffs))

    def add(self,a,b):
        return self.encode([u+v for u,v in zip(self.digits[a],self.digits[b])])

    def neg(self,a):
        return self.encode([-u for u in self.digits[a]])

    def mul(self,a,b):
        coeff=[0]*(2*self.n-1)
        for i,u in enumerate(self.digits[a]):
            for j,v in enumerate(self.digits[b]):
                coeff[i+j] += u*v
        for k in range(len(coeff)-1,self.n-1,-1):
            for j in range(self.n):
                coeff[k-self.n+j] -= coeff[k]*self.m[j]
        return self.encode(coeff[:self.n])

    def power(self,a,n):
        value=1
        while n:
            if n%2: value=self.mul(value,a)
            a=self.mul(a,a)
            n//=2
        return value


class EnumeratingField:
    """Vectorized addition in the polynomial basis, multiplication in a log table."""
    def __init__(self,p,degree):
        self.p=p
        self.degree=degree
        self.q=q=p**degree
        # A primitive modulus lets multiplication by its variable enumerate
        # the whole multiplicative group. This implements field arithmetic
        # independently of Field, which the closed-fibre counter uses.
        for modulus in irreducibles(p,degree):
            sequence=[]
            value=1
            while True:
                sequence.append(value)
                digits=[(value//p**i)%p for i in range(degree)]
                carry=digits[-1]
                value=sum(((old-carry*co)%p)*p**i
                          for i,(old,co) in enumerate(zip([0]+digits[:-1],modulus)))
                if value in (0,1):
                    break
            if value==1 and len(sequence)==q-1:
                break
        else:
            raise ArithmeticError("No primitive finite-field polynomial")
        self.antilog=np.array(sequence,dtype=np.int64)
        self.log=np.zeros(q,dtype=np.int64)
        for exponent,value in enumerate(sequence):
            self.log[value]=exponent

    def mul(self,a,b):
        return np.where((a==0)|(b==0),0,self.antilog[(self.log[a]+self.log[b])%(self.q-1)])

    def add(self,a,b):
        out=0
        for i in range(self.degree):
            q=self.p**i
            out+=((a//q%self.p+b//q%self.p)%self.p)*q
        return out


def enumerate_points(terms,p,n):
    field=EnumeratingField(p,n)
    coords=np.arange(field.q)
    powers=[np.ones(field.q,dtype=np.int64)]
    for _ in range(4): powers.append(field.mul(powers[-1],coords))
    # Pairs [x:y:1], then [x:1:0] and [1:0:0].
    def affine_chart(y_values,z_val):
        count=0
        for yy in y_values:
            co=[0]*5
            for i,j,k,c in terms:
                if z_val or k==0:
                    co[i]=field.add(co[i],field.mul(c%p,powers[j][yy]))
            values=np.zeros(field.q,dtype=np.int64)
            for i in range(5): values=field.add(values,field.mul(co[i],powers[i]))
            count+=int(np.count_nonzero(values==0))
        return count
    return (affine_chart(range(field.q),1)+affine_chart([1],0)
            +int(sum(c for i,j,k,c in terms if i==4)%p==0))


def trim(poly):
    while poly and not poly[-1]: poly.pop()
    return poly


def roots_in_extension(poly,field,power):
    """Number of distinct roots over F_{p^power} by gcd with Y^(p^power)-Y."""
    poly=trim(list(poly))
    if not poly: return field.p**power
    if len(poly)==1: return 0
    def rem(a,b):
        a=a.copy()
        inv=field.power(b[-1],field.q-2)
        while len(a)>=len(b):
            lead=field.neg(field.mul(a[-1],inv))
            offset=len(a)-len(b)
            for i,v in enumerate(b): a[offset+i]=field.add(a[offset+i],field.mul(lead,v))
            trim(a)
        return a
    def product_mod(a,b):
        if not a or not b: return []
        c=[0]*(len(a)+len(b)-1)
        for i,s in enumerate(a):
            for j,t in enumerate(b):
                c[i+j]=field.add(c[i+j],field.mul(s,t))
        return rem(trim(c),poly)
    n=field.p**power
    base=rem([0,1],poly)
    result=[1]
    while n:
        if n%2: result=product_mod(result,base)
        base=product_mod(base,base)
        n//=2
    result += [0]*max(0,2-len(result))
    result[1]=field.add(result[1],field.neg(1))
    a,b=poly,trim(result)
    while b: a,b=b,rem(a,b)
    return len(a)-1


def closed_point_count(terms,p,n):
    """Disjoint closed x-fibres. Every orbit contributes its degree times its count."""
    total=0
    for d in range(1,n+1):
        if n%d: continue
        for modulus in irreducibles(p,d):
            field=Field(p,modulus)
            xval=(-modulus[0])%p if d==1 else p
            vals=[0]*5
            for i,j,k,c in terms:
                vals[j]=field.add(vals[j],field.mul(c%p,field.power(xval,i)))
            total += d*roots_in_extension(vals,field,n)
    line=[sum(c for i,j,k,c in terms if k==0 and i==a)%p for a in range(5)]
    return total+roots_in_extension(line,Field(p,(0,1)),n)+int(line[4]==0)


def power_mod(base,exponent,modulus):
    result=base.one
    while exponent:
        if exponent%2:
            result=(result*base)%modulus
        base=(base*base)%modulus
        exponent//=2
    return result


def classify_polynomial(W,p):
    """Classify geometric elliptic isogeny factors of an actual Jacobian reduction."""
    W=sp.Poly(W,T,domain=sp.ZZ)
    ratio=sp.Poly(sp.resultant(W.as_expr(),W.as_expr().subs(T,U*T),T),U)
    orders=[]
    for factor,multiplicity in ratio.factor_list()[1]:
        if factor.is_cyclotomic:
            d=factor.degree()
            # phi(k)^2 >= k/2, prime power by prime power.
            order=next(k for k in range(1,2*d*d+1)
                       if sp.totient(k)==d and sp.cyclotomic_poly(k,U)==factor.monic().as_expr())
            orders.append(order)
    N=lcm(2,*orders)
    # Replace T^N by its remainder before elimination.
    lifted=power_mod(sp.Poly(T,T,domain=sp.ZZ),N,W)
    extended=sp.Poly(sp.resultant(W.as_expr(),U-lifted.as_expr(),T),U)
    factors=extended.factor_list()[1]
    elliptic=any(A.degree()==1 or (A.degree()==2 and int(A.nth(1))%p and A.nth(0)==p**N)
                 for A,m in factors)
    return {'torsion_orders':sorted(orders),'power':N,'powered_factors':
            [[list(map(int,A.all_coeffs())),m] for A,m in factors],
            'geometric_elliptic_factor_at_p':elliptic}


def analyze(terms,p,independent=True):
    if not sp.isprime(p): raise ValueError('The reduction characteristic must be prime')
    if not terms or any(len(t)!=4 or min(t[:3])<0 or sum(t[:3])!=4
                        or any(type(v)!=int for v in t) for t in terms):
        raise ValueError('Expected integral coefficient terms of a ternary quartic')
    if not smooth(terms,p): return {'prime':p,'smooth':False,'conclusion':'bad reduction'}
    counts=[enumerate_points(terms,p,n) for n in (1,2,3)]
    if independent and counts != [closed_point_count(terms,p,n) for n in (1,2,3)]:
        raise ArithmeticError('Point counters disagree')
    s1,s2,s3=[p**n+1-count for n,count in enumerate(counts,1)]
    if (s1*s1-s2)%2 or (s1**3-3*s1*s2+2*s3)%6:
        raise ArithmeticError('Newton identities failed integrality')
    a1=-s1
    a2=(s1*s1-s2)//2
    a3=-(s1**3-3*s1*s2+2*s3)//6
    coeff=[1,a1,a2,a3,p*a2,p*p*a1,p**3]
    info=classify_polynomial(sp.Poly.from_list(coeff,T),p)
    verdict='inconclusive in characteristic zero' if info['geometric_elliptic_factor_at_p'] else 'no geometric elliptic quotient in characteristic zero'
    return {'prime':p,'smooth':True,'point_counts':counts,'frobenius_coefficients':coeff,
            **info,'conclusion':verdict}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('quartic',type=Path)
    parser.add_argument('primes',type=int,nargs='+')
    args=parser.parse_args()
    terms=json.loads(args.quartic.read_text())['terms']
    print(json.dumps([analyze(terms,p) for p in args.primes],indent=2))
