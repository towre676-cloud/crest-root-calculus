"""Frobenius orbit degrees and traces after removing all possible root-of-unity twists."""
from functools import lru_cache

import numpy as np
import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_gcd, gf_quo, gf_sqf_list


@lru_cache(None)
def twist_power(dimension):
    bound = (2*dimension)**2
    M = 1
    for ell in sp.primerange(2,bound+2):
        q = ell
        while q*(ell-1) <= bound:
            q *= ell
        M *= q
    return M


class ResidueRing:
    """F_p[x]/h for a monic h; elements have ascending, fixed-length coefficients."""
    def __init__(self, h, p):
        self.p, self.n = p, len(h)-1
        self.h = list(map(int,h))
        n = self.n
        self.reductions = np.empty((max(1,n-1),n),dtype=np.int64)
        self.reductions[0] = (-np.asarray(h[:0:-1],dtype=np.int64))%p
        for j in range(1,n-1):
            prev = self.reductions[j-1]
            row = self.reductions[j]
            row[0],row[1:] = 0,prev[:-1]
            row[:] = (row+prev[-1]*self.reductions[0])%p
        self.one = np.zeros(n,dtype=np.int64)
        self.one[0] = 1
        self.x = np.zeros(n,dtype=np.int64)
        if n==1:
            self.x[0] = -int(h[1])%p
        else:
            self.x[1] = 1

    def mul(self,a,b):
        c = np.convolve(a,b)%self.p
        if self.n==1:
            return c
        return (c[:self.n]+c[self.n:]@self.reductions)%self.p

    def power(self,a,e):
        val = self.one
        while e:
            if e&1:
                val = self.mul(val,a)
            e//=2
            if e:
                a = self.mul(a,a)
        return val

    def frobenius(self):
        image = self.power(self.x,self.p)
        phi = np.empty((self.n,self.n),dtype=np.int64)
        phi[0] = self.one
        for i in range(1,self.n):
            phi[i] = self.mul(phi[i-1],image)
        return phi

    def trace_weights(self):
        # Newton sums of x^0,...,x^(n-1), with the multiplicities of h.
        values = [self.n]
        for i in range(1,self.n):
            total = i*self.h[i]+sum(self.h[j]*values[i-j] for j in range(1,i))
            values.append(-total%self.p)
        return np.asarray(values,dtype=np.int64)


def descending(element,p):
    data=list(map(int,(element%p)[::-1]))
    first=next((i for i,v in enumerate(data) if v),len(data))
    return data[first:]


def degree_blocks(squarefree,p):
    ring=ResidueRing(squarefree,p)
    phi=ring.frobenius()
    current=ring.x
    remaining=squarefree
    degree=1
    while 2*degree <= len(remaining)-1:
        current=current@phi%p
        part=gf_gcd(remaining,descending(current-ring.x,p),p,ZZ)
        if len(part)>1:
            yield degree,part
            remaining=gf_quo(remaining,part,p,ZZ)
        degree+=1
    if len(remaining)>1:
        yield len(remaining)-1,remaining


def observations(coefficients,p,dimension):
    """Compute permitted proper dimensions and comparable traces. Inputs are det(1-TF)."""
    if not sp.isprime(p) or not (dimension<p<=10000) or not 1<=dimension<=512:
        raise ValueError("Use a prime above the dimension and within the exact arithmetic bounds")
    if len(coefficients)!=dimension+1 or coefficients[0]!=1 or any(
            type(c) is not int or not 0<=c<p for c in coefficients):
        raise ValueError("Expected a normalized local numerator through the factor dimension")
    last=next(i for i in range(dimension,-1,-1) if coefficients[i])
    traces=np.zeros(3,dtype=np.int64)
    possible=1
    M=twist_power(dimension)
    for squarefree,multiplicity in gf_sqf_list(list(map(int,coefficients[:last+1])),p,ZZ)[1]:
        for n,group in degree_blocks(squarefree,p):
            ring=ResidueRing(group,p)
            y=ring.power(ring.x,M%(p**n-1))
            # Count roots of powered degree exactly s via fixed points of Frobenius^s.
            divs=list(map(int,sp.divisors(n)))
            counts={}
            fy=y
            if len(divs)>1:
                phi=ring.frobenius()
                old=0
                for s in divs[:-1]:
                    for _ in range(s-old):
                        fy=fy@phi%p
                    fixed=len(gf_gcd(group,descending(fy-y,p),p,ZZ))-1
                    counts[s]=fixed-sum(v for k,v in counts.items() if s%k==0)
                    old=s
            counts[n]=ring.n-sum(counts.values())
            for s,total in counts.items():
                if total<0 or total%s:
                    raise ArithmeticError("Invalid orbit multiplicity")
                for _ in range(multiplicity*total//s):
                    possible |= possible<<s
            weights=ring.trace_weights()
            yj=y
            for j in range(len(traces)):
                traces[j] = (traces[j]+multiplicity*int(yj@weights))%p
                if j+1<len(traces):
                    yj=ring.mul(yj,y)
    dims=None if last!=dimension else [i for i in range(1,dimension//2+1) if (possible>>i)&1]
    return dims, [last,*map(int,traces)]
