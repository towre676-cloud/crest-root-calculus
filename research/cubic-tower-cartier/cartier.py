"""Cartier matrices for the full signed support on a cubic-tower stage."""
from itertools import product
from math import comb

import numpy as np
import sympy as sp


def trim(a, p):
    a = np.asarray(a,dtype=np.int64) % p
    nz = np.flatnonzero(a)
    return a[:nz[-1]+1] if len(nz) else np.array([0],dtype=np.int64)


def add(a,b,p):
    n = max(len(a),len(b))
    c = np.zeros(n,dtype=np.int64)
    c[:len(a)] += a
    c[:len(b)] += b
    return trim(c,p)


def mul(a,b,p):
    return trim(np.convolve(a,b),p)


def kernel_basis(M,p):
    M=M.copy() % p
    row=0
    pivots=[]
    for col in range(M.shape[1]):
        candidates=np.flatnonzero(M[row:,col])
        if not len(candidates):
            continue
        j=row+int(candidates[0])
        M[[row,j]]=M[[j,row]]
        M[row]=M[row]*pow(int(M[row,col]),-1,p)%p
        pivot=M[row].copy()
        M=(M-M[:,col:col+1]*pivot)%p
        M[row]=pivot
        pivots.append(col)
        row+=1
    free=[j for j in range(M.shape[1]) if j not in pivots]
    basis=np.zeros((M.shape[1],len(free)),dtype=np.int64)
    basis[free,:]=np.eye(len(free),dtype=np.int64)
    if pivots:
        basis[pivots,:]=(-M[:row,free])%p
    return free,basis


class Cartier:
    """Integer-polynomial work is reduced modulo p at every multiplication."""
    def __init__(self,p,k):
        if not sp.isprime(p) or not (5<=p<=1000) or not (1<=k<=5):
            raise ValueError('Use a prime from 5 to 997, and one through five indices')
        self.p,self.k=p,k
        self.f=np.array([1,1,1,-4],dtype=np.int64)
        self.fpower=[np.array([1],dtype=np.int64)]
        for _ in range(max(k, (p-1)//2)):
            self.fpower.append(mul(self.fpower[-1],self.f,p))
        self._cubic_cache={}
        self.presentations=[self.basis(sign) for sign in (1,-1)]

    def basis(self,sign):
        p,k=self.p,self.k
        labels=[]
        for e in range(1<<k):
            for b in (0,1):
                if (-1)**(e.bit_count()+b+1)!=sign:
                    continue
                for a in range((2*k-e.bit_count()-3*b)//2+1):
                    labels.append((e,b,a))
        constraints={}
        for col,(e,b,a) in enumerate(labels):
            j=e
            while True:
                n=(e^j).bit_count()
                even,odd=divmod(n+b,2)
                for l,coef in enumerate(self.fpower[even]):
                    power=a-n+l
                    if 2*power+3*odd>k:
                        key=(j,odd,power)
                        if key not in constraints:
                            constraints[key]=np.zeros(len(labels),dtype=np.int64)
                        constraints[key][col]+=int(coef)*pow(pow(2,n,p),-1,p)
                if j==0:
                    break
                j=(j-1)&e
        M=np.array(list(constraints.values()),dtype=np.int64).reshape(-1,len(labels))
        free,basis=kernel_basis(M,p)
        if basis.shape[1]!=k*(1<<(k-1)) or np.any(M@basis%p):
            raise ArithmeticError('Incorrect signed holomorphic subspace')
        return labels,free,basis

    def cubic_terms(self,index,e,h):
        key=(index,e,h)
        if key in self._cubic_cache:
            return self._cubic_cache[key]
        p=self.p
        v=(h+1)*p-1-e
        out=[np.array([0],dtype=np.int64) for _ in range(2)]
        for a in range(v//3+1):
            b=v-3*a
            c=p-1-a-b
            if c<0:
                continue
            scale=comb(p-1,a)*comb(p-1-a,b)*pow(3,b,p)*(-1)**c % p
            linear=np.array([comb(b,j)*pow(index,b-j,p)%p for j in range(b+1)],dtype=np.int64)
            value=mul(linear,self.fpower[c//2],p)
            out[c%2]=add(out[c%2],value*scale,p)
        self._cubic_cache[key]=out
        return out

    def product(self,a,b):
        p=self.p
        return [add(mul(a[0],b[0],p),mul(self.f,mul(a[1],b[1],p),p),p),
                add(mul(a[1],b[0],p),mul(a[0],b[1],p),p)]

    def matrices(self,indices):
        k,p=self.k,self.p
        if len(indices)!=k:
            raise ValueError('Wrong number of indices')
        outputs=[]
        for labels,free,basis in self.presentations:
            C=np.zeros((len(labels),len(labels)),dtype=np.int64)
            for h,e in product(range(1<<k),repeat=2):
                T=[np.array([1]),np.array([0])]
                for i,n in enumerate(indices):
                    T=self.product(T,self.cubic_terms(n,(e>>i)&1,(h>>i)&1))
                for b in (0,1):
                    A=T if not b else [mul(T[1],self.f,p),T[0]]
                    B=[mul(A[0],self.fpower[(p-1)//2],p),A[1]]
                    columns=[(j,a) for j,(ee,bb,a) in enumerate(labels) if (ee,bb)==(e,b)]
                    for out,(hh,bit,power) in enumerate(labels):
                        if hh!=h:
                            continue
                        for col,a in columns:
                            idx=p*(power+1)-1-a
                            if 0<=idx<len(B[bit]):
                                C[out,col]=B[bit][idx]
            image=C@basis % p
            small=image[free,:]
            if np.any((basis@small-image)%p):
                raise ArithmeticError('Cartier failed to preserve the holomorphic subspace')
            outputs.append(small)
        return outputs


def numerators(matrices,p):
    return [[int(c)%p for c in sp.Matrix(A).charpoly().all_coeffs()] for A in matrices]
