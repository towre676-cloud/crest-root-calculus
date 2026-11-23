"""Evaluate the signed fibre coefficient matrices in two groups of tensor indices."""
from functools import lru_cache
from itertools import product
from pathlib import Path
import sys

import numpy as np
import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/"cubic-tower-geometric-separation"))
from affine_frobenius import characteristic, f_power, half_column, multiply


def half_rows(p,indices):
    words=list(product((1,2),repeat=len(indices)))
    by_weight=[[] for _ in range(len(indices)+1)]
    for j,src in enumerate(words):
        for i,dst in enumerate(words):
            by_weight[(i^j).bit_count()].append((j,i,half_column(p,indices,src,dst)))
    return by_weight


def packed(rows, p, extra=0):
    polynomials=[multiply(v,f_power(p,extra),p) if extra else v for j,i,v in rows]
    data=np.zeros((len(rows),max(map(len,polynomials))),dtype=np.float64)
    for i,poly in enumerate(polynomials):
        data[i,:len(poly)]=poly
    return data


def assemble(p,S):
    k=len(S)
    split=k//2
    left,right=half_rows(p,S[:split]),half_rows(p,S[split:])
    right_size=1 << (k-split)
    word_weights=np.array([k+i.bit_count() for i in range(1<<k)],dtype=np.int64)
    weights=np.array([(word_weights+3*((word_weights+sign)%2))//2 for sign in (0,1)])
    offsets=np.column_stack((np.zeros(2,dtype=np.int64),np.cumsum(weights,axis=1)))
    matrices=[np.zeros((int(s[-1]),int(s[-1])),dtype=np.int64) for s in offsets]

    @lru_cache(None)
    def right_coefficients(odd,extra):
        # Reverse once so each convolution uses a forward matrix slice.
        return np.ascontiguousarray(packed(right[odd],p,extra)[:,::-1])

    for odd_left,L in enumerate(left):
        lc=packed(L,p)
        jL=np.array([j for j,i,poly in L])[:,None]
        iL=np.array([i for j,i,poly in L])[:,None]
        for odd_right,R in enumerate(right):
            jR=np.array([j for j,i,poly in R])
            iR=np.array([i for j,i,poly in R])
            src = jL*right_size+jR
            dst = iL*right_size+iR
            for d in (0,1):
                sign=(word_weights[src]+d)%2
                rc=right_coefficients(odd_right,d*((p-1)//2)+(d+odd_left+odd_right)//2)
                if min(lc.shape[1],rc.shape[1])*(p-1)**2 >= 2**53:
                    raise OverflowError("The convolution dot product exceeds exact binary64 range")
                nr=weights[sign,dst]
                nc=weights[sign,src]
                row=offsets[sign,dst]; col=offsets[sign,src]
                for b in range(1,int(nc.max())+1):
                    for a in range(1,int(nr.max())+1):
                        q=b*p-a
                        start=max(0,q-rc.shape[1]+1)
                        stop=min(lc.shape[1],q+1)
                        if start>=stop:
                            continue
                        r0=rc.shape[1]-1-q+start
                        dots=lc[:,start:stop]@rc[:,r0:r0+stop-start].T
                        values=dots.astype(np.int64)%p
                        valid=(a<=nr)&(b<=nc)
                        for e in (0,1):
                            use=valid&(sign==e)
                            matrices[e][row[use]+a-1,col[use]+b-1]=values[use]
    right_coefficients.cache_clear()
    return matrices


def signed_blocks(p,indices):
    S=tuple(indices)
    k=len(S); dim=k*(1<<(k-1)) if k else 0
    if len(set(S))!=k or any(type(i) is not int or i<0 for i in S):
        raise ValueError("Support indices must be distinct nonnegative integers")
    if not 1<=k<=7 or not sp.isprime(p) or not max(3,dim)<p<=10000:
        raise ValueError("Use 1 through 7 distinct good indices and a prime exceeding their signed dimension")
    matrices=assemble(p,S)
    chi=[int(sp.legendre_symbol(-(12*i+1)*pow(12,-1,p)%p,p)) for i in S]
    if 0 in chi: raise ValueError("Nonregular fibre at infinity")
    result=[]
    for e,mat in enumerate(matrices):
        coeff=characteristic(mat,p)[:dim+1]
        rank=1<<(k-1)
        trace=(sp.prod(1+c for c in chi)+(-1)**e*sp.prod(1-c for c in chi))//2
        for eigen,mult in ((1,(rank+trace)//2),(-1,(rank-trace)//2)):
            for _ in range(mult):
                for j in range(1,len(coeff)):
                    coeff[j]=(coeff[j]+eigen*coeff[j-1])%p
        result.append(coeff)
    return result
