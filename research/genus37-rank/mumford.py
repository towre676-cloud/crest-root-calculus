"""Cantor addition in the Jacobian of a monic odd-degree hyperelliptic curve."""
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import (gf_add, gf_gcdex, gf_monic, gf_mul, gf_neg, gf_exquo,
                                     gf_rem, gf_sub)

IDENTITY=((1,),())


class Jacobian:
    def __init__(self, polynomial, p):
        self.F=list(polynomial)
        self.p=p
        self.g=(len(self.F)-2)//2
        if self.F[0]!=1 or len(self.F)!=2*self.g+2:
            raise ValueError("A monic odd-degree polynomial is required")

    def add(self,left,right):
        if left==IDENTITY:
            return right
        if right==IDENTITY:
            return left
        u,v=map(list,left)
        w,z=map(list,right)
        p=self.p
        mul=lambda a,b:gf_mul(a,b,p,ZZ)
        add=lambda a,b:gf_add(a,b,p,ZZ)
        quo=lambda a,b:gf_exquo(a,b,p,ZZ)
        rem=lambda a,b:gf_rem(a,b,p,ZZ)
        a,b,d1=gf_gcdex(u,w,p,ZZ)
        c,h,d=gf_gcdex(d1,add(v,z),p,ZZ)
        U=quo(mul(u,w),mul(d,d))
        V=rem(quo(add(mul(c,add(mul(mul(a,u),z),mul(mul(b,w),v))),
                      mul(h,add(mul(v,z),self.F))),d),U)
        while len(U)>self.g+1:
            U=gf_monic(quo(gf_sub(self.F,mul(V,V),p,ZZ),U),p,ZZ)[1]
            V=rem(gf_neg(V,p,ZZ),U)
        return tuple(U),tuple(V)

    def mul(self,n,D):
        if n<0:
            D=(D[0],tuple((-z)%self.p for z in D[1]))
            n=-n
        out=IDENTITY
        while n:
            if n&1:
                out=self.add(out,D)
            n//=2
            if n:
                D=self.add(D,D)
        return out

    def check(self,D):
        u,v=map(list,D)
        if not u or u[0]!=1 or len(u)>self.g+1 or len(v)>=len(u):
            raise ValueError("Unreduced Mumford representation")
        p=self.p
        if gf_rem(gf_sub(self.F,gf_mul(v,v,p,ZZ),p,ZZ),u,p,ZZ):
            raise ValueError("Divisor misses the curve")
