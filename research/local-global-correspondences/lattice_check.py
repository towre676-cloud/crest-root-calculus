"""Direct index-m² Hermite-lattice count of all and primitive genus-two isogeny kernels."""
from itertools import product
from math import gcd
import sympy as sp


def count_lattices(m):
    index=m*m; total=primitive=0
    for a in sp.divisors(index):
        for b in sp.divisors(index//a):
            for c in sp.divisors(index//(a*b)):
                d=index//(a*b*c)
                if a*c%m: continue
                for u,h in product(range(a),range(c)):
                    if a*h%m or u*c%m or (u*h+b*d)%m: continue
                    for v,w,r,s in product(range(a),range(b),range(a),range(b)):
                        if (v*h+w*d-c*r)%m: continue
                        total+=1
                        primitive+=int(gcd(a,b,c,d,u,h,v,w,r,s,m)==1)
    return primitive,total
