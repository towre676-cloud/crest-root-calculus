"""Finite search on the unresolved descending equation; a negative search is not a proof."""
import argparse
from math import gcd,isqrt

def search(H):
    tests=0;solutions=[]
    for c in range(1,H+1):
        for a in range(H+1):
            for b in range(1,H+1):
                if gcd(gcd(a,b),c)!=1:continue
                A=a**3-3*a*b*b;B=3*a*a*b-b**3
                r=7*(A*A-B*B-c**6);j=14*A*B
                norm=r*r+j*j
                z=isqrt(norm); tests+=1
                if z*z!=norm or (r+z)%2:continue
                U=isqrt((r+z)//2);V=isqrt((z-r)//2)
                if 2*U*U==r+z and 2*V*V==z-r and 2*U*V==abs(j):
                    solutions.append((a,b,c))
    return tests,solutions
if __name__=='__main__':
    app=argparse.ArgumentParser(description=__doc__);app.add_argument('bound',nargs='?',default=65,type=int)
    H=app.parse_args().bound
    if H<=0:app.error('supply a positive coefficient bound')
    checks,solutions=search(H)
    print('Primitive nonreal abscissas tested up to sign and conjugation:',checks)
    print('Solutions in this coefficient box:',solutions)
