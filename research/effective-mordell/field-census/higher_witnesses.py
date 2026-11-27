"""Produce exact point witnesses with quartic-base abscissa for the rank-four fields 7 and 11."""
from fractions import Fraction
from itertools import product
from math import gcd
from octic_fields import FIELD_DATA, Z, ROOT, p
def witnesses(d):
    seeds=[Z**0*0]
    for den in range(1,7):
        for a,b,c,e in product(range(-2,3),repeat=4):
            if gcd(gcd(gcd(abs(a),abs(b)),gcd(abs(c),abs(e))),den)!=1:continue
            u=(a+b*Z+c*Z**2+e*Z**3)/den
            v=ROOT(FIELD_DATA,d*(u**6-1))
            if v and u!=0:
                seeds+= [u,1/u]
    found=set()
    for u in seeds:
        for j in range(6):
            x=u*Z**(2*j)
            for rs in ROOT(FIELD_DATA,d*(x**6-1)):
                for y in (rs,-rs):
                    if y*y!=d*(x**6-1):raise ArithmeticError('Invalid point')
                    if y:found.add((str(p.lift(x)),str(p.lift(y))))
    return sorted(found)
if __name__=='__main__':
    import json
    from pathlib import Path
    pts={str(d):witnesses(d) for d in (7,11)}
    print([(d,len(x)) for d,x in pts.items()])
    Path('higher-octic-witnesses.json').write_text(json.dumps(pts,indent=2)+'\n')
