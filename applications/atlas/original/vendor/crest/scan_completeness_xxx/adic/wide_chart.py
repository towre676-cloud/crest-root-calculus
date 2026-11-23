"""Integral lift charts whose two line patches cover neighborhoods of all modular marks."""
from functools import lru_cache
from pathlib import Path
import sys
import sympy as sp

for base in Path(__file__).resolve().parents:
    if (base/'exhaustion_xxviii').is_dir():
        sys.path.insert(0,str(base/'exhaustion_xxviii'))
        break
import free_atlas as fa


@lru_cache(maxsize=1)
def wide_table():
    t,z,s,a,b=fa.t,fa.z,fa.s,fa.a,fa.b
    table={}
    for i in range(5):
        for j in range(5-i):
            if (i,j) in ((0,0),(1,0)):
                continue
            power=i+2*j-2
            assert power>=0
            term=t**i*(1-s*t)**power*(-z-a*t*t-b*t)**(4-i-j)
            for key,coeff in sp.Poly(term,t,z).terms():
                table.setdefault(key,[]).append((i,j,coeff))
    assert all(k<=6 and n<=3 for k,n in table)
    return table


def chart(F,frame,kind):
    if kind=="I":
        variables,eq,target=fa.chart(F,frame,"I")
        return variables,eq
    if kind not in ("II","wide"):
        raise ValueError("Expected transverse I or tangent II.")
    T=fa.translated_coefficients(F,frame)
    _,square=fa.residual_table("II")
    table=wide_table()
    eq=[T[0,0],T[1,0],fa.h*fa.g*T[0,1]-1]
    for key in sorted(set(table)|set(square),reverse=True):
        if key[0]:
            coeff=sum(T[i,j]*value for i,j,value in table.get(key,[]))
            eq.append(coeff-fa.g*square.get(key,0))
    variables=(fa.h,fa.g,fa.e2,fa.e1,fa.f1,fa.e0,fa.f0,fa.s,fa.a,fa.b,fa.d,fa.u,fa.v)
    return variables,[f for f in eq if f!=0]

