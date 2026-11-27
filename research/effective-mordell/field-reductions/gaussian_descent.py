"""Local and rank-zero tests on v^2=7(x^6-1) over Q(i)."""
from itertools import product
from functools import lru_cache
from cypari2 import Pari

p=Pari()
reps=(1,2,3,6,7,14,21,42)
def gaussian_rep(n):
    # 2 and i agree modulo squares, since 2i=(1+i)^2.
    return (0,n//2) if n%2==0 else (n,0)

def local_table(q,m):
    elements=list(product(range(m),repeat=2))
    def mul(u,v):return ((u[0]*v[0]-u[1]*v[1])%m,(u[0]*v[1]+u[1]*v[0])%m)
    def add(u,v):return ((u[0]+v[0])%m,(u[1]+v[1])%m)
    def sub(u,v):return ((u[0]-v[0])%m,(u[1]-v[1])%m)
    def nonunit(u):
        if q==2:return (u[0]+u[1])%2==0
        return u[0]%q==0 and u[1]%q==0
    one=(1,0)
    pairs=[(a,one) for a in elements]+[(one,b) for b in elements if nonunit(b)]
    squares={mul(a,a) for a in elements}
    good=set()
    for u,v in pairs:
        uu,vv,uv=mul(u,u),mul(v,v),mul(u,v)
        A=add(add(uu,vv),uv);B=sub(add(uu,vv),uv);D=mul((7,0),sub(uu,vv))
        aa=[a for a in reps if mul(gaussian_rep(a),A) in squares]
        bb=[b for b in reps if mul(gaussian_rep(b),B) in squares]
        for a,b in product(aa,bb):
            if mul(mul(gaussian_rep(a),gaussian_rep(b)),D) in squares:good.add((a,b))
    return good

def squarefree(n):return int(p.core(n))

@lru_cache(None)
def rankbounds(kind,k):
    if kind=='even':coef=[0,k,0,k*k,0]
    else:coef=[0,9*k,0,30*k*k,36*k**3]
    E=p.ellinit(coef)
    bounds=[]
    for F in (E,E.elltwist(-4)):
        r=F.ellrank()
        bounds.append((int(r[0]),int(r[1])))
    return tuple(map(sum,zip(*bounds)))

@lru_cache(None)
def finite_abscissas(kind,k):
    nf=p.nfinit('w^2+1')
    v=p('x')
    co=[0,k,0,k*k,0] if kind=='even' else [0,9*k,0,30*k*k,36*k**3]
    E=p.ellinit(co,nf)
    order,orders,generators=E.elltors()
    points=[p([0])]
    for n,T in zip(orders,generators):
        points=[E.elladd(P,E.ellmul(T,j)) for P in points for j in range(int(n))]
    if len(set(map(str,points)))!=int(order):raise ArithmeticError('incomplete torsion')
    def roots(c):return list(p.nfroots(nf,v**2-c))
    out=set()
    for P in points:
        if len(P)==1:continue
        if kind=='even':candidates=roots(P[0]/k)
        elif P[0]!=0:candidates=[1+6*k/P[0]]
        else:candidates=[]
        for X in candidates:
            if roots(7*(X**6-1)):out.add(str(p.lift(X)))
    return out
def verify():
    pairs=set(product(reps,repeat=2))
    for q,m in ((2,16),(3,9),(7,49),(11,11),(19,19)):
        pairs &= local_table(q,m)
        print('After Gaussian modulus',m,sorted(pairs))
    survived=[]
    for a,b in sorted(pairs):
        cases=[('even',squarefree(a*b)),('mixed',squarefree(7*b)),('mixed',squarefree(7*a))]
        ranks=[rankbounds(kind,k) for kind,k in cases]
        if all(hi>0 for lo,hi in ranks):survived.append((a,b))
        else:
            j=next(j for j,(_,hi) in enumerate(ranks) if hi==0)
            xs=finite_abscissas(*cases[j])
            if not xs <= {'0','1','-1','2','-2','1/2','-1/2'}:
                raise ArithmeticError(f'unexpected Gaussian abscissa {xs}')
        print((a,b),cases,ranks)
    if set(survived)!={(3,7),(7,3)}:raise ArithmeticError('residual squareclass pair')
    print('Only the pair (7,3) and its reflection (3,7) remain for nonreal abscissas.')
    return pairs,survived
if __name__=='__main__':verify()
