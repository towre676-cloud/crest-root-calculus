"""Classify C6 over L_d=Q(zeta12,sqrt(d)) when the CM rank is at most two.

The main result records coordinates of every ordinate in k=Q(zeta12) with
x in k on v^2=d(x^6-1). Its corresponding original ordinate is v/sqrt(d).
Exceptional quadratic abscissas are recovered by the companion divisor script.
"""
import argparse
from fractions import Fraction
from itertools import product
from math import ceil, isqrt
from cypari2 import Pari

p = Pari()
X, W = p("x"), p("w")
FIELD = W**4-W**2+1
FIELD_DATA = p.nfinit(FIELD)
Z = p.Mod(W, FIELD)
I = Z**3
SQRT3 = Z+1/Z
O = p([0])
DIVISION = p("(E,P,n)->{my(Q);if(ellisdivisible(E,P,n,&Q),Q,[])}")
ROOT = p("(N,a)->{my(b);if(nfeltissquare(N,a,&b),[nfbasistoalg(N,b)],[])}")

def rat_height(value):
    q=Fraction(str(value))
    return max(abs(q.numerator),q.denominator)

def torsion_points(E):
    order,orders,generators=E.elltors()
    out=[O]
    for n,P in zip(orders,generators):
        out=[E.elladd(Q,E.ellmul(P,j)) for Q in out for j in range(int(n))]
    if len({str(Q) for Q in out})!=int(order):
        raise ArithmeticError("Incomplete torsion enumeration")
    return out

def certify_minimum_height(E,t):
    """Every nontorsion point of E_t(Q) has doubled canonical height >1/4."""
    T={str(P) for P in torsion_points(E)}
    candidates={Fraction(a,b) for a in range(-7,8) for b in range(1,8)}
    for q in candidates:
        x=p(str(q))*t
        for y in E.ellordinate(x):
            P=p([x,y])
            if str(P) not in T:
                H=rat_height(E.ellmul(P,4)[0]/t)
                if H<=330:raise ArithmeticError("Uniform rational height lower bound needs refinement")

def excluded_locally(E,t,P,prime):
    # If P=prime*Q, multiplying by #E(F_r)/prime kills its image at a good prime r.
    for r in map(int,p.primes(500)):
        if (6*t)%r==0:
            continue
        if any(int(z.denominator())%r==0 for z in P):
            continue
        F=p.ellinit([0,-t**3],r)
        N=int(F.ellcard())
        if N%prime==0 and len(F.ellmul(P,N//prime))>1:return True
    return False

def primitive(E,t,G):
    if int(E.elltors()[0])!=2 or not E.ellisoncurve(G):
        raise ArithmeticError("Unexpected rational character model")
    certify_minimum_height(E,t)
    T=p([t,0])
    while True:
        bound=isqrt(4*(rat_height(G[0]/t).bit_length()+1))
        divided=False
        for l in map(int,p.primes(max(bound,2))):
            if l>bound:continue
            for P in ([G,E.elladd(G,T)] if l==2 else [G]):
                if excluded_locally(E,t,P,l):
                    continue
                Q=DIVISION(E,P,l)
                if len(Q)==0:
                    continue
                if E.ellmul(Q,l)!=P:raise ArithmeticError("Rational division failed")
                G=Q;divided=True;break
            if divided:break
        if not divided:return G

def companion(G):
    """Degree-three isogeny from E_-d to E_3d."""
    u,v=G
    b=v**2-u**3
    return p([u+4*b/u**2,v*(1-8*b/u**3)])

def lift_coordinates(t,d,P):
    # Q(zeta12)=Q(i,sqrt(3)).  Ratio sqrt(t/d) belongs to this fixed field.
    c={-1:I,3:SQRT3,1:1,-3:I*SQRT3}[t//d]
    return p([P[0]*d/t,P[1]*d/(t*c)])

def halve(E,P,b):
    f=X**4-8*b*X-4*P[0]*(X**3+b)
    for u in p.nfroots(FIELD_DATA,f):
        roots=ROOT(FIELD_DATA,u**3+b)
        for r in roots:
            for v in (r,-r):
                Q=p([u,v])
                if E.ellmul(Q,2)==P:return Q

def saturate(E,d,basis,torsion):
    rows=[(Fraction(1),Fraction(0)),(Fraction(0),Fraction(1))]
    index=1
    while True:
        found=False
        for a,b in ((1,0),(0,1),(1,1)):
            Q=E.elladd(E.ellmul(basis[0],a),E.ellmul(basis[1],b))
            for T in torsion:
                P=halve(E,E.elladd(Q,T),-d**3)
                if P is None:continue
                place=0 if a else 1
                basis[place]=P
                rows[place]=tuple((a*c+b*e)/2 for c,e in zip(*rows))
                index*=2;found=True;break
            if found:break
        if index>16:raise ArithmeticError("Trace index bound exceeded")
        if not found:break
    return rows,index

def classify(d):
    if d<=1 or d%3==0 or int(p.core(d))!=d:
        raise ValueError("Use a squarefree positive radicand >1 prime to 3.")
    p.setrand(d)  # Each receipt is reproducible independently of the worker schedule.
    rational={}
    for t in (d,-d,3*d,-3*d):
        E=p.ellinit([0,-t**3]);r=E.ellrank()
        rational[t]=(E,tuple(map(int,r[:2])),list(r[3]))
    bounds=tuple(sum(r[1][j] for r in rational.values()) for j in (0,1))
    ans=dict(d=d,rank=list(bounds))
    if bounds not in ((0,0),(2,2)):return dict(ans,status='partial')
    E=p.ellinit([0,-d**3],FIELD_DATA)
    torsion=torsion_points(E)
    dirs=[(t,*data) for t,data in rational.items() if data[1][1]]
    def write_twist(P,pts):
        if len(P)==1:return
        roots=ROOT(FIELD_DATA,P[0]/d)
        for r in roots:
            for x in (r,-r):
                y=P[1]/d
                if d*(x**6-1)!=y*y:raise ArithmeticError("Twist fibre equation")
                pts.add((str(p.lift(x)),str(p.lift(y))))
    points=set()
    for T in torsion:write_twist(T,points)
    if not dirs:return dict(ans,status='complete',twist_points=sorted(points),torsion=len(torsion))
    if len(dirs)!=2 or any(r!=(1,1) for t,F,r,G in dirs):
        return dict(ans,status='partial',reason='non-paired character ranks')
    if dirs[0][0]!=-d or dirs[1][0]!=3*d:
        return dict(ans,status='partial',reason='unexpected character signs')
    raw = dirs[0][3]
    if not raw: raise ArithmeticError('No generator of the first rational character')
    first=primitive(dirs[0][1],-d,raw[0])
    # The companion is always present, including when the rational point search returns no generator.
    second=primitive(dirs[1][1],3*d,companion(first))
    lowers=[]
    for P,(t,F,_,_) in zip((first,second),dirs):
        H=rat_height(F.ellmul(P,4)[0]/t)
        lower=max(Fraction(1,4),(Fraction(2,3)*(H.bit_length()-1)-2)/16)
        lowers.append(lower)
    basis=[lift_coordinates(-d,d,first),lift_coordinates(3*d,d,second)]
    if any(not E.ellisoncurve(P) for P in basis):raise ArithmeticError("Character lift")
    rows,index=saturate(E,d,basis,torsion)
    (a,b),(c,e)=rows
    det=abs(a*e-b*c)
    if det!=Fraction(1,index):raise ArithmeticError("Lattice determinant")
    # If H(P)<7, orthogonality implies |coefficient_i| < sqrt(7/lower_i).
    caps=[isqrt((7//L))+1 for L in lowers]
    lims=(ceil((caps[0]*abs(e)+caps[1]*abs(c))/det),
          ceil((caps[0]*abs(b)+caps[1]*abs(a))/det))
    checked=0
    for j,k in product(*(range(-s,s+1) for s in lims)):
        u,v=j*a+k*c,j*b+k*e
        if u*u*lowers[0]+v*v*lowers[1]>=7:continue
        P=E.elladd(E.ellmul(basis[0],j),E.ellmul(basis[1],k))
        for T in torsion:
            write_twist(E.elladd(P,T),points)
            checked+=1
    return dict(ans,status='complete',twist_points=sorted(points),torsion=len(torsion),
                rational_primitive=[[str(c) for c in P] for P in (first,second)],
                rational_height_lower=[str(l) for l in lowers],basis_index=index,
                basis_change=[[str(c) for c in row] for row in rows],lattice_checks=checked)

if __name__=='__main__':
    import json
    cmd=argparse.ArgumentParser(description=__doc__)
    cmd.add_argument('d',type=int)
    print(json.dumps(classify(cmd.parse_args().d),indent=2))
