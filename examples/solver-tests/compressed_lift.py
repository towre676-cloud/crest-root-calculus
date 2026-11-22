"""Lift both kinds of rational marked pencils one p-adic digit at a time.

At a transverse flag the free digits are the residual-point parameter u and
the fibre slope s. At a tangent flag only u varies. The polar conic is forced.
Works over Z/p^n for integral quartics smooth at p>3 and simple local frames.
"""
from itertools import product
from math import comb


def value(F,x,y,derivative):
    e,f=derivative
    return sum(c*comb(i,e)*comb(j,f)*x**(i-e)*y**(j-f)
               for (i,j),c in F.items() if i>=e and j>=f)


def taylor(F,u,v,d,M):
    return {(i,j):sum(comb(a+j,j)*value(F,u,v,(i-a,a+j))*pow(d,a,M)
                      for a in range(i+1))%M for i in range(5) for j in range(5-i)}


def fibre(T,s,M):
    f={}
    for (i,j),t in T.items():
        k=4-i-j
        for n in range(i+1):
            m=n,k
            f[m]=(f.get(m,0)+t*comb(i,n)*(-1)**(n+j)*pow(s,i-n,M))%M
    return f


def deriv(f,row,d,x,M):
    return sum(c*comb(i,d)*pow(x,i-d,M) for (n,i),c in f.items() if n==row and i>=d)%M


def hensel_implicit(x0,level,p,f,df):
    M=level*p
    y=f(x0)%M
    assert y%level==0
    return (x0-(y//level)*pow(df(x0)%p,-1,p)*level)%M


def residual(F,kind,state,digits,level,p):
    """Compute all residual digits and uniquely restored dependent parameters."""
    u0,v0,s0,c0=state
    M=level*p
    u=u0+level*digits[0]
    v=hensel_implicit(v0,level,p,lambda y:value(F,u,y,(0,0)),lambda y:value(F,u,y,(0,1)))
    d=(-value(F,u,v,(1,0))*pow(value(F,u,v,(0,1)),-1,M))%M
    T=taylor(F,u,v,d,M)
    assert T[0,0]==T[1,0]==0
    G={}
    def put(key,v):
        G[key]=(G.get(key,0)+v)%M
    if kind=='I':
        s=s0+level*digits[1]
        f=fibre(T,s,M)
        def div3(c):
            h,h1,h2=(deriv(f,0,m,c,M) for m in range(3))
            return 4*h*h2-h1*h1
        def slope(c):
            return 12*deriv(f,0,0,c,M)*deriv(f,0,3,c,M)
        c=hensel_implicit(c0,level,p,div3,slope)
        h,h1,h2,h3=(deriv(f,0,m,c,M) for m in range(4))
        j,j1,j2=(deriv(f,1,m,c,M) for m in range(3))
        b=-(4*h*h*j2-2*h*h1*j1+h1*h1*j)*pow(12*h*h*h3,-1,M)%M
        for (i,dv),co in f.items():
            for j in range(dv+1):
                for r in range(dv-j+1):
                    put((6-i-r-2*j,j),co*comb(dv,j)*comb(dv-j,r)*pow(c,dv-j-r,M)*pow(b,r,M))
        lifted=(u,v,s,c)
    else:
        A,B,C,D=(T[x] for x in [(0,1),(2,0),(1,1),(3,0)])
        a=4*B*pow(3*A,-1,M)%M
        b=(4*C*pow(3*A,-1,M)-3*D*pow(2*B,-1,M))%M
        for (i,j),co in T.items():
            k=4-i-j
            for vz in range(k+1):
                for r in range(k-vz+1):
                    h=k-vz-r
                    put((i+2*h+r,vz),co*(-1)**k*comb(k,vz)*comb(k-vz,r)*pow(a,h,M)*pow(b,r,M))
        lifted=(u,v,None,None)
    g=G[6,0]
    gm=pow(g,-1,M)
    half=pow(2,-1,M)
    t=lambda i,j:G.get((i,j),0)*gm%M
    e2=t(5,0)*half%M
    e1=(t(4,0)-e2**2)*half%M
    f1=t(4,1)*half%M
    e0=(t(3,0)*half-e2*e1)%M
    f0=(t(3,1)*half-e2*f1)%M
    S={(3,0):1,(2,0):e2,(1,0):e1,(1,1):f1,(0,0):e0,(0,1):f0}
    for (i,j),a in S.items():
        for (k,l),b in S.items():put((i+k,j+l),-g*a*b)
    rem=[]
    for k in range(1,7):
        for j in range(4):
            c=G.get((k,j),0)
            if c%level:
                raise ValueError('Input state fails at the starting precision.')
            rem.append(c//level)
    return rem,lifted


def reduce_rows(A,b,p):
    augmented=[[*row,y] for row,y in zip(A,b)]
    n=len(A[0])
    pivot=[]
    for j in range(n):
        pos=next((i for i in range(len(pivot),len(A)) if augmented[i][j]%p),None)
        if pos is None:continue
        i=len(pivot);pivot.append(j)
        augmented[i],augmented[pos]=augmented[pos],augmented[i]
        scalar=pow(augmented[i][j]%p,-1,p)
        augmented[i]=[v*scalar%p for v in augmented[i]]
        for k,row in enumerate(augmented):
            if k!=i:
                scalar=row[j]%p
                if scalar:augmented[k]=[(v-scalar*w)%p for v,w in zip(row,augmented[i])]
    inconsistent=any(row[-1]%p and not any(row[:-1]) for row in augmented)
    return len(pivot),inconsistent,augmented,pivot


def lift_step(F,kind,state,level,p):
    n=2 if kind=='I' else 1
    z=[0]*n
    e0,_=residual(F,kind,state,z,level,p)
    col=[]
    for j in range(n):
        unit=z.copy();unit[j]=1
        ej,_=residual(F,kind,state,unit,level,p)
        col.append([(a-b)%p for a,b in zip(ej,e0)])
    A=list(zip(*col))
    rank,bad,rref,piv=reduce_rows(A,[(-v)%p for v in e0],p)
    if bad:
        return rank,[]
    free=[i for i in range(n) if i not in piv]
    out=[]
    for vals in product(range(p),repeat=len(free)):
        digits=z.copy()
        for i,v in zip(free,vals):digits[i]=v
        for i,j in enumerate(piv):
            row=rref[i]
            digits[j]=(row[-1]-sum(row[k]*digits[k] for k in free))%p
        check,new=residual(F,kind,state,digits,level,p)
        assert not any(check)
        out.append(new)
    return rank,out


def setup(terms,hit,p):
    i,j,n=hit['frame']
    F={}
    for *m,co in terms:
        key=m[i],m[j]
        F[key]=F.get(key,0)+co
    P=hit['point']
    unit=pow(P[n],-1,p)
    u,v=P[i]*unit%p,P[j]*unit%p
    if hit['kind']=='I':
        s,b,c=hit['parameters']
        return F,(u,v,s,c)
    return F,(u,v,None,None)
