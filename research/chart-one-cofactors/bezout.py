"""Exact evaluation of a univariate Bezout row by bounded modular reconstruction."""
from math import isqrt
import sympy as sp


def mod_bezout_at(a,b,point,p):
    """Evaluate the extended Euclidean updates during modular long division."""
    r0,r1 = [v%p for v in a],[v%p for v in b]
    u0,v0,u1,v1 = 1,0,0,1
    while r1:
        d = len(r0)-len(r1)
        inv = pow(r1[0],-1,p)
        value = 0
        for i in range(d+1):
            q = r0[i]*inv%p
            value = (value*point+q)%p
            if q:
                for k in range(1,len(r1)):
                    r0[i+k] = (r0[i+k]-q*r1[k])%p
        rem = r0[max(0,d+1):]
        first = next((i for i,t in enumerate(rem) if t),len(rem))
        r0,r1 = r1,rem[first:]
        u0,u1 = u1,(u0-value*u1)%p
        v0,v1 = v1,(v0-value*v1)%p
    if len(r0)!=1:
        raise ArithmeticError("Noncoprime modular input at an allegedly good prime.")
    inv = pow(r0[0],-1,p)
    return u0*inv%p,v0*inv%p


def exact_at(a,b,x):
    """Return U(x), V(x) for A*U+B*V=1 and the standard strict degree bounds."""
    a_den,a_num = a.clear_denoms(convert=True)
    b_den,b_num = b.clear_denoms(convert=True)
    a_content,A = a_num.primitive()
    b_content,B = b_num.primitive()
    a_scale,b_scale = sp.QQ.convert(a_content)/sp.QQ.convert(a_den),sp.QQ.convert(b_content)/sp.QQ.convert(b_den)
    R = int(A.resultant(B))
    if not R:
        raise ArithmeticError("The primitive norm quotients must be coprime.")
    ai,bi = [list(map(int,f.all_coeffs())) for f in (A,B)]
    m,n = A.degree(),B.degree()
    # Each numerator coefficient is a signed maximal minor of the Sylvester
    # matrix with one column omitted. Hadamard applies to its shifted columns.
    na,nb = [isqrt(sum(t*t for t in v))+1 for v in (ai,bi)]
    xn,xd = int(x.numerator),int(x.denominator)
    h = max(abs(xn),xd)
    bound = max(n*na**(n-1)*nb**m*h**(n-1),m*na**n*nb**(m-1)*h**(m-1))
    modulus,p,U,V = 1,2**31,0,0
    while modulus<=2*bound:
        p = int(sp.nextprime(p))
        if not (R*xd*ai[0]*bi[0])%p:
            continue
        u,v = mod_bezout_at(ai,bi,xn*pow(xd,-1,p)%p,p)
        u = R%p*u*pow(xd,n-1,p)%p
        v = R%p*v*pow(xd,m-1,p)%p
        w = pow(modulus,-1,p)
        U += modulus*((u-U)*w%p)
        V += modulus*((v-V)*w%p)
        modulus *= p
    U = U-modulus if 2*U>modulus else U
    V = V-modulus if 2*V>modulus else V
    return (sp.QQ(U,R*xd**(n-1))/a_scale,sp.QQ(V,R*xd**(m-1))/b_scale)
