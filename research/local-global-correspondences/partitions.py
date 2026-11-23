"""Recover every rational unordered pairing of the six branch points by one split-prime lift."""
from itertools import combinations

import sympy as sp
from sympy.polys.domains import ZZ
from sympy.polys.galoistools import gf_pow_mod

from frobenius import X

T, Z = sp.symbols('theta z')


def sextic_chart(f):
    """Move infinity away from the branch divisor by a rational determinant-one change up to sign."""
    if f.degree()==6:
        return f, (1,0,0,1)
    n=next(i for i in range(6) if f.eval(i))
    g=sum(f.nth(i)*(n*X+1)**i*X**(6-i) for i in range(7))
    return sp.Poly(g,X,domain=sp.QQ), (n,1,1,0)


def split_prime(integral):
    p=13
    while True:
        p=int(sp.nextprime(p))
        h=[int(x)%p for x in integral.all_coeffs()]
        if gf_pow_mod([1,0],p,h,p,ZZ) != [1,0]:
            continue
        roots=[i for i in range(p) if int(integral.eval(i))%p==0]
        if len(roots)!=6:
            raise ArithmeticError('Split-prime test disagrees with the root count')
        for k in range(min(106,p)):
            if len({(a+b+k*a*b)%p for a,b in combinations(roots,2)})==15:
                return p,k,roots


def pairings(left):
    if not left:
        yield ()
        return
    for j in range(1,len(left)):
        for rest in pairings(left[1:j]+left[j+1:]):
            yield ((left[0],left[j]),*rest)


def centered(v,m):
    return (v+m//2)%m-m//2


def mul_list(a,b):
    out=[0]*(len(a)+len(b)-1)
    for i,x in enumerate(a):
        for j,y in enumerate(b):
            out[i+j]+=x*y
    return out


def pair_algebra(f, m, numerator, k, leading):
    """In Q[theta]/m, reconstruct the pair's sum and product and verify the exact norm."""
    d=sp.Poly(m,T,domain=sp.QQ)
    if d.gcd(d.diff()).degree():
        return None
    t=(sp.Poly(numerator,T,domain=sp.QQ)*d.diff().invert(d)).rem(d)
    s=(sp.Poly(T,T,domain=sp.QQ)-k*t).rem(d)
    # The monic input used for lifting has roots leading times the original branch roots.
    s=s.mul_ground(1/leading); t=t.mul_ground(1/leading**2)
    g=X**2-s.as_expr()*X+t.as_expr()
    norm=sp.Poly(sp.resultant(d.as_expr(),g,T),X,domain=sp.QQ)
    if norm!=f.monic():
        return None
    return d, s, t


def stable_pairings(f):
    f, chart=sextic_chart(f)
    leading=f.LC()
    # Leading coefficient a. The polynomial a^5 f(X/a) is integral and monic.
    integral=sp.Poly.from_list([f.nth(i)*leading**(5-i) for i in range(6,-1,-1)],X,domain=sp.ZZ)
    p,k,roots=split_prime(integral)
    R=1+max(abs(int(c)) for c in integral.all_coeffs()[1:])
    H=2*R+k*R**2
    bound=max((1+H)**3,3*R**2*(1+H)**2)
    modulus=p
    # Double p-adic precision. Each starting root is simple.
    deriv=integral.diff()
    while modulus<=2*bound:
        modulus*=modulus
        roots=[(r-int(integral.eval(r))*pow(int(deriv.eval(r)),-1,modulus))%modulus for r in roots]
    if any(int(integral.eval(r))%modulus for r in roots):
        raise ArithmeticError('Hensel root lost precision')
    result=[]
    for partition in pairings(roots):
        sums=[(a+b+k*a*b)%modulus for a,b in partition]
        vals=[(a*b)%modulus for a,b in partition]
        m=[1]
        for v in sums:
            m=mul_list(m,[-v,1])
        interpol=[0,0,0]
        for i,v in enumerate(vals):
            q=mul_list([-sums[(i+1)%3],1],[-sums[(i+2)%3],1])
            for j in range(3):
                interpol[j]+=v*q[j]
        mc=[centered(c,modulus) for c in m]
        nc=[centered(c,modulus) for c in interpol]
        if any(abs(c)>bound for c in mc+nc):
            continue
        sext=pair_algebra(f,sp.Poly.from_list(mc[::-1],T).as_expr(),
                          sp.Poly.from_list(nc[::-1],T).as_expr(),k,leading)
        if sext is not None:
            result.append(sext)
    return f, {'source_chart':list(chart), 'split_prime':p,'pair_parameter':k,
               'coefficient_bound':bound, 'hensel_modulus':modulus,
               'stable_partitions':len(result)}, result


def rem(expr,m):
    return sp.rem(expr,m.as_expr(),T).expand()


def trace_matrix(m):
    sums=[sp.Integer(3)]
    for i in range(1,5):
        a=sum(m.nth(3-j)*sums[i-j] for j in range(1,min(i,4)))
        a+=i*m.nth(3-i) if i<=3 else 0
        sums.append(-a)
    return sp.Matrix(3,3,lambda i,j:sums[i+j])


def coeff_vector(f):
    return sp.Matrix([f.nth(i) for i in range(3)])


def dual_algebra(source, part):
    """Trace-dual construction of the Richelot target and its three conjugate graphs."""
    m,s,t=part
    unit=sp.Matrix([1,0,0])
    basis=sp.Matrix.hstack(unit,-coeff_vector(s),coeff_vector(t))
    determinant=basis.det()
    if not determinant:
        return None
    V=sp.Matrix([1,T,T**2])
    tr=trace_matrix(m)
    if tr.det()!=m.discriminant():
        raise ArithmeticError('Trace form discriminant mismatch')
    dual=(basis.T*tr).inv()*sp.Matrix([1,-2*Z,Z**2])
    ell=(V.T*dual)[0].expand()
    g=X**2-s.as_expr()*X+t.as_expr()
    e=rem(g*ell,m)
    if sp.expand((sp.Matrix([[3,-m.nth(2),tr[1,1]]])*coeff_vector(sp.Poly(e,T)))[0]-(X-Z)**2):
        raise ArithmeticError('Trace-dual equation failed')
    delta2=determinant**2*m.discriminant()
    c=source.LC()
    fd=sp.Poly(-c*delta2*sp.resultant(m.as_expr(),ell,T),Z,domain=sp.QQ)
    proj=((X-Z)**2-e).expand()
    E=sp.Poly(e,T)
    inv=sp.invert(m.diff(),m).as_expr()
    value=-(c/2)*(X-Z)*determinant*m.discriminant()*(E.nth(1)+E.nth(2)*(-m.nth(2)-T))*inv
    return {'m':m,'g':g,'ell':ell,'delta_squared':delta2,'dual':fd,'Q':proj,'B':rem(value,m)}


def check_graph(source, data):
    """Check in the full cubic etale algebra, retaining every factor of its minimal polynomial."""
    m=data['m']
    for factor,e in m.factor_list()[1]:
        for root in [sp.CRootOf(factor.as_expr(),0)]:
            # One embedded root per rational irreducible factor suffices by conjugacy.
            if factor.degree()==1:
                root=-factor.nth(0)/factor.nth(1)
                field=sp.QQ
            else:
                field=sp.QQ.algebraic_field(root)
            q=sp.Poly(data['Q'].subs(T,root), X,Z,domain=field)
            b=sp.Poly(data['B'].subs(T,root), X,Z,domain=field)
            if q.degree(X)!=2 or q.degree(Z)!=2:
                raise ArithmeticError('Generic bidegree differs from (2,2)')
            lhs=sp.Poly(source.as_expr()*data['dual'].as_expr(),X,Z,domain=field)
            if (lhs-b*b).rem(q):
                raise ArithmeticError('The global graph identity failed')
