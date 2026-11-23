"""Exact recognition of genus-two Richelot correspondences from two binary sextics."""
from itertools import permutations, product
from math import gcd, isqrt, lcm

import sympy as sp

from frobenius import X, polynomial, squareclass

W = sp.Symbol('w')


def vector(q):
    return [q.nth(i) for i in (2, 1, 0)]


def quadratics(f):
    """Unordered rational quadratic-factor triples of the binary sextic, including infinity."""
    unit, entries = f.factor_list()
    if any(e != 1 or q.degree()>2 for q,e in entries):
        return []
    whole = [q for q,e in entries if q.degree()==2]
    linear = [q for q,e in entries if q.degree()==1]
    if f.degree()==5:
        linear.append(sp.Poly(1, X, domain=sp.QQ))  # the binary factor Z
    def match(left):
        if not left:
            yield []
            return
        for j in range(1, len(left)):
            for rest in match(left[1:j]+left[j+1:]):
                yield [left[0]*left[j],*rest]
    triples = []
    for pairs in match(linear):
        g = whole+pairs
        g = [unit*g[0],*g[1:]]
        if sp.prod(g) != f:
            raise ArithmeticError('Factorization lost its scalar')
        triples.append(g)
    return triples


def bracket(f, g):
    return f.diff()*g-f*g.diff()


def dual(quads):
    delta = sp.Matrix([vector(q) for q in quads]).det()
    if not delta:
        return None   # the polarized quotient is a product
    h = [bracket(quads[(i+1)%3],quads[(i+2)%3]) for i in range(3)]
    if any(bracket(h[(i+1)%3],h[(i+2)%3]) != 2*delta*quads[i] for i in range(3)):
        raise ArithmeticError('The double-bracket identity failed')
    # Q(x,x) = -G2(x)H2(x) and sum G_i(x)H_i(z) = delta (x-z)^2.
    pairing = sum(g.as_expr()*r.as_expr().subs(X,W) for g,r in zip(quads,h))
    if sp.expand(pairing-delta*(X-W)**2):
        raise ArithmeticError('The polarized quadratic identity failed')
    return delta, h, sp.prod(h).mul_ground(-1/delta)


def gram(q):
    A = sp.Matrix([vector(f) for f in q])
    return A*sp.Matrix([[0,0,-2],[0,1,0],[-2,0,0]])*A.T


def rational_sqrt(x):
    x=sp.Rational(x)
    if x<0:
        return None
    n, d = isqrt(int(x.p)), isqrt(int(x.q))
    return sp.Rational(n,d) if n*n==x.p and d*d==x.q else None


def primitive(values):
    den = lcm(*(int(sp.denom(x)) for x in values))
    ints = [int(x*den) for x in values]
    g = gcd(*ints)
    if next(x for x in ints if x)<0:
        g=-g
    return tuple(x//g for x in ints)


def transform(q, matrix, degree, variable=X):
    a,b,c,d=matrix
    return sum(q.nth(i)*(a*variable+b)**i*(c*variable+d)**(degree-i)
               for i in range(degree+1)).expand()


def match_frames(h, t, sieve_primes):
    """Recover every rational PGL2 matrix taking this ordered target pairing to h."""
    H, K = gram(h), gram(t)
    ratios = [sp.Rational(H[i,i],K[i,i]) for i in range(3)]
    # Reduction of squared determinant-normalized scales rejects impossible frames.
    for p in sieve_primes:
        for r in ratios:
            n,d=int(r.p)%p,int(r.q)%p
            if n and d and pow(n*d,(p-1)//2,p)==p-1:
                return []
    roots=[rational_sqrt(r) for r in ratios]
    if any(r is None for r in roots):
        return []
    frames=set()
    inv = sp.Matrix([vector(q) for q in h]).inv()
    target = sp.Matrix([vector(q) for q in t])
    for signs in product((1,-1),repeat=2):
        scales = [roots[0], signs[0]*roots[1], signs[1]*roots[2]]
        if any(H[i,j] != scales[i]*scales[j]*K[i,j] for i,j in ((0,1),(0,2),(1,2))):
            continue
        L = inv*sp.diag(*scales)*target
        vv = L*sp.Matrix([X**2,X,1])
        if sp.expand(vv[0]*vv[2]-vv[1]**2):
            continue
        n,d=sp.fraction(sp.cancel(vv[1]/vv[2]))
        n,d=sp.Poly(n,X),sp.Poly(d,X)
        if max(n.degree(),d.degree())>1:
            continue
        matrix = primitive([n.nth(1),n.nth(0),d.nth(1),d.nth(0)])
        a,b,c,e=matrix
        if a*e==b*c:
            continue
        if all(sp.Poly(transform(q,matrix,2),X).monic() == r.monic() for q,r in zip(h,t)):
            frames.add(matrix)
    return sorted(frames)


def terms(expr):
    return [[i,j,int(c.p),int(c.q)] for (i,j),c in sp.Poly(expr,X,W,domain=sp.QQ).terms()]


def build_graph(f, target, quads, delta, h, matrix):
    transformed=[transform(r,matrix,2,W) for r in h]
    target_w=target.as_expr().subs(X,W)
    s=sp.cancel(-sp.prod(transformed)/(delta*target_w))
    if not s or not s.is_Rational:
        raise ArithmeticError('The transported binary sextic is not proportional to the target')
    g=[q.as_expr() for q in quads]
    a,b,c,d=matrix
    graph=(g[0]*transformed[0]+g[1]*transformed[1]).expand()
    product_y=((X*(c*W+d)-(a*W+b))*g[0]*transformed[0]).expand()
    residual=s*f.as_expr()*target_w-product_y**2
    if sp.div(residual,graph,X,W)[1]!=0:
        raise ArithmeticError('The correspondence fails the two curve equations')
    det=a*d-b*c
    return {'source_quadratics_descending':[list(map(str,vector(q))) for q in quads],
            'delta':str(delta), 'target_mobius_matrix':list(matrix),
            'y_scale_squared':str(s), 'field_squareclass':squareclass(int(s.p*s.q)),
            'graph_xw':terms(graph), 'graph_yv_numerator':terms(product_y),
            'pullback_matrix_divided_by_sqrt_scale':
                [[str(sp.Rational(v,det**2)) for v in row] for row in ((a,-c),(-b,d))]}


def correspondences(source, target, primes):
    """All pointwise-rational Richelot kernels with rational projective target identification."""
    sources, targets=quadratics(source),quadratics(target)
    graphs=[]
    products=0
    frames_checked=0
    for quads in sources:
        data=dual(quads)
        if data is None:
            products+=1
            continue
        delta,h,_=data
        matrices=set()
        for triple in targets:
            if not sp.Matrix([vector(q) for q in triple]).det():
                continue
            for ordered in permutations(triple):
                frames_checked+=1
                matrices.update(match_frames(h,ordered,primes))
        for matrix in sorted(matrices):
            graphs.append(build_graph(source,target,quads,delta,h,matrix))
    return {'source_pointwise_rational_kernels':len(sources),
            'product_quotients':products, 'ordered_frames_examined':frames_checked,
            'graphs_up_to_sign':graphs}
