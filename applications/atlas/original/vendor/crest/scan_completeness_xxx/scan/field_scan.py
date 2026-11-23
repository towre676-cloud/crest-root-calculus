"""Exhaust both pencil types at the rational marked points of a finite field.

python field_scan.py quartic.json 5 2
Input coefficients are integers or exact rational strings, in Part XXVIII format.
Field elements encode polynomial-basis coefficients in base p, constant first.
"""
from collections import Counter
from fractions import Fraction
from itertools import permutations, product
from math import comb
import argparse
import json
from pathlib import Path
import sympy as sp


class Field:
    def __init__(self, p, degree, modulus=None):
        if not sp.isprime(p) or p <= 3 or degree < 1:
            raise ValueError("Require prime characteristic at least 5 and positive extension degree.")
        self.p, self.degree, self.q = p, degree, p**degree
        x = sp.Symbol('x')
        if modulus is not None:
            if len(modulus)!=degree+1 or modulus[-1]!=1:
                raise ValueError("Expected a monic polynomial of the specified degree.")
            test = sp.Poly(sum(v*x**i for i,v in enumerate(modulus)),x,modulus=p)
            if not test.is_irreducible:
                raise ValueError("Reducible extension polynomial.")
            modulus = [c%p for c in modulus]
        elif degree == 1:
            modulus = [0, 1]
        elif degree == 2 and p % 3 == 2:
            modulus = [1, 1, 1]
        else:
            for candidate in product(range(p), repeat=degree):
                if not candidate[0]:
                    continue
                P = sp.Poly(x**degree+sum(c*x**i for i,c in enumerate(candidate)), x, modulus=p)
                if P.is_irreducible:
                    modulus = [*candidate, 1]
                    break
        self.modulus = modulus
        digits = [tuple(n//p**i % p for i in range(degree)) for n in range(self.q)]
        def encode(row):
            return sum((v%p)*p**i for i,v in enumerate(row))
        self.neg = [encode(-v for v in a) for a in digits]
        self.plus = [[encode(u+v for u,v in zip(a,b)) for b in digits] for a in digits]
        def mul(a,b):
            coeffs = [0]*(2*degree-1)
            for i,u in enumerate(a):
                for j,v in enumerate(b):
                    coeffs[i+j] += u*v
            for k in reversed(range(degree,len(coeffs))):
                for j in range(degree):
                    coeffs[k-degree+j] -= coeffs[k]*modulus[j]
            return encode(coeffs[:degree])
        self.times = [[mul(a,b) for b in digits] for a in digits]
        self.inv = [None]+[row.index(1) for row in self.times[1:]]

    def power(self, x, n):
        value = 1
        while n:
            if n%2:
                value = self.times[value][x]
            x = self.times[x][x]
            n //= 2
        return value

    def prod(self,*args):
        value = 1
        for a in args:
            value = self.times[value][a]
        return value

    def add(self,*args):
        value = 0
        for a in args:
            value = self.plus[value][a]
        return value

    def put(self, out, key, val):
        if val:
            v = self.plus[out.get(key,0)][val]
            if v:
                out[key] = v
            else:
                out.pop(key,None)

    def polymul(self, A, B):
        out = {}
        for a,u in A.items():
            for b,v in B.items():
                self.put(out,tuple(i+j for i,j in zip(a,b)),self.times[u][v])
        return out

    def scaled_sum(self, terms):
        result = {}
        for c,poly in terms:
            for key,v in poly.items():
                self.put(result,key,self.times[c][v])
        return result


def read_quartic(data,field):
    result = {}
    for i,j,k,raw in data['terms']:
        if any(type(n)!=int or n<0 for n in (i,j,k)) or i+j+k!=4:
            raise ValueError('Invalid quartic exponents.')
        c = Fraction(raw)
        if c.denominator%field.p==0:
            raise ValueError('Denominator divisible by the chosen prime; supply an integral model.')
        field.put(result,(i,j,k),int(c.numerator*pow(c.denominator,-1,field.p) % field.p))
    if not result:
        raise ValueError('Zero reduction.')
    return result


def smooth(F,field):
    X = sp.symbols('X Y Z')
    eq = sum(v*sp.prod(x**i for x,i in zip(X,e)) for e,v in F.items())
    for k,z in enumerate(X):
        affine = [eq.diff(x).subs(z,1) for x in X]
        if list(sp.groebner(affine,*(x for x in X if x!=z),modulus=field.p)) != [1]:
            raise ValueError('The projective reduction must be smooth over the algebraic closure.')


def evaluate(F,vec,field):
    return field.add(*(field.prod(c,*(field.power(v,i) for v,i in zip(vec,m))) for m,c in F.items()))


def points(F,field):
    q = field.q
    for z,yvals,xvals in [(1,range(q),range(q)),(0,[1],range(q)),(0,[0],[1])]:
        for y in yvals:
            for x in xvals:
                if not evaluate(F,(x,y,z),field):
                    yield (x,y,z)


def frame_for(F,P,K):
    p = K.p
    jac = []
    for k in range(3):
        derivative = {tuple(i-(j==k) for j,i in enumerate(m)):K.times[c][m[k]%p]
                      for m,c in F.items() if m[k]}
        jac.append(evaluate(derivative,P,K))
    i,j,n = next(frame for frame in permutations(range(3)) if P[frame[2]] and jac[frame[1]])
    u,v = [K.times[P[k]][K.inv[P[n]]] for k in (i,j)]
    d = K.neg[K.times[jac[i]][K.inv[jac[j]]]]
    # X_i=L+uN, X_j=H+dL+vN, X_n=N. Expand with Hasse coefficients.
    T = {}
    for exps,c in F.items():
        I,J = exps[i],exps[j]
        for a in range(I+1):
            for b in range(J+1):
                for h in range(J-b+1):
                    coeff = comb(I,a)*comb(J,b)*comb(J-b,h) % p
                    K.put(T,(a+b,h),K.prod(c,coeff,K.power(u,I-a),K.power(d,b),K.power(v,J-b-h)))
    assert T.get((0,0),0)==T.get((1,0),0)==0 and T[0,1]
    one = [{tuple(int(k==j) for k in range(3)):1} for j in range(3)]
    N = one[n]
    L = K.scaled_sum([(1,one[i]),(K.neg[u],N)])
    H = K.scaled_sum([(1,one[j]),(K.neg[v],N),(K.neg[d],L)])
    return T,(i,j,n),(L,H,N)


def osculation_data(T,s,K):
    # f(u,v)=F(s-u,-1,v) in the adapted L,H,N frame.
    f = {}
    for (i,j),co in T.items():
        d = 4-i-j
        for r in range(i+1):
            K.put(f,(r,d),K.prod(co,((-1)**(r+j)*comb(i,r))%K.p,K.power(s,i-r)))
    return f


def hasse(f, row, c, n, K):
    return K.add(*(K.prod(v,comb(d,n)%K.p,K.power(c,d-n))
                   for (i,d),v in f.items() if i==row and d>=n))


def transversal_parameters(f,K):
    # h_j and j_j are Hasse coefficients in v at c.
    d0,d1,d2,d3 = (f.get((0,j),0) for j in range(4))
    if not d3:raise ValueError("The transverse fibre polynomial must be cubic.")
    # Coefficients of 2hh''-h'^2, evaluated with Horner's rule.
    division_polynomial=[K.prod(3,d3,d3),K.prod(4%K.p,d3,d2),K.prod(6%K.p,d3,d1),
                         K.prod(12%K.p,d3,d0),K.add(K.prod(4%K.p,d2,d0),K.neg[K.prod(d1,d1)])]
    for c in range(K.q):
        value=division_polynomial[0]
        for coeff in division_polynomial[1:]:
            value=K.plus[K.times[c][value]][coeff]
        if value:continue
        h,h1,h2,h3 = (hasse(f,0,c,n,K) for n in range(4))
        if not h:continue
        assert K.prod(4%K.p,h,h2)==K.times[h1][h1]
        j,j1,j2 = (hasse(f,1,c,n,K) for n in range(3))
        # Taylor coefficient of y^2/r in the square residual,
        # after clearing 4h². The slope of b is 12 h² h3.
        numer = K.add(K.prod(4%K.p,h,h,j2),K.neg[K.prod(2,h,h1,j1)],K.prod(h1,h1,j))
        denom = K.prod(12%K.p,h,h,h3)
        assert denom
        yield K.times[K.neg[numer]][K.inv[denom]],c


def transverse_G(f,b,c,K):
    # G=t^6 f(1/t, c+b/t+z/t²)
    G = {}
    for (i,d),co in f.items():
        for j in range(d+1):
            for r in range(d-j+1):
                key=(6-i-r-2*j,j)
                K.put(G,key,K.prod(co,(comb(d,j)*comb(d-j,r))%K.p,K.power(c,d-j-r),K.power(b,r)))
    assert all(0<=i<=6 and 0<=j<=3 for i,j in G)
    return G


def tangent_parameters(T,K):
    A,B,C,D = [T.get(m,0) for m in [(0,1),(2,0),(1,1),(3,0)]]
    if B:
        a = K.prod(4%K.p,B,K.inv[K.prod(3,A)])
        b = K.add(K.prod(4%K.p,C,K.inv[K.prod(3,A)]),K.neg[K.prod(3,D,K.inv[K.prod(2,B)])])
        return a,b


def tangent_G(T,a,b,K):
    G = {}
    for (i,j),co in T.items():
        k=4-i-j
        for v in range(k+1):
            for r in range(k-v+1):
                h=k-v-r
                K.put(G,(i+2*h+r,v),K.prod(co,((-1)**k*comb(k,v)*comb(k-v,r))%K.p,
                                           K.power(a,h),K.power(b,r)))
    return G


def square_shape(G,K):
    g=G.get((6,0),0)
    if not g:return
    div=lambda i,j: K.times[G.get((i,j),0)][K.inv[g]]
    half=K.inv[2]
    e2=K.times[div(5,0)][half]
    e1=K.prod(half,K.add(div(4,0),K.neg[K.prod(e2,e2)]))
    f1=K.prod(half,div(4,1))
    e0=K.add(K.prod(half,div(3,0)),K.neg[K.prod(e1,e2)])
    f0=K.add(K.prod(half,div(3,1)),K.neg[K.prod(f1,e2)])
    S={(3,0):1,(2,0):e2,(1,0):e1,(1,1):f1,(0,0):e0,(0,1):f0}
    residual=K.scaled_sum([(g,K.polymul(S,S)),(K.neg[1],G)])
    if any(i for i,j in residual):
        return
    k0,k1,k2,k3=[K.times[residual.get((0,j),0)][K.inv[g]] for j in range(4)]
    disc=K.add(K.prod(k1,k1,k2,k2), K.neg[K.prod(4%K.p,k3,k1,k1,k1)],
               K.neg[K.prod(4%K.p,k2,k2,k2,k0)],K.neg[K.prod(27%K.p,k3,k3,k0,k0)],
               K.prod(18%K.p,k3,k2,k1,k0))
    if not k3 or not disc:
        raise ArithmeticError('Degenerate target survived exact square test at a smooth source.')
    return {'g':g,'shape':[e2,e1,f1,e0,f0],'target_constant_first':[k0,k1,k2,k3]}


def scan_point(F,P,K,stats):
    T,frame,(L,H,N)=frame_for(F,P,K)
    result=[]
    def accept(kind,params,U,D,answer):
        result.append((U,D,dict(point=P,frame=frame,kind=kind,parameters=params,**answer)))
    for s in range(K.q):
        f=osculation_data(T,s,K)
        stats['transverse_lines']+=1
        for b,c in transversal_parameters(f,K):
            stats['transverse_conics']+=1
            ans=square_shape(transverse_G(f,b,c,K),K)
            if ans:
                line=K.scaled_sum([(1,L),(s,H)])
                aux=K.scaled_sum([(1,N),(b,line),(c,H)])
                num=K.scaled_sum([(K.neg[1],K.polymul(aux,H))])
                accept('I',[s,b,c],num,K.polymul(line,line),ans)
    params=tangent_parameters(T,K)
    if params:
        stats['tangent_conics']+=1
        a,b=params
        ans=square_shape(tangent_G(T,a,b,K),K)
        if ans:
            num=K.scaled_sum([(K.neg[1],K.polymul(H,N)),(K.neg[a],K.polymul(L,L)),(K.neg[b],K.polymul(L,H))])
            accept('II',[a,b],num,K.polymul(H,H),ans)
    return result


def normal_form(A,F,K):
    A=dict(A)
    n=max(F)
    leadinv=K.inv[F[n]]
    rem={}
    while A:
        m=max(A)
        if any(a<b for a,b in zip(m,n)):
            rem[m]=A.pop(m)
            continue
        coeff=K.neg[K.times[A[m]][leadinv]]
        for s,c in F.items():
            K.put(A,tuple(a-b+k for a,b,k in zip(m,n,s)),K.times[coeff][c])
    return rem


def pencil_relation(pair1,pair2,F,K):
    def span(pair):
        a,b=pair
        return (K.polymul(b,b),K.polymul(a,b),K.polymul(a,a))
    basis={}
    for s in span(pair1):
        for t in span(pair2):
            col=normal_form(K.polymul(s,t),F,K)
            while col and max(col) in basis:
                m=max(col)
                col=K.scaled_sum([(1,col),(K.neg[col[m]],basis[m])])
            if not col:return True
            m=max(col)
            basis[m]=K.scaled_sum([(K.inv[col[m]],col)])
    return False


def classify_marks(triples,F,K):
    classes=[]
    for U,D,hit in triples:
        pair=(U,D)
        for base,marks in classes:
            if pencil_relation(base,pair,F,K):
                marks.append(hit)
                break
        else:
            classes.append((pair,[hit]))
    def frob(pair):
        return tuple({m:K.power(c,K.p) for m,c in h.items()} for h in pair)
    perm=[]
    for base,marks in classes:
        image=frob(base)
        matches=[i for i,(other,_) in enumerate(classes) if pencil_relation(image,other,F,K)]
        if len(matches)!=1:
            raise ArithmeticError('Frobenius image missing or ambiguous.')
        perm.append(matches[0])
    assert sorted(perm)==list(range(len(perm)))
    assert all(len(marks)<=2 for _,marks in classes)
    return {'marked_hits':sum(len(m) for _,m in classes),
            'classes_with_a_rational_mark':len(classes),'p_frobenius_permutation':perm,
            'p_frobenius_fixed_in_scan':sum(i==j for i,j in enumerate(perm)),
            'exhausts_Fp_defined_classes':K.degree%2==0,
            'marks_grouped_by_quotient':[marks for _,marks in classes]}


def scan(data,p,degree):
    K=Field(p,degree)
    F=read_quartic(data,K)
    smooth(F,K)
    stats=Counter({'points':0,'transverse_lines':0,'transverse_conics':0,'tangent_conics':0})
    triples=[]
    for P in points(F,K):
        stats['points']+=1
        triples.extend(scan_point(F,P,K,stats))
    return {'prime':p,'degree':degree,'q':K.q,'modulus_constant_first':K.modulus,
            'operations':dict(stats),**classify_marks(triples,F,K)}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('quartic',type=Path)
    parser.add_argument('p',type=int)
    parser.add_argument('degree',type=int)
    args=parser.parse_args()
    print(json.dumps(scan(json.loads(args.quartic.read_text()),args.p,args.degree),indent=2))
