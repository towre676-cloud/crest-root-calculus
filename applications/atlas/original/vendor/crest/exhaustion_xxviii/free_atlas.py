"""Exact point-free incidence atlas for cubic elliptic quotients of a smooth quartic.

Input JSON follows Part XXVII's format {"terms":[[i,j,k,coefficient],...]}, or
a single dictionary containing that list. 'emit' constructs one or all charts
as textual Singular programs (std over Q). 'solve' uses SymPy's exact basis
algorithm and prints a geometric rational univariate representation for each
chart. A full free-point solve may need substantial time and memory.
"""
import argparse
import itertools
import json
from functools import lru_cache
from pathlib import Path
import sympy as S

X,Y,Z = COORDS = S.symbols('X Y Z')
u,v,d,s,b,c,a = S.symbols('u v d s b c a')
e2,e1,f1,e0,f0,g,h = S.symbols('e2 e1 f1 e0 f0 g hinv')
t,z = S.symbols('t z')
FRAMES = tuple(itertools.permutations(range(3)))
KINDS = ('I','II')


def input_quartic(data):
    value = 0
    for i,j,k,co in data['terms']:
        if any(type(n)!=int or n<0 for n in (i,j,k)) or i+j+k!=4:
            raise ValueError('Expected nonnegative exponents of total degree four.')
        value += S.Rational(co)*X**i*Y**j*Z**k
    poly = S.Poly(value,*COORDS,domain=S.QQ)
    if poly.is_zero: raise ValueError('Zero quartic.')
    return poly


def smooth_check(F):
    for i,V in enumerate(COORDS):
        others = tuple(W for W in COORDS if W!=V)
        jac = [F.diff(W).as_expr().subs(V,1) for W in COORDS]
        if list(S.groebner(jac,*others,domain=S.QQ)) != [1]:
            raise ValueError('The projective quartic is singular over the algebraic closure.')


def translated_coefficients(F,frame):
    """F(L+uZ,H+dL+vZ,Z) in a permuted frame, via the Taylor formula."""
    x,y,n = [COORDS[i] for i in frame]
    f = F.as_expr().subs(n,1)
    result = {}
    for i in range(5):
        for j in range(5-i):
            terms = []
            for k in range(i+1):
                term=S.diff(f,x,i-k,y,k+j).subs({x:u,y:v},simultaneous=True)
                if term:
                    terms.append(term*d**k*S.binomial(i,k))
            result[i,j] = S.Add(*terms)/(S.factorial(i)*S.factorial(j))
    return result


@lru_cache(maxsize=2)
def residual_table(kind):
    """Universal binary parameter substitutions. Keys are (t exponent,z exponent)."""
    if kind=='I':
        lin,tan,norm = -t+s*t**2,-t**2,z+b*t+c*t**2
        shift = -2
    elif kind=='II':
        lin,tan,norm = t,S.Integer(1),-a*t**2-b*t-z
        shift = 0
    else:
        raise ValueError('Chart must be I or II.')
    result={}
    for i in range(5):
        for j in range(5-i):
            if (i,j) in ((0,0),(1,0)):
                continue
            entry=S.Poly(lin**i*tan**j*norm**(4-i-j),t,z)
            for (k,l),coef in entry.terms():
                result.setdefault((k+shift,l),[]).append((i,j,coef))
    assert all(0<=k<=6 and 0<=l<=3 for k,l in result)
    square = t**3+e2*t**2+(e1+f1*z)*t+e0+f0*z
    return result,dict(S.Poly(square**2,t,z).terms())


def chart(F,frame,kind,fix=None):
    """The normalized finite chart, with its target equation g*w² = K(z)."""
    if frame not in FRAMES or kind not in KINDS: raise ValueError('Unknown chart.')
    fix = fix or {}
    T = {key:value.subs(fix,simultaneous=True)
         for key,value in translated_coefficients(F,frame).items()}
    table,square = residual_table(kind)
    G = {(k,l):S.Add(*(T[i,j]*coef for i,j,coef in entries)).subs(fix,simultaneous=True)
         for (k,l),entries in table.items()}
    sq = {key:value.subs(fix,simultaneous=True) for key,value in square.items()}
    get = lambda data,key:data.get(key,S.Integer(0))
    high = sorted({key for key in (G.keys()|sq.keys()) if key[0]>0},reverse=True)
    eq = [T[0,0],T[1,0],(h*g*T[0,1]-1).subs(fix,simultaneous=True)]
    eq += [get(G,key)-g.subs(fix)*get(sq,key) for key in high]
    K = sum((g.subs(fix)*get(sq,(0,j))-get(G,(0,j)))*z**j for j in range(4))
    variables=tuple(x for x in (h,g,e2,e1,f1,e0,f0,*((s,b,c) if kind=='I' else (a,b)),d,u,v) if x not in fix)
    return variables,[entry for entry in eq if entry!=0],K


def singular_text(variables,ideal):
    def expression(q):
        return S.sstr(q).replace('**','^')
    return '\n'.join([
        '// QQ, geometric support of the indicated smooth-source chart.',
        'ring R=0,('+','.join(map(str,variables))+'),dp;',
        'ideal I='+',\n'.join(expression(p) for p in ideal)+';',
        'option(redSB);', 'ideal B=std(I);','print(B);',
        'print("dimension="+string(dim(B)));','quit;',''])


def geometric_rur(polynomials,variables):
    """Trace-pairing RUR of the support, valid also for nonreduced ideals."""
    G = S.groebner(polynomials,*variables,order='grevlex',domain=S.QQ)
    if list(G)==[1]:
        return {'dimension':-1,'length':0,'points':0}
    if not G.is_zero_dimensional:
        raise ValueError('Positive-dimensional ideal; cannot construct a finite RUR.')
    counts = len(variables)
    lms = [p.LM(order=G.order).exponents for p in G.polys]
    def standard(exps): return not any(all(a<=b for a,b in zip(lm,exps)) for lm in lms)
    todo = {(0,)*counts}
    mons = set()
    while todo:
        exps=todo.pop()
        if exps in mons or not standard(exps): continue
        mons.add(exps)
        for i in range(counts):
            vec=list(exps);vec[i]+=1
            todo.add(tuple(vec))
    mons=sorted(mons,key=lambda m:(sum(m),m))
    terms=[S.prod(x**n for x,n in zip(variables,m)) for m in mons]
    length=len(terms)
    def matrix(q):
        cols=[]
        for term in terms:
            reduced = G.reduce(S.expand(term*q))[1]
            coeff=S.Poly(reduced,*variables,domain=S.QQ)
            cols.append([coeff.coeff_monomial(m) for m in mons])
        return S.Matrix(cols).T
    vectors=[matrix(x) for x in variables]
    bases=[S.eye(length)]
    for mon in mons[1:]:
        i=next(i for i in range(counts) if mon[i]>0)
        div=list(mon);div[i]-=1
        bases.append(vectors[i]*bases[mons.index(tuple(div))])
    tr=lambda A,B:sum(A[i,j]*B[j,i] for i in range(A.rows) for j in range(A.cols))
    r=S.Matrix(length,length,lambda i,j: tr(bases[i],bases[j])).rank()
    T=S.Symbol('T')
    for n in itertools.count():
        M=sum((n**i*A for i,A in enumerate(vectors)),S.zeros(length))
        poly=M.charpoly(T).as_poly().sqf_part().monic()
        if poly.degree()==r: break
    powers=[S.eye(length)]
    for _ in range(2*r-2):
        powers.append(powers[-1]*M)
    pairing=S.Matrix(r,r,lambda i,j:powers[i+j].trace())
    rhs=S.Matrix(r,counts,lambda i,j:tr(powers[i],vectors[j]))
    weights=pairing.inv()*rhs
    coords=[sum(weights[k,i]*T**k for k in range(r)) for i in range(counts)]
    # Verify all original equations symbolically at every returned point.
    for p in polynomials:
        remainder=S.rem(p.subs(dict(zip(variables,coords)),simultaneous=True),poly.as_expr(),T)
        if remainder!=0: raise ArithmeticError('RUR residual is nonzero.')
    return {'dimension':0,'length':length,'points':r,'linear_form':str(sum(n**i*x for i,x in enumerate(variables))),
            'minimal_polynomial':str(poly.as_expr()),'coordinates':dict(zip(map(str,variables),map(str,coords)))}


def main():
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('action',choices=('emit','solve','certify','verify'))
    ap.add_argument('input',type=Path)
    ap.add_argument('--frame',choices=tuple(''.join(map(str,f)) for f in FRAMES))
    ap.add_argument('--kind',choices=KINDS)
    ap.add_argument('--output',type=Path,required=True)
    args=ap.parse_args()
    F=input_quartic(json.loads(args.input.read_text()));smooth_check(F)
    args.output.mkdir(exist_ok=True,parents=True)
    for frame in FRAMES:
        code=''.join(map(str,frame))
        for kind in KINDS:
            if args.frame and code!=args.frame or args.kind and kind!=args.kind: continue
            variables,ideal,_=chart(F,frame,kind)
            if args.action=='emit':
                location=args.output/f'chart_{code}_{kind}.sing'
                location.write_text(singular_text(variables,ideal))
            elif args.action in ('certify','verify'):
                from ideal_certificate import generate, verify
                location=args.output/f'chart_{code}_{kind}_certificate.json'
                if args.action=='certify':
                    certificate=generate(ideal,variables)
                    status=verify(ideal,variables,certificate)
                    location.write_text(json.dumps(certificate,separators=(',',':'))+'\n')
                else:
                    status=verify(ideal,variables,json.loads(location.read_text()))
                print(json.dumps(status),flush=True)
            else:
                result=geometric_rur(ideal,variables)
                location=args.output/f'chart_{code}_{kind}.json'
                location.write_text(json.dumps(result,indent=2)+'\n')
            print(location,flush=True)


if __name__=='__main__':
    main()
