"""Evaluate all six original-generator cofactors of a chart-support identity.

The circuit has characteristic-zero polynomial coefficients. This evaluator
also supports modular evaluation of those coefficients. Its divisions are
exact polynomial quotients, evaluated at a point where their denominators are
nonzero. The checker in proof.py verifies the underlying identities globally.
"""
import argparse
import json
from decimal import Decimal
from pathlib import Path
import sympy as sp

from chart import R,S,c,POSITIONS,build
from solve import ResidueField,lifted,specialize
from proof import crt_targets,decode,original_rows,read_u,verify


ORDER = ((2,0),(2,1),(2,2),(1,0),(1,1),(1,2))


def evaluate(data,point,proof,at,prime=None):
    chart = build(data,point)
    summary = verify(data,point,proof)
    if prime is not None and (not sp.isprime(prime) or prime<=3):
        raise ValueError("Choose a prime greater than three.")
    K = sp.QQ if prime is None else sp.GF(prime)
    def val(x):
        rat = sp.Rational(x)
        return K.convert(int(rat.p))/K.convert(int(rat.q))
    sv,bv,cv = map(val,at)
    def ev(poly,coords=(cv,sv)):
        return sum((val(a)*product(x**e for x,e in zip(coords,k)) for k,a in poly.items()),K.zero)
    def product(factors):
        result = K.one
        for a in factors:
            result *= a
        return result
    def evs(poly):
        return sum((val(a)*sv**j for (j,),a in poly.rep.to_dict().items()),K.zero)
    point3 = (bv,cv,sv)
    _,_,_,original = original_rows(chart)
    observations = [ev(original[idx],point3) for idx in ORDER]
    q,h = ev(chart.q),ev(chart.h)
    D,B = ev(chart.denominator),ev(chart.numerator)
    L = D*bv-B
    if not all((h,q,L,D)):
        raise ValueError("Choose an evaluation point away from h*q*(D*b-B)=0.")
    w,u = ev(chart.table[0,3]),ev(chart.table[0,1])
    H = h*h
    unit = lambda k:[K.one if i==k else K.zero for i in range(6)]
    def combine(pairs):
        out = [K.zero]*6
        for a,row in pairs:
            for i,v in enumerate(row): out[i] += a*v
        return out
    # Each vector below represents H times the indicated polynomial.
    Vq = combine([(1/(48*w*w),unit(2))])
    VL = combine([(1/val(16),unit(5)),(bv*u/16,unit(2))])
    Vr = []
    for (i,l),row in zip(POSITIONS,chart.rows):
        n = original[i,l]
        k = 6-i-2*l
        replaced = D**k*ev(n,(B/D,cv,sv))
        u0 = (replaced-ev(row))/q
        v0 = (D**k*ev(n,point3)-replaced)/L
        Vr.append(combine([(H*D**k,unit(ORDER.index((i,l)))),(-u0,Vq),(-v0,VL)]))
    selected,norm_vectors = [],[]
    for recipe,co in zip(proof["norm_rows"],proof["norm_cofactors"]):
        row = sum((i*r for i,r in zip(recipe,chart.rows)),R.zero)
        ad = decode(co)
        n = (row*ad).rem([chart.q])
        selected.append(sp.Poly.from_dict({(j,):a for (i,j),a in n.items()},(S,),domain=sp.QQ))
        quotient = (ev(row)*ev(ad)-ev(n))/q
        norm_vectors.append(combine([(ev(ad)*a,v) for a,v in zip(recipe,Vr)]+[(-quotient,Vq)]))
    factors = [read_u(r["s_factor"]) for r in proof["fibres"]]
    gcd = sp.Poly(1,S,domain=sp.QQ)
    for f,rec in zip(factors,proof["fibres"]):
        gcd *= f**rec["s_multiplicity_in_eliminant"]
    def change(poly):
        return sp.Poly.from_dict({k:K.to_sympy(val(a)) for k,a in poly.rep.to_dict().items()},(S,),domain=K)
    quotients = [n.exquo(gcd) for n in selected]
    quotient_norms = [change(n) for n in quotients]
    if any(a.degree()!=b.degree() for a,b in zip(quotients,quotient_norms)):
        raise ValueError("The coefficient Bezout row lost leading degree at this prime.")
    if len(selected)==1:
        if quotient_norms[0].degree()!=0:
            raise ArithmeticError("Single eliminant cofactor is nonconstant.")
        weights = [1/val(quotient_norms[0].LC())]
    elif any(n.degree()==0 for n in quotient_norms):
        i = next(i for i,n in enumerate(quotient_norms) if n.degree()==0)
        weights = [K.zero,K.zero]
        weights[i] = 1/val(quotient_norms[i].LC())
    elif prime is None:
        from bezout import exact_at
        weights = exact_at(quotients[0],quotients[1],sv)
    else:
        # Equivalently the first column of the inverse Sylvester matrix.
        a,b,d = quotient_norms[0].gcdex(quotient_norms[1])
        if d.degree()!=0:
            raise ValueError("The coefficient Bezout row lost rank at this prime.")
        def evaluate_k(poly):
            return sum((v*sv**i for (i,),v in poly.rep.to_dict().items()),K.zero)
        weights = [evaluate_k(a)/evaluate_k(d),evaluate_k(b)/evaluate_k(d)]
    Ve = combine(zip(weights,norm_vectors))
    targets = crt_targets(proof["fibres"])
    if proof["geometric_marks"]==1:
        positive = next(r for r in proof["fibres"] if r["geometric_marks"])
        b0 = decode(positive["parameters"]["b"])
        if not b0.is_ground:
            raise ArithmeticError("A unique geometric triple should be rational.")
        targets.append(("b_numerator",chart.numerator-chart.denominator*b0))
    outputs = {}
    for name,target in targets:
        delta,tail,wproduct = [K.zero]*6,K.one,K.one
        for f,rec,local in zip(factors,proof["fibres"],proof["fibre_cofactors"]):
            fv = evs(f)
            if not fv:
                raise ValueError("Choose an s outside the projected eliminant.")
            field = ResidueField(f)
            spec = lambda a:specialize(a,f,field)
            gd = spec(decode(rec["fibre_gcd"]))
            tv = target if rec["geometric_marks"] else R.one
            local_multiplier = (spec(chart.h)**4*spec(tv)).exquo(gd)
            coeff = [ev(lifted(local_multiplier*spec(decode(v)))) for v in local]
            X = h**4*ev(tv)
            generators = [q,*[ev(v) for v in chart.rows]]
            W = (X-sum((a*g for a,g in zip(coeff,generators)),K.zero))/fv
            VW = combine(zip(coeff,[Vq,*Vr]))
            for _ in range(rec["s_multiplicity_in_eliminant"]):
                delta = combine([(X,delta),(tail,VW)])
                tail *= fv*W
                wproduct *= W
        v = combine([(K.one,delta),(wproduct,Ve)])
        left = h**summary["h_exponent"]*ev(target)**summary["target_exponent"]
        if name=="b_numerator" and summary["target_exponent"]==1:
            # n-D*b0 + L = D(b-b0), and D=12 w h².
            v = combine([(1/(12*w),v),(h**(summary["h_exponent"]-2)/(12*w),VL)])
            left = h**(summary["h_exponent"]+2)*(bv-ev(b0))
            name = "b_coordinate"
        if sum((x*y for x,y in zip(v,observations)),K.zero)!=left:
            raise ArithmeticError("Original six-generator cofactor identity failed.")
        def output(a):
            if prime is not None:
                return int(a)
            n,den = int(a.numerator),int(a.denominator)
            return str(Decimal(n)) if den==1 else str(Decimal(n))+"/"+str(Decimal(den))
        outputs[name]={"target_value":output(left),"cofactor_values":[output(a) for a in v]}
    return {"original_generator_order":[list(k) for k in ORDER],"evaluation_s_b_c":[str(t) for t in at],
            "prime":prime,"identities":outputs}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("quartic",type=Path)
    parser.add_argument("certificate",type=Path)
    parser.add_argument("point",nargs=3)
    parser.add_argument("--at",nargs=3,default=["11","5","7"],metavar=("s","b","c"))
    parser.add_argument("--prime",type=int)
    args = parser.parse_args()
    print(json.dumps(evaluate(json.loads(args.quartic.read_text()),args.point,
                              json.loads(args.certificate.read_text()),args.at,args.prime),indent=2))
