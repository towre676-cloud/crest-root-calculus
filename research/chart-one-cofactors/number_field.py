"""Solve chart I over a prescribed algebraic mark field, with rational descent."""
import argparse
from itertools import count,permutations
import json
from math import comb
from pathlib import Path

import sympy as sp
from sympy.polys.rings import ring

from chart import R,C,S,X,c,s,build_from_frame
from proof import decode,read_u
from solve import (select_norms,coefficients,lifted,profiled_gcd,sparse,
                   ResidueField,quotient_targets,good_prime_for_gcd)


def expansion(value):
    """Increasing coefficients of the power basis of the mark field."""
    return list(reversed(value.to_list()))


def coordinates(values,field):
    v = field.zero
    for coeff in reversed(values):
        v = v*field.unit+field.convert(sp.Rational(coeff))
    return v


def frame(data):
    f = read_u(data["field"])
    if f.degree()<2 or f.LC()!=1 or not f.is_irreducible:
        raise ValueError("The specified field requires a monic irreducible polynomial of degree at least two.")
    K = sp.QQ.alg_field_from_poly(f,alias="theta")
    F = sp.Poly(0,*X,domain=sp.QQ)
    for i,j,k,v in data["terms"]:
        if any(type(n)!=int or n<0 for n in (i,j,k)) or i+j+k!=4:
            raise ValueError("Expected rational quartic terms.")
        F += sp.Poly(sp.Rational(v)*X[0]**i*X[1]**j*X[2]**k,*X,domain=sp.QQ)
    for n in range(3):
        aff = [x for i,x in enumerate(X) if i!=n]
        if list(sp.groebner([F.diff(x).as_expr().subs(X[n],1) for x in X],*aff,domain=sp.QQ))!=[1]:
            raise ValueError("The projective curve must be smooth.")
    point = [coordinates(row,K) for row in data["point"]]
    if len(point)!=3 or not any(point):
        raise ValueError("Expected a nonzero projective triple.")
    def ev(poly):
        return sum((K.convert(a)*point[0]**i*point[1]**j*point[2]**k
                    for (i,j,k),a in poly.rep.to_dict().items()),K.zero)
    jac = [ev(F.diff(x)) for x in X]
    if ev(F) or not any(jac):
        raise ValueError("The mark must be a smooth point of the quartic in the given field.")
    i,j,n = next(row for row in permutations(range(3)) if point[row[2]] and jac[row[1]])
    u,v,d = point[i]/point[n],point[j]/point[n],-jac[i]/jac[j]
    adapted = {}
    for exp,co in F.rep.to_dict().items():
        I,J = exp[i],exp[j]
        for a in range(I+1):
            for b in range(J+1):
                for h in range(J-b+1):
                    key = a+b,h,4-a-b-h
                    term = K.convert(co*comb(I,a)*comb(J,b)*comb(J-b,h))
                    term *= u**(I-a)*d**b*v**(J-b-h)
                    adapted[key] = adapted.get(key,K.zero)+term
    L,H,N = sp.symbols("L H N")
    T = sp.Poly.from_dict(adapted,(L,H,N),domain=K)
    out = {"axes":[i,j,n],"u":list(map(str,expansion(u))),"v":list(map(str,expansion(v))),
           "d":list(map(str,expansion(d)))}
    A,*_ = ring("c,s",K)
    return build_from_frame(T,out,A),f


def shift_norm(e,j,m):
    # Here the first variable of R holds theta and the second holds the new slope parameter U.
    mod = R.from_dict({(i,0):a for (i,),a in m.rep.to_dict().items()})
    shifted = R.zero
    for a in e.rep.to_list():
        coeff = R.from_dict({(i,0):t for i,t in enumerate(expansion(a)) if t})
        shifted = (shifted*(s-j*c)+coeff).rem([mod])
    n = mod.resultant(shifted)
    rational = sp.Poly.from_dict(n.to_dict(),(S,),domain=sp.QQ)
    return shifted,rational


def parts(p):
    out = {}
    for (i,j),v in p.items():
        out.setdefault(i,{})[j,] = v
    return out


def conjugate_fibres(chart,m,selected):
    for j in count():
        norm_rows = [shift_norm(n,j,m) for _,_,n,_ in selected]
        norms = [n for t,n in norm_rows]
        E = norms[0].monic()
        for n in norms[1:]:
            E = E.gcd(n).monic()
        witnesses = []
        for f,mult in E.factor_list()[1]:
            f = f.monic()
            field = ResidueField(f)
            # Determine whether the projection separates the conjugate fields on this fibre.
            original = parts(norm_rows[0][0])
            ep = sp.Poly.from_dict({(i,):sp.Poly.from_dict(v,(S,),domain=sp.QQ).rem(f).as_expr()
                                    for i,v in original.items()},(C,),domain=field)
            field_eq = sp.Poly.from_dict(m.as_dict(),(C,),domain=field)
            root = field_eq.gcd(ep)
            if root.degree()>1:
                break
            if root.degree()!=1:
                raise ArithmeticError("An eliminant root lost its prescribed embedding.")
            theta = sp.Poly(-root.nth(0),S,domain=sp.QQ).rem(f)
            if f.degree()%m.degree():
                raise ArithmeticError("The slope residue field cannot contain the mark field.")
            witnesses.append((f,mult,theta))
        else:
            return j,norms,E,witnesses


def fiber_specialization(f,theta,slope):
    field = ResidueField(f)
    one = sp.Poly(1,S,domain=sp.QQ)
    theta_powers,slope_powers = [one],[one]
    def spec(p):
        if not p:return sp.Poly(0,C,domain=field)
        while len(slope_powers)<=p.degree(1):
            slope_powers.append((slope_powers[-1]*slope).rem(f))
        out = {}
        for (i,j),co in p.items():
            a = expansion(co)
            while len(theta_powers)<len(a):
                theta_powers.append((theta_powers[-1]*theta).rem(f))
            v = sum((q*b for q,b in zip(a,theta_powers)),one.zero)
            out[i,] = out.get((i,),one.zero)+v*slope_powers[j]
        return sp.Poly.from_dict({i:a.rem(f).as_expr() for i,a in out.items()},(C,),domain=field)
    return field,spec


def run(data,details=False):
    chart,m = frame(data)
    selected = select_norms(chart)
    shift,norms,E,witnesses = conjugate_fibres(chart,m,selected)
    result = []
    for f,mult,theta in witnesses:
        slope = (sp.Poly(S,S,domain=sp.QQ)-shift*theta).rem(f)
        field,spec = fiber_specialization(f,theta,slope)
        inputs = [spec(x) for x in (chart.q,*chart.rows)]
        h = spec(chart.h)
        d0,bezout = profiled_gcd(inputs)
        good = d0.exquo(d0.gcd(h**4)).monic()
        if good.gcd(h).degree()!=0 or good.gcd(good.diff()).degree()!=0:
            raise ArithmeticError("Boundary removal or separability failed.")
        record = {"factor":coefficients(f),"multiplicity":mult,"slope":coefficients(slope),
                  "field_generator":coefficients(theta),"c":sparse(lifted(good)),
                  "geometric_marks_per_embedding":f.degree()*good.degree()//m.degree(),
                  "parameters":quotient_targets(chart,f,field,good,spec)}
        if details:
            record["fibre_gcd"] = sparse(lifted(d0))
            record["fibre_cofactors"] = [sparse(lifted(b)) for b in bezout]
        result.append(record)
    out = {"frame":chart.frame,"slope_shift":shift,"norm_rows":[r for _,r,_,_ in selected],
           "rational_norm_degrees":[p.degree() for p in norms],"eliminant_degree":E.degree(),
           "coprimality_prime":good_prime_for_gcd(norms,E),"fibres":result,
           "geometric_marks_per_embedding":sum(x["geometric_marks_per_embedding"] for x in result)}
    if details:
        out["mark_field_norm_cofactors"] = [
            [[i,j,list(map(str,expansion(v)))] for (i,j),v in sorted(u.items())]
            for _,_,_,u in selected]
    return out


def verify(data,proof):
    """Replay supplied cofactor rows without finding the chart's solutions anew."""
    chart,m = frame(data)
    if proof["frame"]!=chart.frame:
        raise ValueError("Unexpected adapted frame.")
    j = proof["slope_shift"]
    if type(j)!=int or j<0:
        raise ValueError("Invalid separating shift.")
    norms,newrows = [],[]
    recipes,cofs = proof["norm_rows"],proof["mark_field_norm_cofactors"]
    if not 1<=len(cofs)==len(recipes)<=2:
        raise ValueError("Incomplete set of generic norm rows.")
    K = chart.q.ring.domain
    for t,co in zip(recipes,cofs):
        if len(t)!=4 or any(type(k)!=int for k in t):
            raise ValueError("Malformed residual combination.")
        entries = {}
        for i,k,a in co:
            if type(i)!=int or type(k)!=int or min(i,k)<0 or (i,k) in entries or len(a)>m.degree():
                raise ValueError("Invalid mark-field cofactor support.")
            entries[i,k] = coordinates(a,K)
        u = chart.q.ring.from_dict(entries)
        r = sum((v*row for v,row in zip(t,chart.rows)),chart.q.ring.zero)
        value = (u*r).rem([chart.q])
        if not value or value.degree(0)!=0:
            raise ArithmeticError("The mark-field cofactor fails its norm identity.")
        e = sp.Poly.from_dict({(v,):a for (i,v),a in value.items()},(S,),domain=K)
        new,n = shift_norm(e,j,m)
        newrows.append(new)
        norms.append(n)
    records = proof["fibres"]
    fs = [read_u(rec["factor"]) for rec in records]
    E = sp.Poly(1,S,domain=sp.QQ)
    if any(f.LC()!=1 or not f.is_irreducible for f in fs) or len(fs)!=len(set(fs)):
        raise ValueError("Invalid rational eliminant factorization.")
    for rec,f in zip(records,fs):
        mult = rec["multiplicity"]
        if type(mult)!=int or mult<1:
            raise ValueError("Invalid component multiplicity.")
        E *= f**mult
    if [p.degree() for p in norms]!=proof["rational_norm_degrees"] or E.degree()!=proof["eliminant_degree"]:
        raise ArithmeticError("Norm degrees disagree.")
    p = proof["coprimality_prime"]
    if type(p)!=int or not sp.isprime(p):
        raise ValueError("Invalid rational coprimality prime.")
    quotients = [n.exquo(E).clear_denoms(convert=True)[1].primitive()[1] for n in norms]
    reduced = [sp.Poly(n.as_expr(),S,modulus=p) for n in quotients]
    if any(a.degree()!=b.degree() for a,b in zip(quotients,reduced)):
        raise ArithmeticError("A norm quotient lost leading degree.")
    g = reduced[0]
    for other in reduced[1:]:
        g = g.gcd(other)
    if g.degree()!=0:
        raise ArithmeticError("The projected univariate cofactors do not match.")
    total = 0
    for record,f in zip(records,fs):
        theta,slope = [read_u(record[k]) for k in ("field_generator","slope")]
        if theta.degree()>=f.degree() or slope!=(sp.Poly(S,S,domain=sp.QQ)-j*theta).rem(f):
            raise ValueError("Invalid recovered embedding or slope.")
        field,spec = fiber_specialization(f,theta,slope)
        field_eq = sp.Poly.from_dict(m.as_dict(),(C,),domain=field)
        ep = sp.Poly.from_dict({(i,):sp.Poly.from_dict(v,(S,),domain=sp.QQ).rem(f).as_expr()
                                for i,v in parts(newrows[0]).items()},(C,),domain=field)
        if field_eq.gcd(ep)!=sp.Poly(C-theta.as_expr(),C,domain=field):
            raise ArithmeticError("The slope fibre does not select the stated embedding.")
        inputs = [spec(x) for x in [chart.q,*chart.rows]]
        def residue(a):
            num = decode(a)
            return sp.Poly.from_dict({(i,):sp.Poly.from_dict(v,(S,),domain=sp.QQ).rem(f).as_expr()
                                      for i,v in parts(num).items()},(C,),domain=field)
        d = residue(record["fibre_gcd"])
        b = [residue(a) for a in record["fibre_cofactors"]]
        if d.is_zero or d.LC()!=1 or len(b)!=len(inputs) or any(p.rem(d) for p in inputs):
            raise ArithmeticError("The supplied monic polynomial does not divide every fibre equation.")
        if d!=sum((u*p for u,p in zip(b,inputs)),d.zero):
            raise ArithmeticError("Incorrect specialized Bezout row.")
        h,good = spec(chart.h),residue(record["c"])
        if good!=d.exquo(d.gcd(h**4)).monic() or good.gcd(good.diff()).degree()!=0 or good.gcd(h).degree()!=0:
            raise ArithmeticError("Incorrect boundary division.")
        if record["parameters"]!=quotient_targets(chart,f,field,good,spec):
            raise ArithmeticError("Incorrect elliptic map data.")
        if f.degree()%m.degree():
            raise ArithmeticError("The field of the mark did not embed.")
        n = f.degree()*good.degree()//m.degree()
        if n!=record["geometric_marks_per_embedding"]:
            raise ArithmeticError("Incorrect branch count per embedding.")
        total += n
    if total!=proof["geometric_marks_per_embedding"]:
        raise ArithmeticError("Incorrect geometric count.")
    return total


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input",type=Path)
    parser.add_argument("--cofactors",action="store_true")
    parser.add_argument("--output",type=Path)
    parser.add_argument("--verify",type=Path,help="replay an exported certificate with its polynomial cofactor rows")
    args = parser.parse_args()
    data = json.loads(args.input.read_text())
    if args.verify:
        verified = verify(data,json.loads(args.verify.read_text()))
        result = json.dumps({"verified_marks_per_embedding":verified})+"\n"
    else:
        result = json.dumps(run(data,args.cofactors),indent=2)+"\n"
    if args.output:args.output.write_text(result)
    else:print(result,end="")
