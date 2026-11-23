"""Replay a submitted cofactor certificate by polynomial identities.

Input cofactors are never taken as new generators. Verification begins at the
six square-completion equations of the given smooth marked plane quartic.
"""
from math import comb
import sympy as sp
from sympy.polys.rings import ring
from chart import R, S, C, c, s, POSITIONS, adapted_frame, build_from_frame
from solve import ResidueField, coefficients, lifted, quotient_targets, specialize, sparse


def decode(data):
    cols = {}
    for i,j,a in data:
        if type(i)!=int or type(j)!=int or min(i,j)<0 or (i,j) in cols:
            raise ValueError("Invalid cofactor support.")
        val = sp.QQ.convert(sp.Rational(a))
        if not val:
            raise ValueError("Omit zero terms.")
        cols[i,j] = val
    return R.from_dict(cols)


def read_u(vector):
    return sp.Poly.from_dict({(i,):sp.Rational(a) for i,a in enumerate(vector)},(S,),domain=sp.QQ)


def embed(f):
    return R.from_dict({(0,j):a for (j,),a in f.rep.to_dict().items()})


def crt_targets(fibres):
    """The two polynomial generators for the geometric support on h!=0."""
    good = [(read_u(f["s_factor"]),decode(f["unit_fibre"])) for f in fibres if f["geometric_marks"]]
    P = sp.Poly(1,S,domain=sp.QQ)
    for f,g in good:
        P *= f
    if not good:
        return [("unit",R.one)]
    eq = R.zero
    for f,g in good:
        u = P.exquo(f)
        eq += embed(u*u.invert(f))*g
    return [("s_support",embed(P)),("c_support",eq.rem([embed(P)]))]


def original_rows(chart,source=None):
    """Expand the sextic literally in the adapted plane, before imposing q or L."""
    A,b,c3,s3 = ring("b,c,s",sp.QQ)
    cast = lambda p:A.from_dict({(0,i,j):co for (i,j),co in p.items()})
    t = {(i,j):cast(value) for (i,j),value in chart.table.items()}
    def g(i,l):
        weight = 6-i-2*l
        return sum((t[j,l+weight-j]*comb(l+weight-j,l)*b**(weight-j)
                    for j in range(weight+1) if j<=4 and l+weight-j<=3),A.zero)
    if source is not None:
        full,t0,z0,b0,c0,s0 = ring("t,z,b,c,s",sp.QQ)
        substitution = (s0*t0**2-t0,-t0**2,z0+b0*t0+c0*t0**2)
        pullback = sum((co*substitution[0]**i*substitution[1]**j*substitution[2]**k
                        for (i,j,k),co in source.rep.to_dict().items()),full.zero).exquo(t0**2)
        pieces = {}
        for (i,j,k,l,m),co in pullback.items():
            pieces.setdefault((i,j),{})[k,l,m] = co
        if any(i>6 or j>3 for i,j in pieces) or any(
                A.from_dict(pieces.get((i,j),{}))!=g(i,j) for i in range(7) for j in range(4)):
            raise ArithmeticError("The literal projective pullback disagrees with the chart sextic.")
    h,a = g(6,0),g(5,0)
    # Independent coefficient comparison using S^2. The numerator of S is
    # 16 h^3 S. Thus [t^2] of 256 h^5 (G-h S²) is 4h² N2,
    # and its [t] coefficient is 4h N1.
    square = {(3,0):16*h**3,(2,0):8*h*h*a}
    for ell in (0,1):
        square[1,ell] = 2*h*(4*h*g(4,ell)-(0 if ell else a*a))
        square[0,ell] = 8*h*h*g(3,ell)-a*square[1,ell].exquo(2*h)
    rows = {}
    for i in (1,2):
        for l in range(3):
            sq = sum((x*y for (k,j),x in square.items() for (u,v),y in square.items()
                      if k+u==i and j+v==l),A.zero)
            rows[i,l] = (256*h**5*g(i,l)-sq).exquo(4*h**i)
    uq = chart.table[0,1]
    w = chart.table[0,3]
    Q = 3*w*w*chart.q
    if rows[2,2]!=16*h*h*cast(Q) or rows[1,2]!=16*h*h*(cast(chart.denominator)*b-cast(chart.numerator)-b*cast(uq*Q)):
        raise ArithmeticError("The square-coefficient and osculation cofactors disagree.")
    return A,b,cast,rows


def check_raw_reduction(chart,source):
    """Check the four equalities r_j = d^k N_j(n/d) modulo q independently."""
    _,_,_,raw = original_rows(chart,source)
    q = chart.q
    red = lambda p:p.rem([q])
    n,d = red(chart.numerator),red(chart.denominator)
    powers = lambda x:[red(x**i) for i in range(6)]
    np,dp = powers(n),powers(d)
    for (i,l),r in zip(POSITIONS,chart.rows):
        k = 6-i-2*l
        parts = {}
        for (e,j,u),v in raw[i,l].items():
            parts.setdefault(e,{})[j,u] = v
        if any(e>k for e in parts):
            raise ArithmeticError("A weighted residual exceeded its predicted degree in b.")
        observed = sum((red(R.from_dict(p)*red(np[e]*dp[k-e])) for e,p in parts.items()),R.zero)
        if observed!=r:
            raise ArithmeticError("Cofactor substitution failed on an original generator.")


def verify(data,point,result,raw=False):
    source,frame = adapted_frame(data,point)
    chart = build_from_frame(source,frame,R)
    if raw:
        check_raw_reduction(chart,source)
    if chart.frame != result["frame"]:
        raise ValueError("Wrong normalized frame.")
    factors = [read_u(f["s_factor"]) for f in result["fibres"]]
    if any(f.LC()!=1 or not f.is_irreducible for f in factors) or len(set(factors))!=len(factors):
        raise ValueError("Invalid factorization over QQ.")
    G = sp.Poly(1,S,domain=sp.QQ)
    for f,rec in zip(factors,result["fibres"]):
        m = rec["s_multiplicity_in_eliminant"]
        if type(m)!=int or m<=0:
            raise ValueError("Invalid multiplicity.")
        G *= f**m
    p = result["coprimality_prime"]
    if type(p)!=int or not sp.isprime(p):
        raise ValueError("Invalid prime.")
    recipes,us = result["norm_rows"],result["norm_cofactors"]
    if not 1<=len(us)==len(recipes)<=2:
        raise ValueError("Missing generic Bezout rows.")
    gcdp = None
    degrees = []
    for row,co in zip(recipes,us):
        if len(row)!=4 or any(type(v)!=int for v in row):
            raise ValueError("Invalid residual combination.")
        observed = (sum((k*r for k,r in zip(row,chart.rows)),R.zero)*decode(co)).rem([chart.q])
        if not observed or observed.degree(0)>0:
            raise ArithmeticError("The proposed cofactor fails to produce a nonzero eliminant.")
        e = sp.Poly.from_dict({(j,):a for (i,j),a in observed.items()},(S,),domain=sp.QQ)
        degrees.append(e.degree())
        primitive = e.exquo(G).clear_denoms(convert=True)[1].primitive()[1]
        reduced = sp.Poly(primitive.as_expr(),S,modulus=p)
        if reduced.degree()!=primitive.degree():
            raise ArithmeticError("A polynomial lost leading degree at the witness prime.")
        gcdp = reduced if gcdp is None else gcdp.gcd(reduced)
    if gcdp.degree()!=0 or degrees!=result["norm_degrees"] or G.degree()!=result["gcd_degree"]:
        raise ArithmeticError("Elimination certificate failed.")
    bounds = []
    if len(factors)!=len(result["fibre_cofactors"]):
        raise ValueError("Wrong number of local cofactor rows.")
    for f,record,cofs in zip(factors,result["fibres"],result["fibre_cofactors"]):
        field = ResidueField(f)
        spec = lambda x:specialize(x,f,field)
        d = spec(decode(record["fibre_gcd"]))
        good = spec(decode(record["unit_fibre"]))
        q,h = spec(chart.q),spec(chart.h)
        inputs = [q,*map(spec,chart.rows)]
        if d.is_zero or d.LC()!=1 or len(cofs)!=len(inputs):
            raise ValueError("Malformed monic fibre polynomial.")
        if any(poly.rem(d) for poly in inputs) or d!=sum((spec(decode(a))*v for a,v in zip(cofs,inputs)),d.zero):
            raise ArithmeticError("The local Bezout row fails.")
        if good!=d.exquo(d.gcd(h**4)).monic() or good.gcd(good.diff()).degree()!=0 or good.gcd(h).degree()!=0:
            raise ArithmeticError("Boundary or multiplicity error.")
        if f.degree()*good.degree()!=record["geometric_marks"]:
            raise ArithmeticError("Incorrect geometric count.")
        if quotient_targets(chart,f,field,good)!=record["parameters"]:
            raise ArithmeticError("Incorrect map formulas.")
        bounds.append((f,field,d,h,record))
    targets = crt_targets(result["fibres"])
    for _,target in targets:
        for f,K,d,h,record in bounds:
            value = specialize(target,f,K) if record["geometric_marks"] else d.one
            if (h**4*value).rem(d):
                raise ArithmeticError("The claimed global target has no local cofactor.")
    if result["geometric_marks"]!=sum(r["geometric_marks"] for r in result["fibres"]):
        raise ArithmeticError("Incorrect geometric total.")
    length = sum(f["s_multiplicity_in_eliminant"] for f in result["fibres"])
    good_exp = sum(f["s_multiplicity_in_eliminant"] for f in result["fibres"] if f["geometric_marks"])
    return {"geometric_marks":result["geometric_marks"],"cofactor_targets":{k:sparse(v) for k,v in targets},
            "h_exponent":4*length+2,"target_exponent":max(1,good_exp)}
