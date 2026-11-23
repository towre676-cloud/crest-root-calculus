"""Certify the degree and pullback differential of a plane cubic net.

Inputs are exact rational ternary forms, or forms over an explicit SymPy
algebraic field. The underlying curve must be a smooth plane quartic.
"""
import json
import sympy as sp
from sympy.polys.matrices import DomainMatrix

x,y,z=sp.symbols('x y z')
u,v,w=sp.symbols('U V W')
xyz=(x,y,z)


def terms(degree):
    return [(i,j,degree-i-j) for i in range(degree,-1,-1)
            for j in range(degree-i,-1,-1)]


def monomial(coords,exps):
    return sp.prod(a**b for a,b in zip(coords,exps))


def polynomial(expr,coords,domain):
    return sp.Poly(expr,*coords,domain=domain)


def form(expr,coords,degree,domain):
    F=polynomial(expr,coords,domain)
    if F.is_zero or any(sum(m)!=degree for m in F.monoms()):
        raise ValueError(f'Expected a nonzero homogeneous form of degree {degree}.')
    return F


def columns_matrix(columns,monoms,domain):
    dictionaries=[col.rep.to_dict() for col in columns]
    rows=[[col.get(m,domain.zero) for col in dictionaries] for m in monoms]
    return DomainMatrix(rows,(len(rows),len(columns)),domain)


def null_vectors(matrix):
    if matrix.domain.is_AlgebraicField:
        reduced,pivots=matrix.rref(method='GJ')
        return reduced.nullspace_from_rref(pivots).to_list()
    return matrix.nullspace().to_list()


def normal_form(poly,divisor):
    """Full lexicographic division by one polynomial, in its native field."""
    base=divisor.rep.to_dict()
    lead=max(base)
    lc=base[lead]
    zero=divisor.domain.zero
    pending=poly.rep.to_dict().copy()
    remainder={}
    while pending:
        m=max(pending)
        if any(i<j for i,j in zip(m,lead)):
            remainder[m]=pending.pop(m)
            continue
        shift=tuple(i-j for i,j in zip(m,lead))
        scalar=pending[m]/lc
        for n,c in base.items():
            key=tuple(i+j for i,j in zip(n,shift))
            value=pending.get(key,zero)-scalar*c
            if value:pending[key]=value
            else:pending.pop(key,None)
    return sp.Poly.from_dict(remainder,poly.gens,domain=poly.domain)


def smooth(F):
    for i,coord in enumerate(F.gens):
        others=F.gens[:i]+F.gens[i+1:]
        equations=[F.diff(a).as_expr().subs(coord,1) for a in F.gens]
        if list(sp.groebner(equations,*others,domain=F.domain))!=[1]:
            return False
    return True


def base_length(F,net):
    cols=[]
    for source in [F,*net]:
        for m in terms(8-source.total_degree()):
            cols.append(source*polynomial(monomial(xyz,m),xyz,F.domain))
    matrix=columns_matrix(cols,terms(8),F.domain)
    rank=matrix.rank()
    return 45-rank,rank


def relation(F,net):
    keys=terms(3)
    col=[normal_form(monomial(net,m),F) for m in keys]
    null=null_vectors(columns_matrix(col,terms(9),F.domain))
    if len(null)!=1:
        raise ValueError('The cubic relation space must have dimension one.')
    lead=next(c for c in null[0] if c)
    eq=sum(F.domain.to_sympy(c/lead)*monomial((u,v,w),m)
           for c,m in zip(null[0],keys))
    Q=form(eq,(u,v,w),3,F.domain)
    if not smooth(Q):
        raise ValueError('The image cubic is singular.')
    return Q


def differential(F,net,Q,raw=False):
    domain=F.domain
    affine=lambda e:polynomial(e.subs(z,1), (x,y),domain)
    source=affine(F.as_expr())
    A,B,C=[affine(g.as_expr()) for g in net]
    deriv=(C*A.diff(x)-A*C.diff(x))*source.diff(y)-(C*A.diff(y)-A*C.diff(y))*source.diff(x)
    def reduced(poly):
        return normal_form(poly,source)
    deriv=reduced(deriv)
    denom=sum((monomial((A,B,C),m).mul_ground(co) for m,co in Q.diff(v).rep.to_dict().items()),
              polynomial(0,(x,y),domain))
    cols=[reduced(denom*affine(a)) for a in xyz]
    assert not reduced(denom).is_zero
    monoms=sorted(set().union(*(col.monoms() for col in [*cols,deriv])))
    kernel=null_vectors(columns_matrix([*cols,-deriv],monoms,domain))
    if len(kernel)!=1 or not kernel[0][-1]:
        raise ValueError('Pullback differential identity failed.')
    coef=[c/kernel[0][-1] for c in kernel[0][:3]]
    first=next(v for v in coef if v)
    normalized=[v/first for v in coef]
    return normalized if raw else [str(domain.to_sympy(v)) for v in normalized]


def certify(quartic,cubics,domain=sp.QQ,raw=False):
    F=form(quartic,xyz,4,domain)
    if not smooth(F):
        raise ValueError('Singular source quartic.')
    if len(cubics)!=3:
        raise ValueError('Provide exactly three cubics.')
    net=[form(f,xyz,3,domain) for f in cubics]
    if columns_matrix(net,terms(3),domain).rank()!=3:
        raise ValueError('Linearly dependent net.')
    Q=relation(F,net)
    length,rank=base_length(F,net)
    degree,extra=divmod(12-length,3)
    if extra or degree<1:
        raise ArithmeticError('Inconsistent moving degree.')
    return {'base_length':length,'degree_eight_rank':rank,'map_degree':degree,
            'image_cubic':str(Q.as_expr()),'differential_line':differential(F,net,Q,raw=raw)}


def run():
    # A tangent-chart map in which all three fixed base points coincide.
    F=x**4+3*x**2*y*z+x*y**3-8*x*z**3-3*y*y*z*z
    conic=x*z+y*y/2
    net=(-x*conic, y**3-6*y*conic-x**3,x**3)
    triple=certify(F,net)
    assert triple['map_degree']==3 and triple['base_length']==3
    # Target PGL_3 leaves both the differential line and map degree unchanged.
    a,b,c=net
    alt=certify(F,(a+c,a+b,b+c))
    assert alt['map_degree']==3 and alt['differential_line']==triple['differential_line']
    # Pullback to an elliptic quotient of degree two with a common line in the net.
    F2=x**4-y*z*(y-z)*(y-3*z)
    double=certify(F2,(y*z*z,x*x*z,z**3))
    assert double['map_degree']==2 and double['base_length']==6
    # Smooth cubic source-map hypothesis must reject a false candidate.
    try:
        certify(F,(a+b+x*y*z,b,c))
    except ValueError as err:
        rejected=str(err)
    else:
        raise AssertionError('Perturbed cubic net was accepted.')
    return {'tangent_map':triple,'changed_target_coordinates':alt,
            'bielliptic_control':double,'perturbed_net_rejected':rejected}


if __name__=='__main__':
    print(json.dumps(run(),indent=2))
