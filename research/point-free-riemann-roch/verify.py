"""Replay divisor spaces, independent local and projective controls, and cubic reconstruction."""
import argparse
import json
from pathlib import Path

import sympy as sp

from cyclic_quotients import recover
from embedded_curve import EmbeddedCurve
from fermat_divisors import parameterize
from family_spectrum import verify_family
from riemann_roch import PlaneCurve, T, U, V, X, Y, Z, terms
from symmetric_power import A, B, P, affine_chart, affine_pieces, chart_equations, chords, encode, from_projective, moment_frame

HERE = Path(__file__).resolve().parent


def require(test,message):
    if not test:
        raise ArithmeticError(message)


def rejected(action):
    try:
        action()
    except ValueError:
        return
    raise ArithmeticError('An inadmissible input was accepted')


def is_multiple(curve, function, result):
    """For a one-dimensional space check its function differs by a rational scalar."""
    numerator,denominator = map(sp.expand,sp.fraction(sp.cancel(function)))
    a = sp.reduced(numerator*result['denominator'],[curve.form],*P)[1]
    b = sp.reduced(denominator*result['numerators'][0],[curve.form],*P)[1]
    return result['dimension']==1 and sp.cancel(a/b).is_rational is True


def identify_fermat_classes():
    curve = PlaneCurve(X**4+Y**4+Z**4)
    ds = list(chords(curve,1))
    require(len(ds)==4,'Four rational bitangents expected')
    O = sp.S.One
    coords = {X:U,Y:V,Z:O}
    ideals = [[L.xreplace(coords),q.xreplace(coords)] for _,L,q in ds]
    dims = []
    for m in range(7):
        pos = [(m,ideals[0])] if m else []
        r = curve.space(pos)
        require(r['degree']==2*m,'Nonreduced divisor length differs')
        curve.check_duality(pos,[],r)
        dims.append(r['dimension'])
    require(dims==[1,1,3,4,6,8,10],'Incorrect sequence for a nonrational closed divisor')
    for i in range(1,4):
        require(curve.space([(1,ideals[i])],[(1,ideals[0])])['dimension']==0,'Distinct bitangent pairs coincide in Pic^2')
        r = curve.space([(2,ideals[i])],[(2,ideals[0])])
        require(is_multiple(curve,ds[0][1]/ds[i][1],r),'Divisor-equivalence function was not recovered')
    overlap = curve.space([(3,ideals[1]),(1,ideals[2])],[(2,ideals[1]),(1,ideals[2])])
    require(overlap['dimension']==1 and overlap['degree']==2,'Overlapping positive and negative pieces did not cancel')
    powers = {str(m):[p['projective_dimension'] for p in parameterize(m)] for m in range(9)}
    for n in range(1,5):
        for i in range(4):
            pieces = [(1,ideals[i])]+([(n-1,ideals[0])] if n>1 else [])
            require(curve.space(pieces)['dimension']==powers[str(2*n)][i]+1,'Independent complete-system dimension differs')
    return dims,powers


def local_and_symmetric_controls():
    curve=PlaneCurve(Y*Z**3+Y**4-X**4)
    expected=[1,1,1,2,3,3,4,5,6]
    dims=[]
    for m in range(9):
        pos=[(m,[U,V])] if m else []
        r=curve.space(pos)
        curve.check_duality(pos,[],r)
        dims.append(r['dimension'])
    require(dims==expected,'Hyperflex semigroup differs from <3,4>')
    # At the affine origin y+y^4=x^4. Uniformizer x gives y=x^4 mod x^16.
    for f,order in [(X/Y,3),(Z/Y,4)]:
        a,b=sp.fraction(f.xreplace({X:T,Y:T**4,Z:sp.S.One}))
        require(sp.degree(b,T)-sp.degree(a,T)==order,'Independent leading Laurent order differs')
    conic=PlaneCurve(X*X+Y*Y+Z*Z)
    for m in range(1,5):
        r=conic.space([(m,[U])]);conic.check_duality([(m,[U])],[],r)
        require(r['dimension']==2*m+1,'Conic Riemann--Roch failed')
    # Rotate the line at infinity automatically for divisor Z=0 on a pointless conic.
    data={'terms':terms(conic.form,P),'projective_positive':[{'generators':[terms(Z,P)]}]}
    model,pieces,_=from_projective(data)
    require(model.frame!=sp.eye(3) and model.space(pieces)['dimension']==3,'Projective support was lost')
    rejected(lambda:from_projective({**data,'frame':sp.eye(3).tolist()}))
    s,h,g=encode(conic,[(1,[U])]);require((s,h,g)==(1,T*T+1,T),'Separating-slope fallback failed')
    # Interpolation must also retain a threefold thickening, including its second jet.
    _,h3,g3=encode(curve,[(3,[U,V])]);require((h3,g3)==(T**3,0),'Nonreduced interpolation error')
    for m in (2,3):
        fidx,slope=m,1
        eq=chart_equations(curve,m,fidx,slope)
        a,b=sp.symbols(f'a0:{m}'),sp.symbols(f'b0:{m}')
        sub={q:n-1 for n,q in enumerate((*a,*b))}
        # Cross-check all coefficients by direct ordinary polynomial division after specialization.
        hp=T**m+sum(sub[a[i]]*T**i for i in range(m));gp=sum(sub[b[i]]*T**i for i in range(m))
        coords=moment_frame(fidx).inv()*sp.Matrix([T-slope*gp,gp,1])
        rem=sp.Poly(curve.form.xreplace(dict(zip(P,coords))),T).rem(sp.Poly(hp,T))
        require([sp.sympify(p).subs(sub) for p in eq['equations']]==[rem.nth(i) for i in range(m)],'Symmetric-power equations failed')
    quintic=PlaneCurve(X**5+Y**5+Z**5)
    for n in range(1,5):
        r=quintic.space([(n,[U])]);quintic.check_duality([(n,[U])],[],r)
        require(r['dimension']==(n+1)*(n+2)//2,'Genus-six plane test failed')
    rejected(lambda: PlaneCurve(X*X*Y*Y))
    rejected(lambda: conic.space([(1,[U**2+V**2+1])]))
    return dims


def embedded_controls():
    # A smooth rational quartic in P^3 with one missing ambient section in degree one.
    data=json.loads((HERE/'space-quartic.json').read_text(encoding='utf-8'))
    curve=EmbeddedCurve.from_data(data);a,b,c,d=curve.variables
    r=curve.space(curve.pieces(data['projective_positive']))
    require((curve.degree,curve.genus,r['dimension'],r['denominator_degree'])==(4,0,5,2),'Incomplete linear series was missed')
    require(len(curve.ambient_sections(1))==4 and len(curve.ambient_sections(2))==9,'Incorrect Hilbert function')
    s,t=sp.symbols('s t');param={a:s**4,b:s**3*t,c:s*t**3,d:t**4}
    rational=[sp.cancel(f.xreplace(param)/r['denominator'].xreplace(param)) for f in r['numerators']]
    require(set(rational)=={(t/s)**k for k in range(5)},'Independent rational parametrization check failed')
    require(curve.space([(2,[b,c,d])],[(1,[b,c,d])])['dimension']==2,'Space-curve overlapping divisors failed')
    qs=EmbeddedCurve((a,b,c,d),[a*a+b*b+c*c+d*d,a*a+2*b*b+3*c*c+4*d*d])
    q=qs.space([(1,[a])]);require((qs.genus,q['degree'],q['dimension'])==(1,4,4),'Intersection of quadrics failed')
    # Smooth curve plus an isolated reduced point has nonnegative arithmetic genus here.
    # The tangent upper bound alone cannot remove it; geometric connectedness does.
    F=a**4+b**4+c**4
    rejected(lambda: EmbeddedCurve((a,b,c),[a*F,b*F]))
    # The nonreduced degree-two subscheme at parametric t=0 has coordinates [1,t,0,0] mod t^2.
    eq=curve.symmetric_chart(2,0,0)
    sub={k:0 for row in eq['equations'] for f in row for k in sp.sympify(f).free_symbols}
    require(all(sp.sympify(f).subs(sub)==0 for row in eq['equations'] for f in row),'Embedded symmetric chart lost the double jet')
    line=EmbeddedCurve((a,b),[])
    line_hits=list(line.search(2,1));binary=set()
    for hit in line_hits:
        A,B=hit['coefficient_numerators'];den=hit['denominator'];L=a+hit['chart_index']*b
        binary.add(sp.Poly(den*b*b+B*b*L+A*L*L,a,b,domain=sp.QQ).monic().as_expr())
    require(len(line_hits)==len(binary)==11,'Projective-line atlas search missed or duplicated a divisor')
    conic=EmbeddedCurve((a,b,c),[a*a+b*b+c*c])
    require(not list(conic.search(1,2)) and len(list(conic.search(2,1)))==6,'Anisotropic conic symmetric-power search failed')
    return {'bounded_line_symmetric_square':len(binary),'bounded_conic_points':0,'bounded_conic_symmetric_square':6,
            'space_quartic_degree_one_image':4,'complete_section_dimension':5,'certificate_denominator_degree':2,
            'genus_one_intersection_of_quadrics_dimension':4}


def quotient_controls():
    F=PlaneCurve(X**4+Y**4+Z**4)
    r=recover(F)
    require(len(r['cyclic_cubic_quotients'])==4,'Fermat cyclic quotient count failed')
    quartic=sp.Poly(sp.sympify(r['cyclic_cubic_quotients'][0]['target_quartic']),T)
    a,b,c,d,e=quartic.all_coeffs();I=12*a*e-3*b*d+c*c;J=72*a*c*e+9*b*c*d-27*a*d*d-27*b*b*e-2*c**3
    j=sp.cancel(6912*I**3/(4*I**3-J**2));require(j==287496,'Fermat target invariant differs')
    # Mixed no-real-point quartic. Every deck generator is defined over Q(sqrt(6)).
    N=PlaneCurve.from_data(json.loads((HERE/'nonsplit.json').read_text(encoding='utf-8')))
    q=recover(N);maps=q['cyclic_cubic_quotients']
    require(len(q['rational_bitangents'])==len(maps)==1 and maps[0]['line']==[1,-1,1], 'Non-split positive missed')
    deck=sp.Rational(maps[0]['deck_extension_squareclass']);require(sp.sqrt(6*deck).is_rational,'Wrong quadratic deck field')
    # The contacts are rational on this input, but the deck generator needs Q(sqrt(-3)).
    special=Z**4+Z**2*X*Y+Z*(X**3+2*Y**3)+X**2*Y**2
    special_maps=recover(PlaneCurve(special))['cyclic_cubic_quotients']
    require(len(special_maps)==1 and special_maps[0]['line']==[0,0,1] and special_maps[0]['deck_extension_squareclass']=='-3','Last dual line chart failed')
    for phi,want in [({Y:Z,Z:Y},[0,1,0]),({X:Z,Z:X},[1,0,0])]:
        found=recover(PlaneCurve(special.xreplace(phi)))['cyclic_cubic_quotients']
        require(len(found)==1 and found[0]['line']==want,'Changing projective frames lost a quotient')
    flex=recover(PlaneCurve(Y*Z**3+Y**4-X**4));require(len(flex['rational_bitangents'])==2 and not flex['cyclic_cubic_quotients'],'Hyperflex or noncyclic bitangent falsely accepted')
    zero=recover(PlaneCurve(F.form+X*X*Y*Z+X*Z**3));require(zero['rational_bitangents']==[],'Empty bitangent control failed')
    q0=X*X+2*Y*Y
    miss=recover(PlaneCurve(q0*q0+Z*(X**3-6*X*Y*Y)+3*Z*Z*q0+10*Z**4+Z**3*X))
    require(miss['rational_bitangents']==[[0,0,1]] and not miss['cyclic_cubic_quotients'],
            'A broken cyclic symmetry passed with its ordinary bitangent intact')
    return {'fermat_rational_maps':4,'fermat_target_j':int(j),'mixed_positive':maps[0],
            'rational_contact_control':special_maps[0], 'ordinary_bitangent_near_miss_maps':0, 'hyperflex_curve_bitangents':2,'hyperflex_curve_cyclic_quotients':0}


def run():
    dims,powers=identify_fermat_classes()
    return {'fermat_dimensions_LmD':dims,'fermat_rational_symmetric_powers':powers,
            'hyperflex_dimensions_LmP':local_and_symmetric_controls(),'space_curves':embedded_controls(),
            'cubic_quotients':quotient_controls(), 'arithmetic_family':verify_family()}


if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--json',action='store_true')
    args=parser.parse_args();actual=run()
    if args.json:
        print(json.dumps(actual,indent=2))
    else:
        if actual!=json.loads((HERE/'receipt.json').read_text(encoding='utf-8')):
            raise ArithmeticError('The computed research receipt differs')
        print('Verified point-free divisor spaces, symmetric-power charts and rational cyclic cubic quotients.')
