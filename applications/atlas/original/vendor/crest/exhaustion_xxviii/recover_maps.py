"""Recover certified cubic nets from free_atlas RUR files and group Galois orbits.

A replay of supplied files checks the individual points. For a completeness
statement use --solve-all, which computes all twelve unrestricted charts.
"""
import argparse
import ast
import json
from pathlib import Path
import sympy as sp
import cubic_net as cn
import free_atlas as fa

T=sp.Symbol('T')
QCOORDS=sp.symbols('A B C')


def parse_poly(text):
    """Parse only rational polynomial arithmetic in the RUR parameter."""
    def convert(n):
        if isinstance(n,ast.Constant) and type(n.value)==int:
            return sp.Integer(n.value)
        if isinstance(n,ast.Name) and n.id=='T':
            return T
        if isinstance(n,ast.UnaryOp) and isinstance(n.op,(ast.UAdd,ast.USub)):
            v=convert(n.operand)
            return -v if isinstance(n.op,ast.USub) else v
        if isinstance(n,ast.BinOp):
            a,b=convert(n.left),convert(n.right)
            if isinstance(n.op,ast.Add): return a+b
            if isinstance(n.op,ast.Sub): return a-b
            if isinstance(n.op,ast.Mult): return a*b
            if isinstance(n.op,ast.Div) and b.is_Rational and b!=0: return a/b
            if isinstance(n.op,ast.Pow) and b.is_Integer and b>=0: return a**b
        raise ValueError("Expected a rational polynomial in T")
    return sp.Poly(convert(ast.parse(str(text),mode='eval').body),T,domain=sp.QQ)


def field_for(P):
    if P.degree()==1:
        return sp.QQ, -P.nth(0)/P.nth(1)
    if P.degree()==2:
        a,b,c=P.all_coeffs()
        alpha=(-b+sp.sqrt(b*b-4*a*c))/(2*a)
    else:
        alpha=sp.CRootOf(P.as_expr(),0)
    return sp.QQ.algebraic_field(alpha),alpha


def coefficient_as_polynomial(coefficient,domain):
    if domain==sp.QQ:
        return domain.to_sympy(coefficient)
    # The ANP coefficients are in descending powers of the specified field root.
    return sp.Poly.from_list(coefficient.to_list() or [0],T,domain=sp.QQ).as_expr()


def orbit_chow(line,domain,minimal):
    """Squarefree norm of the normalized hyperplane in dual coefficient variables."""
    first=next(a for a in line if a)
    linear=sum(coefficient_as_polynomial(a/first,domain)*c for a,c in zip(line,QCOORDS))
    norm=sp.resultant(minimal.as_expr(),linear,T)
    return sp.Poly(norm,*QCOORDS,domain=sp.QQ).sqf_part().monic()


def printable(poly,domain):
    return str(sum(coefficient_as_polynomial(co,domain)*sp.prod(v**a for v,a in zip(poly.gens,m))
                   for m,co in poly.rep.to_dict().items()))


def net_for_solution(frame,kind,point):
    i,j,k=frame
    n=cn.xyz[k]
    lin=cn.xyz[i]-point['u']*n
    tan=cn.xyz[j]-point['v']*n-point['d']*lin
    if kind=='I':
        den=lin+point['s']*tan
        tnum=tan
        unum=-tan*(n+point['b']*den+point['c']*tan)
    elif kind=='II':
        den=tan
        tnum=lin
        unum=-tan*n-point['a']*lin**2-point['b']*lin*tan
    else:
        raise ValueError('Expected chart I or II')
    wnum=(tnum**3+point['e2']*tnum**2*den+point['e1']*tnum*den**2
          +point['f1']*unum*tnum+point['e0']*den**3+point['f0']*unum*den)
    return (unum*den,wnum,den**3)


def extract_chart(F,frame,kind,rur,fixed=None):
    if frame not in fa.FRAMES or kind not in fa.KINDS:
        raise ValueError('Unrecognized frame or chart')
    fixed=fixed or {}
    allvars,ideal,shape_target=fa.chart(F,frame,kind)
    names={str(v):v for v in allvars}
    if rur.get('dimension')==-1 and rur.get('points')==0 and rur.get('length')==0:
        return []
    if rur.get('dimension')!=0:
        raise ValueError('Expected a finite chart RUR')
    P=parse_poly(rur['minimal_polynomial']).monic()
    if P.degree()<=0 or P.gcd(P.diff()).degree() or rur['points']!=P.degree() or rur['length']<P.degree():
        raise ValueError('Inconsistent squarefree support polynomial or counts')
    coord={str(k):parse_poly(v) for k,v in rur['coordinates'].items()}
    for name,value in fixed.items():
        if name in coord: raise ValueError('Fixed and recovered coordinate repeated')
        coord[name]=parse_poly(value)
    if set(names)!=set(coord):
        raise ValueError('Incomplete or extraneous chart coordinates')
    # Check the RUR directly in QQ[T]/(P), including every full chart equation.
    subs={v:coord[name].as_expr() for name,v in names.items()}
    for eq in ideal:
        if not sp.Poly(eq.xreplace(subs),T,domain=sp.QQ).rem(P).is_zero:
            raise ValueError('RUR point fails a chart equation')
    quartic=F.as_expr().xreplace(dict(zip(fa.COORDS,cn.xyz)))
    result=[]
    for factor,_ in P.factor_list()[1]:
        M=factor.monic()
        domain,alpha=field_for(M)
        point={name:poly.rem(M).eval(alpha) for name,poly in coord.items()}
        cubics=net_for_solution(frame,kind,point)
        verified=cn.certify(quartic,cubics,domain,raw=True)
        if verified['map_degree']!=3:
            raise ValueError('A recovered net has incorrect map degree')
        line=verified['differential_line']
        chow=orbit_chow(line,domain,M)
        inverse_g=sp.invert(coord['g'].rem(M),M)
        kappa=[]
        target_coefficients=sp.Poly(shape_target,fa.z)
        for exponent in range(4):
            coeff=target_coefficients.nth(exponent)
            reduced=sp.Poly(coeff.xreplace(subs),T,domain=sp.QQ).rem(M)
            kappa.append((reduced*inverse_g).rem(M).as_expr())
        target=cn.v**2*cn.w-sum(k*cn.u**j*cn.w**(3-j) for j,k in enumerate(kappa))
        result.append({'frame':''.join(map(str,frame)),'kind':kind,'fixed':fixed,
                       'point_polynomial':str(M.as_expr()),'rur_point_degree':M.degree(),
                       'net_in_Q_T_mod_point_polynomial':[printable(sp.Poly(h,*cn.xyz,domain=domain),domain) for h in cubics],
                       'image_cubic_in_Q_T_mod_point_polynomial':str(target),
                       'base_length':verified['base_length'],'map_degree':3,
                       'differential_line':[str(coefficient_as_polynomial(a,domain)) for a in line],
                       'chow_polynomial':str(chow.as_expr()),'orbit_size':chow.total_degree()})
    return result


def summarize(representatives):
    orbits={}
    for item in representatives:
        label=item['chow_polynomial']
        if label not in orbits: orbits[label]={'geometric_size':item['orbit_size'],'chart_representatives':[]}
        orbits[label]['chart_representatives'].append(item)
    return {'galois_orbits':len(orbits),'geometric_classes_found':sum(g['geometric_size'] for g in orbits.values()),
            'orbits':orbits}


def solve_all(F):
    fa.smooth_check(F)
    representatives=[]
    details=[]
    for frame in fa.FRAMES:
        for kind in fa.KINDS:
            variables,ideal,_=fa.chart(F,frame,kind)
            rur=fa.geometric_rur(ideal,variables)
            representatives.extend(extract_chart(F,frame,kind,rur))
            details.append({'frame':''.join(map(str,frame)),'kind':kind,'rur':rur})
    summary=summarize(representatives)
    return {'scope':'all twelve unrestricted charts solved in this run','complete':True,
            'geometric_classes':summary['geometric_classes_found'],**summary,'charts':details}


def read_supplied(F,bundles):
    reps=[]
    for report in bundles:
        frame=tuple(map(int,report['frame']))
        reps.extend(extract_chart(F,frame,report['kind'],report['rur'],report.get('fixed')))
    return {'scope':'supplied chart support only','complete':False,**summarize(reps)}


if __name__=='__main__':
    ap=argparse.ArgumentParser(description=__doc__)
    ap.add_argument('quartic',type=Path)
    group=ap.add_mutually_exclusive_group(required=True)
    group.add_argument('--reports',type=Path,help='JSON list of {frame,kind,rur,fixed?} records')
    group.add_argument('--rur-directory',type=Path,help='directory of unrestricted chart_012_I.json etc from free_atlas solve')
    group.add_argument('--solve-all',action='store_true',help='run every unrestricted solve; may be expensive')
    ap.add_argument('--output',type=Path)
    args=ap.parse_args()
    F=fa.input_quartic(json.loads(args.quartic.read_text()))
    if args.solve_all: result=solve_all(F)
    else:
        if args.reports:
            reports=json.loads(args.reports.read_text())
        else:
            reports=[]
            for frame_tuple in fa.FRAMES:
                frame=''.join(map(str,frame_tuple))
                for kind in fa.KINDS:
                    file=args.rur_directory/f'chart_{frame}_{kind}.json'
                    if file.is_file():
                        reports.append({'frame':frame,'kind':kind,'rur':json.loads(file.read_text())})
            if not reports: ap.error('No unrestricted chart RUR files found')
        result=read_supplied(F,reports)
    text=json.dumps(result,indent=2)+'\n'
    if args.output: args.output.write_text(text)
    else: print(text,end='')
