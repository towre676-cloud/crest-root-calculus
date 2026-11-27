"""Run the point classifications in this package on one number-field model."""
import argparse
import json
from pathlib import Path
from fixed_curve import classify as baseline
from cypari2 import Pari

p=Pari()
w,x=p('w'),p('x')

def coordinates(archive,nf):
    src=p(archive['field'])
    dst=nf.nf_get_pol()
    if src==dst:return {tuple(pt) for pt in archive.get('affine_points',archive.get('verified_affine_points',[]))}
    if src.poldegree()!=dst.poldegree():raise ArithmeticError('Field degrees disagree.')
    roots=p.nfroots(nf,p.subst(src,w,x))
    if not len(roots):raise ArithmeticError('Expected field embedding missing.')
    root=roots[0]
    def embed(c):return p.subst(p(c),w,root)
    out=set()
    for a,b in archive.get('affine_points',archive.get('verified_affine_points',[])):
        u,v=embed(a),embed(b)
        if v*v!=u**6-1:raise ArithmeticError('Transported point did not verify.')
        out.add((str(p.lift(u)),str(p.lift(v))))
    return out

def octic_panel(nf,subfields):
    """Reconstruct an archived octic list using one consistent embedding in the input model."""
    if not {-1,-3,3}.issubset(subfields) or len(subfields)!=7:
        return None
    radicands=[d for d in subfields if d>1 and d%3]
    if len(radicands)!=1:
        return None
    d=radicands[0]
    def read(name):
        return json.loads(Path(__file__).with_name(name).read_text())
    row=next((r for r in read('octic-panel.json') if r['d']==d),None)
    if row is None:
        return None
    zeta=p.nfroots(nf,p.polcyclo(12,x))[0]
    radical=p.nfroots(nf,x*x-d)[0]
    def embed(value):
        return p.subst(p(value),w,zeta)
    points=set()
    def put(a,b):
        if b*b!=a**6-1:
            raise ArithmeticError('Invalid octic fibre')
        points.add((str(p.lift(a)),str(p.lift(b))))
    for j in range(6):
        put(zeta**(2*j),0)
    for sign in (-1,1):
        put(0,sign*zeta**3)
    twists=row.get('twist_points',[]) if row['status']=='complete' else read('higher-octic-witnesses.json').get(str(d),[])
    for a,v in twists:
        put(embed(a),embed(v)/radical)
    for divisor in read('octic-quadratic-divisors.json'):
        if divisor['rational_class']!=d:
            continue
        U=p.Polrev([embed(a) for a in divisor['u']],x)
        V=p.Polrev([embed(a) for a in divisor['v']],x)
        roots=p.nfroots(nf,U)
        if len(roots)!=2:
            raise ArithmeticError('Certified divisor did not split')
        for t in roots:
            put(1+1/t,p.subst(V,x,t)/t**3)
    if len(points)+2!=row['projective_points_verified']:
        raise ArithmeticError('Octic panel count mismatch')
    return dict(field=str(nf.nf_get_pol()),affine_points=sorted(points),status=row['status'],
                rank=row['rank'])

def classify(co):
    result=baseline(co)
    if result['status']=='complete':return result
    f=p.Polrev(co,w)
    nf=p.nfinit(f)
    degree=len(co)-1
    def merge(record,method):
        prior=set(map(tuple,result['verified_affine_points']))
        current=coordinates(record,nf)
        complete=record['status']=='complete'
        if complete and not prior<=current:raise ArithmeticError('Complete list omitted a verified point.')
        result.update(verified_affine_points=sorted(prior|current),method=method,
                      status='complete' if complete else 'partial')
        return complete
    if degree==8 and len(p.nfroots(nf,p.polcyclo(24,x))):
        archive=json.loads(Path(__file__).with_name('zeta24-points.json').read_text())
        merge(archive,'CM rank-two height bound and saturated cyclotomic lattice')
        return result
    if degree==8:
        proof=octic_panel(nf,result['quadratic_subfields'])
        if proof is not None:
            result['elliptic_rank_interval']=proof['rank']
            merge(proof,'cyclotomic quadratic divisors and elliptic twist fibres')
            return result
    from gaussian_tower_fields import classify as tower
    theorem=tower(co)
    if theorem['status']=='complete':
        merge(theorem,theorem['method'])
        return result
    if degree==4 and -1 in result['quadratic_subfields'] and len(result['quadratic_subfields'])==3:
        d=next(t for t in result['quadratic_subfields'] if t>0)
        from rank_one_height import classify_biquadratic
        if merge(classify_biquadratic(d),'rank-one Gaussian height test and quadratic divisors'):
            return result
        from gaussian_quartic import solve
        eligible=all(q==2 or q%4==3 for q in map(int,p.factorint(d)[0]))
        if eligible and d!=7:
            proof=solve(d)
            result['remaining_gaussian_pairs']=proof['remaining_pairs']
            if merge(proof,'Gaussian squareclass descent and quadratic divisors'):return result
    from classgroup_descent import classify as general
    proof=general(co)
    result['elliptic_rank_upper_bounds']=proof['rank_upper_bounds']
    if proof['status']=='complete':
        merge(proof,'exact number-field 2-isogeny Selmer groups')
    elif degree==4 and -1 in result['quadratic_subfields'] and proof['rank_upper_bounds']['-1']<=1:
        from rank_one_general import points_low_degree
        imag=p.nfroots(nf,x*x+1)[0]
        record=dict(field=str(f),status='complete',affine_points=points_low_degree(nf,imag,degree))
        merge(record,'rank-one Gaussian height bound and exact divisor fibres')
    return result

if __name__=='__main__':
    app=argparse.ArgumentParser(description=__doc__)
    app.add_argument('ascending_polynomial')
    args=app.parse_args()
    print(json.dumps(classify(json.loads(args.ascending_polynomial)),indent=2))
