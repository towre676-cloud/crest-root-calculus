"""Replay the 101..200 extension, including the three-prime sieve for d=118."""
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path

from descent import check, direct_rank, lifts, pari, source, sqrtq
from quadratic_descent import sieve

HERE=Path(__file__).resolve().parent


def reduce_key(P):
    return None if len(P)==1 else tuple(int(t.lift()) for t in P)


def sieve118():
    a4=324*118**2
    coeff=[0,3304,0,a4,0]
    E=pari.ellinit(coeff)
    G=pari(['177/4','114873/8'])
    T=pari([0,0])
    check(E.ellisoncurve(G) and E.elltors()[0]==2 and list(E.ellrank()[:2])==[1,1],
          'rational Mordell-Weil data for d=118')
    check(sqrtq(a4) is not None and sqrtq(Fraction(177,4)) is None,'2-saturation')
    for ell,p,size,gorder in [(3,13,12,12),(5,43,40,20),(7,53,56,56)]:
        F=pari.ellinit(coeff,p)
        check(int(F.ellcard())==size and int(F.ellorder(G))==gorder and
              size%ell==0 and (size//gorder)%ell!=0,'odd saturation witness')
    # N kills all three finite groups and its prime factors have been saturated.
    N=3360
    allowed=set(range(1,48,2))
    records=[]
    for p,card,order in ((5,8,8),(61,60,12),(661,672,48)):
        F=pari.ellinit(coeff,p)
        check(int(E.disc())%p!=0 and int(F.ellcard())==card and N%card==0 and
              int(F.ellorder(G))==order and 48%order==0,'reduction or period')
        source_x={708*pow(t,4,p)%p for t in range(p)}
        candidates={None} # rational source poles reduce to O
        for xx in source_x:
            candidates.update(reduce_key(pari([xx,yy])) for yy in F.ellordinate(xx))
        allowed_mod={n for n in range(1,order,2)
                     if reduce_key(F.elladd(F.ellmul(G,n),T)) in candidates}
        allowed={n for n in allowed if n%order in allowed_mod}
        records.append(dict(p=p,order=order,odd_n=sorted(allowed_mod),
                            remaining_mod_48=sorted(allowed)))
    check(not allowed,'unresolved torsion coset')
    s=pari('Mod(w,w^2+2)')
    root=(43+7*s)/4
    beta=pari('3/16')+7+4*s
    check(root**2==(13-2*s)*beta,'squareclass of the odd-G coset')
    return records


def extended_panel():
    rows={}
    for d in range(101,201):
        if any(d%(i*i)==0 for i in range(2,isqrt(d)+1)):
            continue
        r=direct_rank(d)
        check(r['torsion']==[2,2] and all(l==h for l,h in r['ranks']),'direct rank and torsion')
        rec={'direct':r}
        if all(hi for lo,hi in r['ranks']):
            classes=lifts(d)
            check(all(not c.get('rational_x') for c in classes),'unexpected rank-zero fibre')
            rec['squareclasses']=classes
            rec['quadratic']={str(c['a']):sieve(d,c['a']) for c in classes
                              if c.get('ranks',[0,0])[1]}
        rows[str(d)]=rec
    expected=[(107,2,[-17,4]),(107,2,[-11,10]),(118,1,[7,-8]),(118,1,[13,-2]),
              (131,2,[1,-14]),(131,2,[19,4]),(166,1,[-11,-8])]
    surviving=[(int(d),int(a),z['delta']) for d,rec in rows.items()
               for a,branches in rec.get('quadratic',{}).items()
               for z in branches if z['obstruction'] is None]
    check(surviving==expected,'remaining quadratic branches')
    finite_x={107:Fraction(50),118:Fraction(1,4),131:Fraction(2)}
    for d,rec in rows.items():
        affine=[['0','0']]
        if int(d) in finite_x:
            x=finite_x[int(d)]
            y=sqrtq(source(x,int(d)))
            check(y is not None,'extra rational point')
            affine.extend(([str(x),str(-y)],[str(x),str(y)]))
        rec['affine']=affine
        rec['rational_infinities']=1
    check(len(rows)==61,'second interval count')
    return rows


if __name__=='__main__':
    rec={'mordell_weil_118':sieve118(),'parameters':extended_panel()}
    path=HERE/'extension-receipt.json'
    if path.exists():
        check(json.loads(path.read_text())==rec,'extension receipt')
    else:
        path.write_text(json.dumps(rec,indent=2)+'\n')
    print('Three-prime sieve:',rec['mordell_weil_118'])
    print('61 additional parameters checked. Extra points at 107, 118, 131.')
