"""Run the preserved Part XXX decision algorithm on bounded atlas cases."""
from pathlib import Path
import sys

ROOT=Path(__file__).resolve().parent
ENGINE=ROOT/'vendor'/'crest'/'scan_completeness_xxx'/'decision'

def solve(data):
    try:
        import numpy
        import sympy
    except ImportError as error:
        raise ValueError("Part XXX needs NumPy and SymPy. Install the packages in atlas/original/requirements.txt, or select the standard-library case.") from error
    if str(ENGINE) not in sys.path: sys.path.insert(0,str(ENGINE))
    from decide import choose_prime,decide,fa
    from algebra import parse_form
    entries=data.get('curves')
    if not isinstance(entries,list) or not 1<=len(entries)<=32:
        raise ValueError("Supply one to 32 curves")
    primes=data.get('primes',[5,7,11])
    cap=data.get('max_field_size',64)
    depth=data.get('max_precision',4)
    if not isinstance(primes,list) or len(set(primes))!=len(primes) or not 1<=len(primes)<=5 or any(type(p) is not int or p not in [5,7,11,13,17] for p in primes):
        raise ValueError("Choose distinct primes from 5, 7, 11, 13 and 17")
    if type(cap) is not int or not 1<=cap<=256 or type(depth) is not int or not 1<=depth<=12:
        raise ValueError("Field-size budget must be 1–256 and precision 1–12")
    if any(not isinstance(e,dict) or not isinstance(e.get('id'),str) or not e['id'] for e in entries):
        raise ValueError("Give each curve a name")
    if len({e['id'] for e in entries})!=len(entries):
        raise ValueError("Curve names must be unique")
    curves=[]
    for e in entries:
        poly={'terms':[[*m,c] for m,c in parse_form(e['quartic']).items()]}
        F=fa.input_quartic(poly);fa.smooth_check(F)
        row={'id':e['id'],'status':'unresolved','prime_checks':[]}
        for p in primes:
            try: _,_,degree,_,_=choose_prime(F,p,'geometric')
            except ValueError as error:
                if 'singular reduction' not in str(error): raise
                row['prime_checks'].append({'p':p,'smooth':False});continue
            q=p**degree if degree else 1
            row['prime_checks'].append({'p':p,'smooth':True,'scan_degree':degree,'field_size':q})
            if q>cap:continue
            proof=decide(poly,p,'geometric',max_precision=depth)
            if proof['complete']:
                row['proof']=proof
                row['status']='degree_three_map_count_certified' if proof['lower'] else 'no_degree_three_map_over_Qbar'
                break
            row.setdefault('partial_proofs', []).append(proof)
        else:
            partial=row.get('partial_proofs', [])
            if partial:
                row['proof']=min(partial,key=lambda entry:entry['upper']-entry['lower'])
                row['reason']='No tested prime closed the bounds within the precision budget.'
            else:
                row['reason']='No allowed prime supplied a smooth reduction within the field-size budget.'
        curves.append(row)
    return {'engine':'part_xxx','counts':{k:sum(r['status']==k for r in curves) for k in sorted({r['status'] for r in curves})},'curves':curves}
