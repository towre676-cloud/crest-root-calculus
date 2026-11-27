"""Run from the extracted quartic reconstruction source directory."""
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json,os,sys
review=Path(__file__).resolve().parent
sys.path.insert(0,str(Path.cwd()))
sys.path.insert(0,str(review))
def row_check(i):
    from verify_number_fields import polynomial, recover, summary, pari, sp, x
    from verify_solver_records import check_maps, check_orbits, poly_in_x
    r=json.loads(Path('receipt.json').read_text())['sweep'][i]
    expected=json.loads(Path('number-field-receipt.json').read_text())['sweep'][i]
    data=summary(*recover(polynomial(r['coefficients'])))
    rational=sorted(m['half_poles'] for m in data['maps'] if m['point_field']=='Q')
    old=sorted(str(pari(str(sp.Poly.from_list([sp.Rational(c) for c in m['half_poles'][::-1]],x).as_expr()).replace('**','^'))) for m in r['maps'])
    assert rational==old and data['geometric_maps']==expected['geometric_maps']
    assert [m['point_field'] for m in data['maps']]==expected['point_fields']
    assert r['pairing']==expected['pairing'] and r['source']==expected['source']
    from verify_solver_records import x as px
    a=data['reciprocal_chart_center']
    f=poly_in_x(r['coefficients']) if a is None else sum(c*(a*px+1)**j*px**(6-j) for j,c in enumerate(r['coefficients']))
    check_orbits(data,f)
    check_maps(data,f)
    print('verified',i,data['geometric_maps'],flush=True)
    return (i,data['geometric_maps'])
if __name__=='__main__':
    with ProcessPoolExecutor(max_workers=int(os.environ.get("SWEEP_WORKERS", "4"))) as executor:
        results=list(executor.map(row_check,[i for i in range(36) ]))
    print('verified all sweeps',results,flush=True)
