"""Replay the quadratic field classification and check the archived point table."""
import argparse
import json
from collections import Counter
from pathlib import Path
from fixed_curve import bounded_quadratic_abscissas, classify, signed_squarefree

def panel():
    rows=[]
    for d in range(-200,201):
        if d in (0,1) or signed_squarefree(d)!=d:
            continue
        result=classify([-d,0,1])
        known_x={x for x,y in result['verified_affine_points']}
        for x in bounded_quadratic_abscissas().get(d,[]):
            if str(x) not in known_x or str(-x) not in known_x:
                raise ArithmeticError(f'missing bounded rational abscissa in field {d}')
        rows.append(dict(d=d,status=result['status'],method=result['method'],
                         verified_points=2+len(result['verified_affine_points']),
                         remaining_squareclass_pairs=result.get('remaining_squareclass_pairs',[])))
    return rows

if __name__=='__main__':
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--write',action='store_true',help='Regenerate the receipt.')
    opt=parser.parse_args()
    rows=json.loads(json.dumps(panel()))
    target=Path(__file__).with_name('quadratic-field-panel.json')
    if opt.write:
        target.write_text(json.dumps(rows,indent=2)+'\n')
    elif json.loads(target.read_text())!=rows:
        raise ArithmeticError('Quadratic panel disagrees with the receipt.')
    print(len(rows),'quadratic fields:',dict(Counter(r['status'] for r in rows)))
