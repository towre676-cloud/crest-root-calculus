"""Replay both signs for every squarefree radicand from 201 to 1000."""
import argparse
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from sympy import mobius
from fixed_curve import classify


def classify_signs(d):
    rows=[]
    for t in (d,-d):
        result=classify([-t,0,1])
        rows.append(dict(d=t,status=result['status'],method=result['method'],
                         point_count=2+len(result['verified_affine_points']),
                         remaining_squareclass_pairs=result.get('remaining_squareclass_pairs',[])))
    return rows


if __name__=='__main__':
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--regenerate',action='store_true')
    args=cli.parse_args()
    numbers=[d for d in range(201,1001) if mobius(d)]
    with ProcessPoolExecutor(max_workers=3) as pool:
        rows=[row for pair in pool.map(classify_signs,numbers) for row in pair]
    target=Path(__file__).with_name('extended-quadratic.json')
    if args.regenerate:
        target.write_text(json.dumps(rows,indent=2)+'\n')
    elif json.loads(target.read_text())!=json.loads(json.dumps(rows)):
        raise ArithmeticError('Extended signed quadratic census disagreement')
    print(len(rows),dict(Counter(row['status'] for row in rows)))
    print('Counts in completed fields',dict(Counter(row['point_count'] for row in rows if row['status']=='complete')))
