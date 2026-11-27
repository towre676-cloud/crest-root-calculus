"""Replay the 145 archived higher-degree field models and the cyclotomic cross-check."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
import json
from pathlib import Path
from isogeny_over_field import classify

def run(row):
    polynomial,expected=row
    result=classify(polynomial,512)
    # JSON normalization reflects the stored representation of coordinate pairs.
    if json.loads(json.dumps(result))!=expected:
        raise ArithmeticError(f'Field descent disagreed with receipt for {polynomial}.')
    return (len(polynomial)-1,result['status'])

if __name__=='__main__':
    rows=json.loads(Path(__file__).with_name('higher-field-panel.json').read_text())
    with ProcessPoolExecutor(max_workers=3) as executor:
        counts=Counter(executor.map(run,rows))
    if counts!={(3,'complete'):19,(3,'partial'):56,(4,'complete'):2,(4,'partial'):33,
               (5,'complete'):2,(5,'partial'):33}:
        raise ArithmeticError('Unexpected field panel totals.')
    from fixed_curve import classify as classify_field
    certificate=classify([1,0,-1,0,1],512)
    tower=classify_field([1,0,-1,0,1])
    if (certificate['status']!='complete' or len(certificate['verified_affine_points'])!=8 or
        set(map(tuple,certificate['verified_affine_points']))!=set(map(tuple,tower['verified_affine_points']))):
        raise ArithmeticError('Twelfth cyclotomic field classification.')
    print('Verified',dict(sorted(counts.items())),'and both proofs of the ten cyclotomic points.')
