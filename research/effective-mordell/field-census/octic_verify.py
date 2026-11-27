"""Replay the complete rank-at-most-two panel, the divisor classification and higher-rank witnesses."""
from collections import Counter
from concurrent.futures import ProcessPoolExecutor
from pathlib import Path
import json
from octic_fields import classify,p
from octic_divisors import certify
from higher_witnesses import witnesses

def verify(regenerate=False):
    root=Path(__file__).resolve().parent
    radicands=[d for d in range(2,201) if d%3 and int(p.core(d))==d]
    divisors=certify()
    pairs=Counter(row['rational_class'] for row in divisors)
    if pairs!=Counter({2:12,7:36,11:24,None:24}):raise ArithmeticError("Exceptional divisor census")
    known={str(d):witnesses(d) for d in (7,11)}
    with ProcessPoolExecutor(max_workers=3) as pool:
        panel=list(pool.map(classify,radicands))
    for row in panel:
        d=row['d']
        extra=len([1 for x,v in row.get('twist_points',[]) if v!='0'])
        if row['status']=='partial':
            extra=len(known.get(str(d),[]))
        row['outside_quartic_base']=2*pairs[d]
        row['projective_points_verified']=10+extra+2*pairs[d]
    counts=Counter(row['projective_points_verified'] for row in panel if row['status']=='complete')
    if counts!=Counter({10:48,46:1}):raise ArithmeticError('Complete point counts')
    partial={r['d']:r['projective_points_verified'] for r in panel if r['status']=='partial'}
    if partial[7]!=130 or partial[11]!=58 or len(partial)!=42:raise ArithmeticError('Partial panel')
    for filename,data in [('octic-panel.json',panel),('octic-quadratic-divisors.json',divisors),('higher-octic-witnesses.json',known)]:
        path=root/filename
        normalized=json.loads(json.dumps(data))
        if regenerate:
            path.write_text(json.dumps(data,indent=2)+'\n')
        elif json.loads(path.read_text())!=normalized:
            raise ArithmeticError("Receipt disagreement: "+filename)
    print("49 complete octic fields, 42 partial, whole rational-radicand divisor classification, and 130/58 higher-rank lower bounds passed.")
if __name__=='__main__':
    import argparse
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument('--regenerate',action='store_true')
    verify(cli.parse_args().regenerate)
