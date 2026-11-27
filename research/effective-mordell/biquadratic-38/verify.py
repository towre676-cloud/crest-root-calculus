"""Replay the field, squareclass, Mordell-Weil, local analytic, and point certificates."""
import json
from pathlib import Path
import subprocess
import sys


def main():
    if not __debug__:
        raise RuntimeError('Use Python with assertions enabled.')
    root=Path(__file__).resolve().parent
    for part in ('prior/quadratic_pairs.py','prior/gaussian_descent.py','upstream/verify_squareclasses.py','upstream/verify_descent.py'):
        subprocess.run([sys.executable,str(root/part)],cwd=root,check=True)
    from biquadratic37 import verify
    if verify()!=json.loads((root/'localization-certificate.json').read_text()):
        raise ArithmeticError('Local certificate mismatch')
    from points import classify,DEFAULT,p,x
    for twist in (-1,1):
        assert list(p.ellinit([0,twist]).ellrank()[:2])==[0,0]
    field=p.nfinit('w^2+1');e=p.ellinit([0,-1],field)
    order,orders,gs=e.elltors()
    tors=[p([0])]
    for n,g in zip(orders,gs):
        tors=[e.elladd(t,e.ellmul(g,j)) for t in tors for j in range(int(n))]
    assert len({str(t) for t in tors})==int(order)
    xs={str(p.lift(a)) for t in tors if len(t)>1 for a in p.nfroots(field,x*x-t[0])}
    assert xs=={'-1','0','1'}
    pts=classify(DEFAULT)
    assert json.loads(json.dumps(pts))==json.loads((root/'biquadratic-38-points.json').read_text())
    assert classify(list(map(int,p.polredabs(p.Polrev(DEFAULT)).Vecrev())))['projective_points']==38
    from higher_genus import classify_covers
    assert json.loads(json.dumps(classify_covers()))==json.loads((root/'higher-genus-points.json').read_text())
    print('All certificates passed. Exactly 38 projective points.')
if __name__=='__main__':
    main()
