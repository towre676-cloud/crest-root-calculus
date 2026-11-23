"""A separately seeded test panel for constant-vertical hyperelliptic maps."""
import hashlib
import json
from pathlib import Path
import random
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'bounded-correspondences'))
from maps import curve, X
from finite import coefficients


def create():
    rng = random.Random(2026112319)
    rows = []
    for genus, degree in [(1,1),(1,3),(1,5),(2,3),(3,5),(4,7)]:
        while True:
            G = X**(2*genus+1) + sum(rng.randint(-4,4)*X**j for j in range(2*genus+1))
            H = sp.Rational(rng.choice([2,3,4]),rng.choice([1,3]))*X**degree
            H += sum(sp.Rational(rng.randint(-3,3),rng.choice([1,2]))*X**j for j in range(degree))
            vertical = sp.Rational(3,2)
            F = sp.Poly(G.subs(X,H)/vertical**2,X,domain=sp.QQ)
            variants = [F,2*F,F+sp.Poly(X,X)/7]
            try:
                curve(coefficients(sp.Poly(G,X)))
                for f in variants:
                    curve(coefficients(f))
                break
            except ValueError:
                pass
        for suffix, source in zip(('cover','twist','perturbed'),variants):
            rows.append({'id':f'g{genus}-d{degree}-{suffix}','target_genus':genus,'degree':degree,
                         'construction':suffix,'source':coefficients(source),
                         'target':coefficients(sp.Poly(G,X))})
    return {'schema':1,'seed':2026112319,
            'task':'all rational maps u=U(x), v=c*y with deg U=deg F/deg G, plus exact geometric map count',
            'selection_rule':'Seeded smooth targets and inner polynomials. Six genus/degree pairs, each with a constructed cover, quadratic twist, and x/7 coefficient perturbation. Reject singular polynomials only. No solver outcome selected an input.',
            'cases':rows}


if __name__=='__main__':
    output = HERE/'map-panel.json'
    if output.exists():
        raise SystemExit('A frozen panel already exists.')
    text = json.dumps(create(),indent=2)+'\n'
    output.write_text(text)
    print(hashlib.sha256(text.encode()).hexdigest())
