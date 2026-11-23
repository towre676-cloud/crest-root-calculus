"""Deterministic inputs frozen before benchmarking the marked-quartic solvers."""
import hashlib
import json
from pathlib import Path
import random
import sys

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'chart-one-cofactors'))
from chart import X, adapted_frame


def terms(expression):
    return {'terms': [[*power, str(coefficient)]
                      for power, coefficient in sp.Poly(expression, *X).terms()]}


def square_family(parameters):
    a, b, c, d, k = parameters
    t, z = sp.symbols('t z')
    square = t**3 + a*t**2 + b*t*z + c + d*z
    cubic = z**3 + k*z**2 + 2*c*d*z + c*c
    F = sp.S.Zero
    for (r, n), coefficient in sp.Poly(square*square-cubic, t, z).terms():
        j = r+n-2
        i = 4-j-n
        if min(i, j, n) < 0:
            raise ArithmeticError('The selected square shape does not define a quartic.')
        F += coefficient*(-1)**(i+j)*X[0]**i*X[1]**j*X[2]**n
    return terms(F)


def make_panel():
    x, y, z = X
    cases = []
    # Separate parameter choices change the source moduli, not just its coordinates.
    tuples = [(2,-1,1,3,2),(0,2,3,1,-1),(-2,3,2,-1,5),
              (1,0,2,1,3),(3,1,2,2,-2),(2,2,-1,0,7)]
    for i, parameters in enumerate(tuples, 1):
        data = square_family(parameters)
        adapted_frame(data, (0,0,1))
        cases.append({'id': f'square-{i}', 'family': 'constructed square identity',
                      'selection': list(parameters), 'quartic': data, 'mark': [0,0,1]})
    # Uniform coefficient sampling among all homogeneous quartics through [0:0:1].
    # Reject singular sources only; no solver outcome is used for input selection.
    rng = random.Random(2026112302)
    monomials = [(i,j,4-i-j) for i in range(5) for j in range(5-i)
                 if (i,j)!=(0,0)]
    for i, height in enumerate([1,1,2,2,4,4,8,8], 1):
        while True:
            data = {'terms': [[a,b,c,rng.randint(-height,height)] for a,b,c in monomials]}
            try:
                adapted_frame(data, (0,0,1))
                break
            except ValueError:
                pass
        cases.append({'id': f'random-{i}', 'family': 'uniform marked quartic',
                      'selection': {'seed':2026112302, 'height':height},
                      'quartic': data, 'mark': [0,0,1]})
    classical = [
        ('fermat', x**4+y**4-z**4, (1,0,1)),
        ('klein', x**3*y+y**3*z+z**3*x, (1,0,0)),
        ('cyclic', x**4+x*y**3+x*z**3+y*y*z*z, (0,1,0)),
        ('diagonal-2', x**4+2*y**4-z**4, (1,0,1)),
    ]
    for name, F, point in classical:
        data = terms(F)
        adapted_frame(data, point)
        cases.append({'id':name, 'family':'classical equation', 'selection':'fixed equation',
                      'quartic':data, 'mark':list(point)})
    return {'schema':1, 'task':'all geometric transverse degree-three elliptic pencils at a rational mark',
            'selection_rule':'Six square constructions, eight seeded generic marked quartics, four fixed classical curves. Reject only singular sources. No timed result used for selection.',
            'cases':cases}


if __name__ == '__main__':
    output = HERE / 'quartic-panel.json'
    if output.exists():
        raise SystemExit('The frozen panel already exists. Choose a new study to change its inputs.')
    text = json.dumps(make_panel(), indent=2) + '\n'
    output.write_text(text)
    print(hashlib.sha256(text.encode()).hexdigest())
