"""Solve the same marked-quartic support with the project or a general solver."""
import argparse
import json
import os
from pathlib import Path
import resource
import shutil
import subprocess
import sys
import tempfile
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / 'chart-one-cofactors'))
from chart import C, S, R, Chart, adapted_frame, build_from_frame, hasse_table
from proof import crt_targets, original_rows
from solve import describe, fibre


def raw_system(data, point):
    """Equivalent generators on h != 0, before the project's norm elimination."""
    source, frame = adapted_frame(data, point)
    table = hasse_table(source, R)
    h, u, v, w = [table[0,j] for j in range(4)]
    q = (4*h*v-u*u).monic()
    denominator = 12*h*h*w
    numerator = -(4*h*h*table[1,2]-2*h*u*table[1,1]+u*u*table[1,0])
    chart = Chart(frame, table, q, h, denominator, numerator, [])
    algebra, b, cast, original = original_rows(chart, source)
    generators = [cast(q), cast(denominator)*b-cast(numerator)]
    for position in ((2,0),(2,1),(1,0),(1,1)):
        polynomial = original[position]
        while polynomial:
            quotients, remainder = polynomial.div([cast(h)])
            if remainder:
                break
            polynomial = quotients[0]
        generators.append(polynomial)
    return source, frame, generators, cast(h)


def singular_polynomial(polynomial):
    terms = []
    for powers, coefficient in polynomial.items():
        monomial = '*'.join(f'{variable}^{exponent}'
                            for variable, exponent in zip(polynomial.ring.gens, powers)
                            if exponent)
        terms.append(f'({coefficient})' + ('*'+monomial if monomial else ''))
    return '+'.join(terms) or '0'


def singular_support(generators, h):
    executable = os.environ.get('SINGULAR') or shutil.which('Singular')
    if not executable:
        raise RuntimeError('Set SINGULAR to a Singular executable or install Singular.')
    script = '\n'.join([
        'LIB "primdec.lib";',
        'ring r=0,(qinv,b,c,s),dp;',
        'short=0; option(redSB); system("random",76123);',
        'ideal I='+','.join(singular_polynomial(p) for p in generators)
        +',1-qinv*('+singular_polynomial(h)+');',
        'ideal G=groebner(I);',
        'if (vdim(G)<0) {print("FAIL dimension");quit;}',
        'if (vdim(G)>0) {G=groebner(radical(G));}',
        'int count=vdim(G); print("COUNT "+string(count));',
        'ring lexring=0,(qinv,b,c,s),lp;',
        'ideal L;',
        'if (count==0) {L=1;} else {L=fglm(r,G);}',
        'poly p; int k;',
        'for(k=1;k<=size(L);k++) {',
        ' p=L[k];',
        ' if (subst(subst(p,qinv,0),b,0)==p) {print("POLY "+string(p));}',
        '}',
        'quit;',
    ])
    with tempfile.TemporaryDirectory(prefix='crest-comparison-') as directory:
        path = Path(directory)/'problem.sing'
        path.write_text(script)
        result = subprocess.run([executable, '-q', str(path)], text=True,
                                capture_output=True, check=True)
    if '?' in result.stdout or 'FAIL' in result.stdout or 'much slower' in result.stdout:
        raise ArithmeticError(result.stdout[-2000:] + result.stderr[-500:])
    rows = [line[5:] for line in result.stdout.splitlines() if line.startswith('POLY ')]
    counts = [int(line[6:]) for line in result.stdout.splitlines() if line.startswith('COUNT ')]
    if len(counts)!=1 or not rows:
        raise ArithmeticError('Singular did not emit a complete support.')
    polynomials = [sp.Poly(row.replace('^','**'), C, S, domain=sp.QQ).as_expr() for row in rows]
    return polynomials, counts[0]


def sympy_support(generators, h):
    inverse, b, c, s = sp.symbols('qinv b c s')
    equations = [p.as_expr() for p in generators]+[1-inverse*h.as_expr()]
    basis = sp.groebner(equations, inverse, b, c, s, order='grevlex',
                       domain=sp.QQ, method='buchberger')
    if list(basis)==[1]:
        return [sp.S.One], 0
    if not basis.is_zero_dimensional:
        raise ArithmeticError('The saturated marked chart unexpectedly has positive dimension.')
    lex = basis.fglm('lex')
    return [p.as_expr() for p in lex.polys if not p.as_expr().has(inverse,b)], None


def canonical_support(polynomials):
    """Return a unique rational lex basis for the reduced geometric support."""
    basis = sp.groebner(polynomials, C, S, order='lex', domain=sp.QQ)
    return [[[*powers,str(coefficient)] for powers,coefficient in p.terms()]
            for p in basis.polys]


def recover_support(source, frame, polynomials):
    """Recover every target and remove nilpotents using the original fibre equations."""
    basis = sp.groebner(polynomials, C, S, order='lex', domain=sp.QQ)
    if list(basis)==[1]:
        return {'basis':canonical_support([1]), 'geometric_marks':0, 'fibres':[]}
    univariate = next((p.as_expr() for p in basis.polys if not p.as_expr().has(C)),None)
    if univariate is None:
        raise ArithmeticError('A finite projected support needs a univariate eliminant.')
    chart = build_from_frame(source, frame, R)
    records = [fibre(chart,p.monic(),1)[0]
               for p,m in sp.Poly(univariate,S,domain=sp.QQ).factor_list()[1]]
    support = [p.as_expr() for label,p in crt_targets(records)]
    return {'basis':canonical_support(support),
            'geometric_marks':sum(record['geometric_marks'] for record in records),
            'fibres':[record for record in records if record['geometric_marks']]}


def solve_case(case, method):
    data, point = case['quartic'], case['mark']
    if method=='crest':
        report = describe(data,point,details=True)
        support = [polynomial.as_expr() for label,polynomial in crt_targets(report['fibres'])]
        return {'basis':canonical_support(support), 'geometric_marks':report['geometric_marks'],
                'fibres':[record for record in report['fibres'] if record['geometric_marks']],
                'certificate':report}
    source,frame,generators,h = raw_system(data,point)
    if method=='singular':
        support, count = singular_support(generators,h)
    elif method=='sympy-buchberger':
        support, count = sympy_support(generators,h)
    else:
        raise ValueError('Unknown method.')
    result = recover_support(source,frame,support)
    if count is not None and result['geometric_marks']!=count:
        raise ArithmeticError('The reconstructed map count disagrees with the radical algebra.')
    return result


def cpu_time():
    own = resource.getrusage(resource.RUSAGE_SELF)
    child = resource.getrusage(resource.RUSAGE_CHILDREN)
    return own.ru_utime+own.ru_stime+child.ru_utime+child.ru_stime


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('panel', type=Path)
    parser.add_argument('case')
    parser.add_argument('method', choices=['crest','singular','sympy-buchberger'])
    parser.add_argument('output', type=Path)
    args = parser.parse_args()
    case = next(case for case in json.loads(args.panel.read_text())['cases'] if case['id']==args.case)
    start_wall, start_cpu = time.perf_counter(), cpu_time()
    result = solve_case(case,args.method)
    measurements = {'solve_wall_seconds':time.perf_counter()-start_wall,
                    'solve_cpu_seconds':cpu_time()-start_cpu}
    args.output.write_text(json.dumps({'result':result,**measurements},indent=2)+'\n')
