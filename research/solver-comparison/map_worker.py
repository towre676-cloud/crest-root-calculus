"""Matched constant-vertical maps from structural decomposition or raw elimination."""
import argparse
import json
import os
from pathlib import Path
import shutil
import subprocess
import sys
import tempfile
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent/'bounded-correspondences'))
from composition import recover_constant
from finite import coefficients, rational_points
from maps import X, curve, incidence
from worker import cpu_time


def singular_expression(expression,variables):
    polynomial = sp.Poly(expression,*variables,domain=sp.QQ)
    terms = []
    for powers,coefficient in polynomial.terms():
        monomial = '*'.join(f'{variable}^{exponent}' for variable,exponent in zip(variables,powers)
                            if exponent)
        terms.append(f'({coefficient})'+('*'+monomial if monomial else ''))
    return '+'.join(terms) or '0'


def singular_basis(equations,variables):
    executable = os.environ.get('SINGULAR') or shutil.which('Singular')
    if not executable:
        raise RuntimeError('Set SINGULAR or install Singular.')
    names = ','.join(map(str,variables))
    script = '\n'.join([
        'LIB "primdec.lib";',f'ring r=0,({names}),dp;',
        'short=0; option(redSB); system("random",76123);',
        'ideal I='+','.join(singular_expression(p,variables) for p in equations if p!=0)+';',
        'ideal G=groebner(I);',
        'if (vdim(G)<0) {print("FAIL dimension");quit;}',
        'if (vdim(G)>0) {G=groebner(radical(G));}',
        'int count=vdim(G); print("COUNT "+string(count));',
        f'ring lexring=0,({names}),lp;',
        'ideal L;',
        'if(count==0) {L=1;} else {L=fglm(r,G);}',
        'int k; for(k=1;k<=size(L);k++) {print("POLY "+string(L[k]));}',
        'quit;',
    ])
    with tempfile.TemporaryDirectory(prefix='crest-map-comparison-') as directory:
        path = Path(directory)/'problem.sing'; path.write_text(script)
        result = subprocess.run([executable,'-q',str(path)],text=True,capture_output=True,check=True)
    if '?' in result.stdout or 'FAIL' in result.stdout or 'much slower' in result.stdout:
        raise ArithmeticError(result.stdout[-2000:]+result.stderr[-500:])
    rows = [line[5:] for line in result.stdout.splitlines() if line.startswith('POLY ')]
    counts = [int(line[6:]) for line in result.stdout.splitlines() if line.startswith('COUNT ')]
    if len(counts)!=1 or not rows:
        raise ArithmeticError('Singular did not return a finite basis.')
    context = dict(zip(map(str,variables),variables))
    return [sp.sympify(row.replace('^','**'),locals=context) for row in rows],counts[0]


def corrected_right_decomposition(dense, degree, domain):
    """Recover the monic zero-constant inner polynomial by the exact recurrence.

    If n=m*d and f=K(H), the terms above degree n-1 of f'*H-m*f*H'
    vanish. Their coefficients give the factor i-(m+1)*j used below.
    This independent adapter does not modify the installed SymPy routines.
    """
    n = len(dense) - 1
    if degree < 1 or n % degree:
        raise ValueError('The inner degree must be a positive divisor of the source degree.')
    outer_degree = n // degree
    inner = [domain.one]
    for i in range(1, degree):
        numerator = sum(((i-(outer_degree+1)*j)*dense[i-j]*inner[j]
                         for j in range(i)), domain.zero)
        inner.append(domain.quo(numerator, i*outer_degree*dense[0]))
    return inner + [domain.zero]


def sympy_decomposition(source_values,target_values,corrected=False):
    """Compare stock SymPy or the explicit corrected recurrence, with shared extraction."""
    from sympy.polys.densetools import _dup_left_decompose, _dup_right_decompose
    source,target = curve(source_values),curve(target_values)
    degree = source.degree()//target.degree()
    dense = source.monic().rep.to_list()
    right = (corrected_right_decomposition(dense,degree,sp.QQ) if corrected
             else _dup_right_decompose(dense,degree,sp.QQ))
    left = _dup_left_decompose(dense,right,sp.QQ)
    if left is None:
        return {'geometric_maps':0,'rational_maps':[]}
    inner,outer = sp.Poly.from_list(right,X,domain=sp.QQ),sp.Poly.from_list(left,X,domain=sp.QQ)
    if sp.expand(outer.as_expr().subs(X,inner.as_expr())-source.monic().as_expr())!=0:
        raise ArithmeticError('Conventional decomposition failed its polynomial identity.')
    a,s,inverse = sp.symbols('a s inverse')
    m = target.degree()
    normalized = target.monic()
    shift = (a*outer.nth(m-1)-normalized.nth(m-1))/m
    difference = sp.Poly(normalized.as_expr().subs(X,a*X+shift)-a**m*outer.as_expr(),X)
    constraints = [sp.Poly(q,a,domain=sp.QQ) for q in difference.all_coeffs() if q!=0]
    eliminant = constraints[0]
    for constraint in constraints[1:]:
        eliminant = sp.gcd(eliminant,constraint)
    while eliminant.nth(0)==0:
        eliminant = eliminant.exquo(sp.Poly(a,a))
    equations = [eliminant.monic().as_expr(),source.LC()*s*s-target.LC()*a**m,a*inverse-1]
    variables = [a,s,inverse]
    certificate = rational_points(equations,variables)
    result = []
    for point in certificate['rational_points']:
        assignment = dict(zip(variables,map(sp.Rational,point)))
        horizontal = sp.Poly(assignment[a]*inner.as_expr()+shift.subs(assignment),X)
        vertical = assignment[s]
        if sp.expand(target.as_expr().subs(X,horizontal.as_expr())-source.as_expr()*vertical**2)!=0:
            raise ArithmeticError('A conventional-decomposition map failed its original identity.')
        result.append({'u':coefficients(horizontal),'v_over_y':[str(vertical)],'degree':degree})
    return {'geometric_maps':certificate['geometric_points'],
            'rational_maps':sorted(result,key=lambda value:json.dumps(value,sort_keys=True))}


def solve_case(case,method):
    if method in {'sympy-decomposition','corrected-dense-decomposition'}:
        return sympy_decomposition(case['source'],case['target'],
                                   corrected=method=='corrected-dense-decomposition')
    if method=='crest':
        result = recover_constant(case['source'],case['target'])
        return {'geometric_maps':result['geometric_maps'],
                'rational_maps':sorted(result['rational_maps'],key=lambda v:json.dumps(v,sort_keys=True))}
    source,target = curve(case['source']),curve(case['target'])
    degree = source.degree()//target.degree()
    if source.degree()!=degree*target.degree():
        raise ValueError('Constant-vertical inputs require divisible branch degrees.')
    equations,variables,U,V = incidence(source,target,degree)
    if method=='singular':
        polynomials,count = singular_basis(equations,variables)
    elif method=='sympy-buchberger':
        basis = sp.groebner(equations,*variables,order='grevlex',domain=sp.QQ,method='buchberger')
        polynomials = list(basis)
        count = None
    else:
        raise ValueError('Unknown method.')
    certificate = rational_points(polynomials,variables)
    if count is not None and certificate['geometric_points']!=count:
        raise ArithmeticError('The geometric support count changed in rational extraction.')
    maps = []
    for point in certificate['rational_points']:
        assignment = dict(zip(variables,map(sp.Rational,point)))
        if any(sp.expand(equation).subs(assignment)!=0 for equation in equations):
            raise ArithmeticError('A rational map failed the original coefficient equations.')
        horizontal,vertical = sp.Poly(U.subs(assignment),X),sp.Poly(V.subs(assignment),X)
        maps.append({'u':coefficients(horizontal),'v_over_y':coefficients(vertical),'degree':degree})
    return {'geometric_maps':certificate['geometric_points'],
            'rational_maps':sorted(maps,key=lambda v:json.dumps(v,sort_keys=True))}


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('panel',type=Path)
    parser.add_argument('case')
    parser.add_argument('method',choices=['crest','singular','sympy-buchberger','sympy-decomposition','corrected-dense-decomposition'])
    parser.add_argument('output',type=Path)
    args = parser.parse_args()
    case = next(row for row in json.loads(args.panel.read_text())['cases'] if row['id']==args.case)
    start_wall,start_cpu = time.perf_counter(),cpu_time()
    result = solve_case(case,args.method)
    result['geometric_marks'] = result['geometric_maps']
    result['basis'] = result['rational_maps']
    args.output.write_text(json.dumps({'result':result,
        'solve_wall_seconds':time.perf_counter()-start_wall,
        'solve_cpu_seconds':cpu_time()-start_cpu},indent=2)+'\n')
