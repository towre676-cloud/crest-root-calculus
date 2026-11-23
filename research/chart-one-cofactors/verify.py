"""Replay the exact marked chart proofs, original-cofactor checks and degree spectra."""
import argparse
from concurrent.futures import ProcessPoolExecutor
from copy import deepcopy
from pathlib import Path
import itertools
import json

import sympy as sp
from chart import X,R,c,adapted_frame
from solve import describe,select_norms
from proof import verify
from cofactors import evaluate
from family import verify_family,source
import number_field

HERE = Path(__file__).resolve().parent
CSTAR = HERE.parents[1]/"applications/atlas/original/C-star.json"
P = (0,0,1)


def terms(F):
    return {"terms":[[*i,int(c)] for i,c in sp.Poly(F,*X).terms()]}


def expression(data):
    return sum(a*X[0]**i*X[1]**j*X[2]**k for i,j,k,a in data["terms"])


def cases(panel):
    x,y,z = X
    ref = json.loads(CSTAR.read_text())
    positive = source(6)
    yield "cubic-family",positive,P,1
    yield "old-first-mark",ref,(1,0,1),1
    yield "old-second-mark",ref,(-3,1,0),1
    yield "old-negative",ref,(3,1,0),0
    yield "cyclic-quotient",terms(x**4+x*y**3+x*z**3+y*y*z*z),(0,1,0),1
    yield "diagonal-negative",terms(x**4+y**4-z**4),(1,0,1),0
    if not panel:
        return
    # A coordinate change mixes all three input axes and moves the rational mark.
    A = sp.Matrix([[2,1,1],[1,1,1],[1,2,1]])
    move = dict(zip(X,A*sp.Matrix(X)))
    for k in (1,6,11,16):
        pt = tuple(map(int,A.inv()*sp.Matrix(P)))
        yield "mixed-"+str(k),terms(expression(source(k)).xreplace(move)),pt,1
    for k in (21,26,31,36):
        yield "new-k-"+str(k),source(k),P,1
    for i in range(1,9):
        # The full chart calculation itself establishes these negative controls.
        F = x**4+i*x**3*y+(i+1)*y**3*z+z**3*x+(2*i-3)*x*y*y*z+y**4
        yield "perturbed-"+str(i),terms(F),(0,0,1),0


    from math import gcd
    catalogue = json.loads((CSTAR.parent/"catalogue.json").read_text())
    for row in catalogue:
        points = itertools.chain(((1,0,0),(0,1,0),(0,0,1)),itertools.product(range(-3,4),repeat=3))
        point = next(pt for pt in points if gcd(*pt)==1 and
                     not sum(a*pt[0]**i*pt[1]**j*pt[2]**k for i,j,k,a in row["terms"]))
        yield row["id"],{"terms":row["terms"]},point,0


def run_case(job):
    name,data,point,expected = job
    report = describe(data,point,details=True)
    summary = verify(data,point,report,raw=True)
    if summary["geometric_marks"]!=expected:
        raise ArithmeticError(f"Unexpected number of marked pencils on {name}.")
    if name=="cyclic-quotient":
        # The three six-generator circuit outputs are evaluated over QQ here.
        exact = evaluate(data,point,report,[11,5,7])
        if set(exact["identities"])!={"s_support","c_support","b_coordinate"}:
            raise ArithmeticError("Missing rational first-power cofactor identity.")
    if name in ("old-negative","cubic-family"):
        mod = evaluate(data,point,report,[11,5,7],1009)
        if len(mod["identities"])!=(1 if not expected else 3):
            raise ArithmeticError("Incomplete six-generator cofactor replay.")
    return {"case":name,"geometric_marks":expected,"norm_degrees":report["norm_degrees"],
            "eliminant_degree":report["gcd_degree"],"h_exponent":summary["h_exponent"],
            "target_exponent":summary["target_exponent"]}


def algebraic_case(data):
    report = number_field.run(data,details=True)
    count = number_field.verify(data,report)
    if count!=data["expected_marks_per_embedding"]:
        raise ArithmeticError("The computed number of pencils at the algebraic mark differs.")
    return {"case":data["name"],"field_degree":len(data["field"])-1,"geometric_marks_per_embedding":count,
            "rational_norm_degrees":report["rational_norm_degrees"]}


def main(panel=False,workers=1):
    if workers<1:
        raise ValueError("Choose at least one worker.")
    from types import SimpleNamespace
    roots = [c-i for i in range(4)]
    q = roots[0]*roots[1]*roots[2]*roots[3]
    # Every individual row is a zero divisor, but together they generate the unit ideal.
    rows = [q.exquo(root) for root in roots]
    selected = select_norms(SimpleNamespace(q=q,rows=rows))
    if len(selected)!=1 or selected[0][1]!=[1,1,1,1] or selected[0][2].degree()!=0:
        raise ArithmeticError("The fallback missed the first nonzero combined norm.")
    jobs = list(cases(panel))
    if workers==1:
        rows = list(map(run_case,jobs))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            rows = list(pool.map(run_case,jobs))
    stored = json.loads((HERE/"certificate.json").read_text())
    fresh = json.loads((HERE/"quartic.json").read_text())
    if fresh!=source(6):
        raise ArithmeticError("Published positive equation differs from the family.")
    verify(fresh,P,stored,raw=True)
    for error in ("cofactor","omit-direction","drop-boundary","singular-source"):
        broken = deepcopy(stored)
        try:
            if error=="cofactor":
                entry = broken["norm_cofactors"][0][0]
                entry[2] = str(sp.Rational(entry[2])+1)
            elif error=="omit-direction":
                next(f for f in broken["fibres"] if f["geometric_marks"])["geometric_marks"] = 0
            elif error=="drop-boundary":
                broken["fibres"].pop()
            else:
                adapted_frame(terms(X[0]**4+X[1]*X[2]**3),P)
                raise AssertionError("Singular source was accepted.")
            verify(fresh,P,broken)
        except (ValueError,ArithmeticError):
            pass
        else:
            raise ArithmeticError(f"Alteration {error} was accepted.")
    alg = json.loads((HERE/"quadratic.json").read_text())
    stored_alg = json.loads((HERE/"quadratic-certificate.json").read_text())
    if number_field.verify(alg,stored_alg)!=1:
        raise ArithmeticError("The quadratic point should have one geometric transverse pencil.")
    altered = deepcopy(stored_alg)
    altered["fibres"][0]["field_generator"] = []
    try:
        number_field.verify(alg,altered)
    except (ArithmeticError,ValueError):
        pass
    else:
        raise ArithmeticError("A false mark-field embedding passed.")
    more = json.loads((HERE/"algebraic-marks.json").read_text())
    if not panel:
        more = more[:1]
    if workers==1:
        fields = list(map(algebraic_case,more))
    else:
        with ProcessPoolExecutor(max_workers=workers) as pool:
            fields = list(pool.map(algebraic_case,more))
    return {"cases":rows,"extra_algebraic_mark_cases":fields,"quadratic_certificate_replayed":True,
            "family":verify_family(),"vanished_individual_norm_control":True,"invalid_certificates_or_inputs_rejected":5}


if __name__=="__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--panel",action="store_true")
    parser.add_argument("--workers",type=int,default=1)
    args = parser.parse_args()
    print(json.dumps(main(args.panel,args.workers),indent=2))
