"""Certify all-coprime square kernels through the depleted Shimura Sturm bound."""
import argparse
from fractions import Fraction
import json
from math import gcd, lcm
from pathlib import Path
import sys

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "forced-square-kernel"))
from star import ambient_genus, fiber, gram_matrix, pari, residue_field
from observe import form_gram, scalar_orbits


LEVELS = [(-8,13),(-3,17)]


def sturm_bound(level, D, ell):
    index = level
    for p in sorted({2,ell,abs(D)} - {4,8}):
        index = index//p*(p+1)
    return index//6


def depleted_sturm(D, ell):
    level = lcm(2*abs(D)*ell*ell, 4)
    return level, sturm_bound(level,D,ell)


def kernel(forms, D, ell, bound):
    indices = [m*m-1 for m in range(1,bound+1) if gcd(m,2*D*ell)==1]
    cols = []
    for form in forms:
        counts = form.qfrep(bound**2)
        cols.append(pari.Col([int(counts[i]) for i in indices]))
    return pari.Mat(cols).matker()


def coefficient_projector(basis, aut):
    if len(basis):
        transpose = basis.mattranspose()
        weight = pari.matdiagonal(aut)
        return basis * (transpose*weight*basis)**(-1) * transpose * weight
    return pari.matrix(len(aut),len(aut))


def projectors(basis, aut, forced):
    projector = coefficient_projector(basis,aut)
    total = [Fraction(str(projector[i,i])) for i in range(len(aut))]
    residual = [v-s for v,s in zip(total,forced)]
    if sum(total)!=len(basis) or min(residual)<0:
        raise ArithmeticError("The projector comparison failed.")
    return total,residual


def exact_genus(records, D, ell, bound):
    forms = [gram_matrix(Q) for Q,_ in records]
    aut = [int(form.qfauto()[0]) for form in forms]
    basis = kernel(forms,D,ell,bound)
    forced = [s for _,s in records]
    total,residual = projectors(basis,aut,forced)
    maximum = max(total)
    at_maximum = []
    if maximum:
        field = residue_field(D,ell)
        labels = [(a,gram_matrix(form_gram(field,a))) for a in scalar_orbits(field)]
        for i,value in enumerate(total):
            if value!=maximum:
                continue
            a = next((a for a,Q in labels if forms[i].qfisom(Q)),None)
            at_maximum.append({"gram":records[i][0],"automorphisms":aut[i],"scalar_parameter":a,
                               "forced":str(forced[i]),"residual":str(residual[i])})
    return {
        "classes":len(forms),"exact_kernel_dimension":len(basis),"forced_dimension":int(sum(forced)),
        "max_total":str(maximum),"max_residual":str(max(residual)),
        "at_maximum":at_maximum,
    }


def run(D, ell, ambient):
    records = {"split":[],"nonsplit":[]}
    field = residue_field(D,ell)
    for Q,automorphisms in ambient:
        fiber(Q,automorphisms,field,records)
    level,bound = depleted_sturm(D,ell)
    return {"D":D,"ell":ell,"depleted_weight_two_level":level,"sturm_bound":bound,
            **{g:exact_genus(values,D,ell,bound) for g,values in records.items()}}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("levels",nargs="*",help="D:ell")
    parser.add_argument("--verify",type=Path)
    args = parser.parse_args()
    if args.verify and args.levels:
        parser.error("Use a saved receipt or a list of levels.")
    levels = [tuple(map(int,level.split(":"))) for level in args.levels] if args.levels else LEVELS
    expected = json.loads(args.verify.read_text()) if args.verify else None
    if expected is not None:
        levels = [(r["D"],r["ell"]) for r in expected]
    pari.allocatemem(32_000_000,768_000_000,silent=True)
    ambient = {D:ambient_genus(D) for D in sorted({D for D,_ in levels})}
    out = []
    for D,ell in levels:
        r = run(D,ell,ambient[D])
        out.append(r)
        print(D,ell,r["sturm_bound"],[(r[g]["exact_kernel_dimension"],r[g]["max_total"],r[g]["max_residual"])
                                      for g in ("split","nonsplit")],file=sys.stderr,flush=True)
    rendered = json.dumps(out,indent=2)
    if expected is not None and json.loads(rendered)!=expected:
        raise ArithmeticError("The recomputed kernel differs from the saved receipt.")
    print(rendered)


if __name__ == "__main__":
    main()
