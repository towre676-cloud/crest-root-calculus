"""Bound the all-coprime kernel from above and compare with its forced subspace."""
import argparse
from fractions import Fraction
import json
from math import gcd
from pathlib import Path
import sys

from flint import nmod_mat

from star import ambient_genus, fiber, gram_matrix, pari, residue_field


MODULUS = 65521

LEVELS = [
    (-3,5),(-4,7),(-7,5),(-8,5),(-8,7),(-11,7),(-3,11),(-4,11),(-7,13),
    (-8,13),(-3,17),(-4,19),(-3,23),(-3,41),(-3,59),(-4,23),(-4,31),(-4,47),
    (-7,19),(-7,31),(-7,41),(-8,23),(-8,31),(-8,47),(-19,13),(-43,5),(-67,5),(-163,5),
]


def bound_on_genus(records, D, ell, lower, cutoff):
    width = len(records)
    good = [m*m - 1 for m in range(1,cutoff+1) if gcd(m,2*D*ell)==1]
    forms = [gram_matrix(Q) for Q,_ in records]
    columns = []
    for Q in forms:
        data = Q.qfrep(cutoff**2)
        columns.append([int(data[n]) for n in good])
    entries = [x for column in columns for x in column]
    rank_mod = nmod_mat(width,len(good),entries,MODULUS).rank()
    upper = width-rank_mod
    result = {"classes":width,"cutoff":cutoff,"actual_kernel_lower":lower,"actual_kernel_upper":upper}
    if width <= 150:
        aut = [int(Q.qfauto()[0]) for Q in forms]
        M = pari.matrix(width,len(good),entries).mattranspose()
        K = M.matker()
        k = len(K)
        if k>upper:
            raise ArithmeticError("Rational nullity exceeds its modular upper bound.")
        result["actual_kernel_upper"] = k
        if k:
            KT = K.mattranspose()
            W = pari.matdiagonal(aut)
            unweighted = K * (KT*W*K)**(-1) * KT
            diagonals = [aut[j]*Fraction(str(unweighted[j,j])) for j in range(width)]
        else:
            diagonals = [Fraction(0)]*width
        excess = [v-record[1] for v,record in zip(diagonals,records)]
        if sum(diagonals)!=k or sum(excess)!=k-lower or min(excess)<0:
            raise ArithmeticError("Forced subspace trace or projector comparison failed.")
        result["diag_total_upper"] = str(max(diagonals))
        result["diag_residual_upper"] = str(max(excess))
    if result["actual_kernel_upper"]<lower:
        raise ArithmeticError("The incidence lower bound exceeds the finite upper bound.")
    return result


def run(D, ell, ambient):
    records = {"split":[],"nonsplit":[]}
    field = residue_field(D,ell)
    forced = sum(fiber(Q,g,field,records)["forced_kernel"] for Q,g in ambient)
    cutoff = max(300,3*len(records["split"]))
    split = bound_on_genus(records["split"],D,ell,forced,cutoff)
    nonsplit = bound_on_genus(records["nonsplit"],D,ell,0,cutoff)
    if split["classes"]-nonsplit["classes"] != forced:
        raise ArithmeticError("The two finite incidence spaces have the wrong dimension difference.")
    shared_upper = min(split["actual_kernel_upper"]-forced,nonsplit["actual_kernel_upper"])
    return {"D":D,"ell":ell,"split":split,"nonsplit":nonsplit,"residual_kernel_upper":shared_upper}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("levels",nargs="*",help="D:ell")
    parser.add_argument("--verify",type=Path)
    args = parser.parse_args()
    if args.verify and args.levels:
        parser.error("Use a saved receipt or a list of levels.")
    recorded = json.loads(args.verify.read_text()) if args.verify else None
    levels = ([(r["D"],r["ell"]) for r in recorded] if recorded is not None else
              [tuple(map(int,level.split(":"))) for level in args.levels] if args.levels else LEVELS)
    pari.allocatemem(32_000_000,768_000_000,silent=True)
    ambient = {D:ambient_genus(D) for D in sorted({D for D,_ in levels})}
    result = []
    for D,ell in levels:
        r = run(D,ell,ambient[D])
        result.append(r)
        print(r,file=sys.stderr,flush=True)
    if recorded is not None and recorded!=result:
        raise ArithmeticError("The replay differs from the saved receipt.")
    print(json.dumps(result,indent=2))


if __name__ == "__main__":
    main()
