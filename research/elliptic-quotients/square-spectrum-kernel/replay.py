"""Verify the imported linear kernel and recover it with additional square classes."""
import argparse
from fractions import Fraction
import json
from math import isqrt
from pathlib import Path
from subprocess import run
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "inverse-spectrum"))
from census import residue_field
from verify import sturm_index


def center(field, D, psi):
    p = field.p
    shift = (D-field.trace_omega)//2

    def standard(x, y):
        return field.reduce((x+shift*y, y))

    image_one = standard(*(r[0] for r in psi))
    image_tau = standard(*(2*r[1]-D*r[0] for r in psi))
    plus = field.mul(image_tau, field.inv((-field.trace_omega, 2)))
    a = field.mul((pow(2, -1, p), 0), field.add(image_one, plus))
    b = field.add(image_one, field.neg(a))
    if b == (0, 0) or (field.norm(a)-field.norm(b)+1)%p:
        raise ArithmeticError("Invalid graph representative.")
    return a


def coefficient_rows(binary, D, l, forms, indices):
    field = residue_field(D, l)
    centers = [center(field, D, form["psi"]) for form in forms]
    traces = [(field.trace(a), field.trace(field.mul(a, (0, 1)))) for a in centers]
    content = f"{-D} {field.trace_omega} {l} {len(traces)} {len(indices)}\n"
    content += "".join(f"{s} {t}\n" for s,t in traces)
    content += "\n".join(map(str, indices)) + "\n"
    out = run([str(binary)], input=content, text=True, capture_output=True, check=True)
    values = [list(map(int, row.split())) for row in out.stdout.splitlines()]
    if len(values) != len(indices) or any(len(r) != len(forms) for r in values):
        raise ArithmeticError("Coefficient counter returned malformed data.")
    return dict(zip(indices, values))


def pivot_indices(rows):
    pivots = {}
    for i, v in enumerate(rows):
        w = list(map(Fraction, v))
        for j, value in enumerate(w):
            if not value:
                continue
            if j not in pivots:
                pivots[j] = [x/value for x in w]
                yield i
                if len(pivots) == len(w):
                    return
                break
            for k in range(j, len(w)):
                w[k] -= value*pivots[j][k]


def image(row, relations):
    return [sum(x*y for x,y in zip(row,c)) for c in relations]


def squarefree(n):
    return all(n % (k*k) for k in range(2,isqrt(n)+1))


def audit(binary):
    source = json.loads((HERE.parent/"inverse-spectrum"/"imported-census.json").read_text())
    relations = {(r["d"],r["l"],r["genus"]):r for r in json.loads((HERE/"relations.json").read_text())}
    imported = {}
    for r in json.loads((HERE/"imported-ranks.json").read_text()):
        imported.setdefault((r["d"],r["l"]),[]).append(r)
    results = []
    common_twists = None
    for level in source["levels"]:
        D,l = level["D"],level["ell"]
        rank_records = iter(imported.get((D,l),[]))
        for gi,g in enumerate(level["genera"]):
            forms = g["forms"]
            size = len(forms)
            if size == 1:
                # The pure first coordinate gives a representation of ell^2.
                singleton = coefficient_rows(binary, D, l, forms, [l*l])
                if singleton[l*l][0] < 2:
                    raise ArithmeticError("Missing singleton square representation.")
                results.append({"D":D,"ell":l,"genus":gi,"forms":1,"full_theta_rank":1,
                                "square_rank":1,"square_degree_checked":l})
                continue
            reported = next(rank_records)
            key = (D,l,gi)
            kernel = relations.get(key)
            cutoff = min(reported["sturm"],300)
            if kernel:
                cutoff = sturm_index(g["level"])
                if kernel["sturm"] != cutoff or kernel["forms"] != size:
                    raise ArithmeticError("Wrong Sturm cutoff or number of forms.")
            square_indices = [m*m for m in range(1,cutoff+1)]
            small_indices = list(range(1,3001))
            ts = [t for t in range(2,49) if squarefree(t)]
            extra = {t*m*m for t in ts for m in range(1,31)} if kernel else set()
            matrix = coefficient_rows(binary,D,l,forms,sorted(set(square_indices+small_indices)|extra))
            full_rank = sum(1 for _ in pivot_indices(matrix[n] for n in small_indices))
            sr = sum(1 for _ in pivot_indices(matrix[n] for n in square_indices[:300]))
            if full_rank != size or full_rank != reported["full_rank"] or sr != reported["square_rank"]:
                raise ArithmeticError(f"Imported ranks disagree at {key}.")
            row = {"D":D,"ell":l,"genus":gi,"forms":size,"full_theta_rank":full_rank,
                   "square_rank":sr,"square_degree_checked":cutoff}
            if kernel:
                K = kernel["kernel"]
                if any(len(v)!=size or sum(v) for v in K) or len(K)!=size-sr:
                    raise ArithmeticError("Wrong relation length, constant term, or dimension.")
                if any(any(image(matrix[n],K)) for n in square_indices):
                    raise ArithmeticError(f"Nonzero square-index combination at {key}.")
                if sum(1 for _ in pivot_indices(K))!=len(K):
                    raise ArithmeticError("Kernel vectors are dependent.")
                essential = list(pivot_indices(image(matrix[n],K) for n in small_indices))
                for recovery_indices in ([5,8],[5,20]):
                    if len(list(pivot_indices(image(matrix[n],K) for n in recovery_indices))) != len(K):
                        raise ArithmeticError(f"Additional coefficients do not resolve {key}.")
                resolving = {}
                for t in ts:
                    ms = list(pivot_indices(image(matrix[t*m*m],K) for m in range(1,31)))
                    if len(ms)==len(K):
                        resolving[t] = [m+1 for m in ms]
                common_twists = set(resolving) if common_twists is None else common_twists & set(resolving)
                row.update({"kernel_dimension":len(K),
                            "first_full_theta_indices":[i+1 for i in essential],
                            "first_full_theta_columns":[image(matrix[i+1],K) for i in essential],
                            "resolving_square_classes":resolving})
                print(f"{key}: {len(K)} independent relations verified through degree {cutoff}.",
                      file=sys.stderr,flush=True)
            elif sr!=size:
                raise ArithmeticError("Missing reported kernel.")
            results.append(row)
        if next(rank_records,None) is not None:
            raise ArithmeticError("Extraneous rank record.")
    return {"source_commit":"f340eba16bdd9760739ef6d66a1db74dad045e0b",
            "universal_recovery_indices":[[5,8],[5,20]],
            "common_resolving_square_classes":sorted(common_twists),"genera":results}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compiler", default="g++")
    args = parser.parse_args()
    with TemporaryDirectory() as td:
        binary = Path(td)/"count"
        run([args.compiler,"-std=c++17","-O3",str(HERE/"count.cpp"),"-o",str(binary)],check=True)
        print(json.dumps(audit(binary),indent=2))


if __name__ == "__main__":
    main()
