"""Rank of the square-coefficient matrix of each genus versus the rank of the full theta matrix.

Rows: forms of one genus.  Full matrix: r_q(n), n <= M.  Square matrix: r_q(m^2), m <= sqrt(M').
Both ranks are lower bounds; when a rank equals the number of forms it is exact.  The full-theta rank
equals the number of forms (Schiemann: ternary theta series determine the form; distinct series may
still be dependent, so we compute it).  A square rank below the full rank exhibits a nonzero theta
combination whose coefficients vanish at all squares up to the bound (in all of them once the bound
reaches the Sturm bound of the Shimura lift): the Waldspurger kernel, eigencomponents with L(F,1)=0.
"""
import glob, json, sys
import cypari2
from inverse import mat, level, sturm
from discform import isomorphic as disc_iso

pari = cypari2.Pari()
pari.allocatemem(8 * 10 ** 9)


def rank_rows(rows):
    M = pari.matrix(len(rows), len(rows[0]), [x for r in rows for x in r])
    return int(pari.matrank(M))


def study_genus(forms, Bsq=300, Mfull=3000):
    full = [[int(x) for x in pari.qfrep(mat(q), Mfull)] for q in forms]
    sqr = []
    for q in forms:
        r = pari.qfrep(mat(q), Bsq * Bsq)
        sqr.append([int(r[m * m - 1]) for m in range(1, Bsq + 1)])
    return rank_rows(full), rank_rows(sqr)


if __name__ == "__main__":
    out = []
    for f in sorted(glob.glob("inert_*.json")):
        for R in json.load(open(f)):
            for g in R["genera"]:
                forms = [x["q"] for x in g["forms"]]
                if len(forms) < 2:
                    continue
                N = g["level_2q"]
                B = min(sturm(N), 300)
                rf, rs = study_genus(forms, Bsq=B)
                row = {"d": R["d"], "l": R["l"], "forms": len(forms), "full_rank": rf, "square_rank": rs,
                       "square_bound_used": B, "sturm": sturm(N)}
                out.append(row)
                print(row, flush=True)
    json.dump(out, open("square_rank.json", "w"), indent=1)
