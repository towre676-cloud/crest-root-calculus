"""Models over quadratic extensions of the archived quartic invariant field.

Exact inputs: the absolute invariants z1 = J2^5/J10, z2 = J2^3 J4/J10, z3 = J2^2 J6/J10 in F (invariants_F_<d>.json).
Normalise J2 = 1, J10 = 1/z1, J4 = z2/z1, J6 = z3/z1; the weight-zero Mestre coefficients a = A / I2^(deg/2)
come from mestre_formulas.json. The conic is diagonalised over F and split
by adjoining a square root of a ratio of its diagonal entries. The model is
parametrised (mestre.parametrise) and its absolute Igusa invariants are checked EXACTLY against z1, z2, z3."""
import itertools, json, sys
sys.path.insert(0, ".")
import cypari2
import mestre as Me
from igusa_exact import absolute_invariants

pari = Me.pari
pari.allocatemem(128_000_000, 1_000_000_000)
FORM = json.load(open("mestre_formulas.json"))
DEG = {1: 3, 2: 5, 3: 7}


def absinv_exact(S):
    return absolute_invariants(pari, [pari.polcoef(S, k, "x") for k in range(7)])


def build(Fw, z1, z2, z3):
    md = lambda e: pari.Mod(pari(e).subst("x", "w"), Fw)
    z1, z2, z3 = md(z1), md(z2), md(z3)
    Jv = {"J2": pari.Mod(1, Fw), "J4": z2 / z1, "J6": z3 / z1, "J10": 1 / z1}
    ev = lambda k: pari.substvec(pari(FORM[k]), ["J2", "J4", "J6", "J10"], [Jv["J2"], Jv["J4"], Jv["J6"], Jv["J10"]])
    I2 = ev("I2")
    a2, a3 = {}, {}
    for p, q in itertools.combinations_with_replacement((1, 2, 3), 2):
        a2[(p, q)] = ev(f"A{p}{q}") / I2 ** ((DEG[p] + DEG[q]) // 2)
    for p, q, r in itertools.combinations_with_replacement((1, 2, 3), 3):
        e = 1 + DEG[p] + DEG[q] + DEG[r]
        a3[(p, q, r)] = ev(f"A{p}{q}{r}") / I2 ** (e // 2)
    T = [pari("X1"), pari("X2"), pari("X3")]
    L, M = Me.conic_cubic(a2, a3, T)
    # Gram matrix of the conic and diagonalisation over F
    G = [[pari.deriv(pari.deriv(L, vi), vj) / 2 for vj in ("X1", "X2", "X3")] for vi in ("X1", "X2", "X3")]
    G = [[pari.simplify(x) for x in row] for row in G]
    n = 3
    P = [[pari.Mod(1 if i == j else 0, Fw) for j in range(n)] for i in range(n)]      # columns: new basis
    A = [row[:] for row in G]
    for k in range(n):
        if pari.lift(A[k][k]) == 0:
            for j in range(k + 1, n):
                if pari.lift(A[j][j]) != 0 or pari.lift(A[k][j]) != 0:
                    # replace e_k by e_k + e_j
                    for r in range(n):
                        P[r][k] = P[r][k] + P[r][j]
                    A = [[sum(P[a][i] * G[a][b] * P[b][j2] for a in range(n) for b in range(n)) for j2 in range(n)] for i in range(n)]
                    break
        for j in range(k + 1, n):
            c = A[k][j] / A[k][k]
            for r in range(n):
                P[r][j] = P[r][j] - c * P[r][k]
        A = [[sum(P[a][i] * G[a][b] * P[b][j2] for a in range(n) for b in range(n)) for j2 in range(n)] for i in range(n)]
    d = [A[k][k] for k in range(3)]
    assert all(pari.lift(A[i][j]) == 0 for i in range(3) for j in range(3) if i != j)
    # split the conic over K2 = F(sqrt(delta)), delta = -d2/d1: diagonal point (sqrt(delta), 1, 0).
    # (A point over F itself needs the norm equation X^2 + (d2/d1) Y^2 = -d3/d1, whose targets are too tall here.)
    nfF = pari.nfinit(Fw.subst("w", "y"))
    delta = pari.lift(-d[1] / d[0]).subst("w", "y")
    Kabs, a, k = pari.rnfequation(nfF, pari("x^2") - delta, 1)
    Kpol = Kabs.subst("x", "w")
    a = pari.Mod(pari.lift(a).subst("x", "w"), Kpol)                     # image of the F-generator
    root = pari.Mod(pari("w"), Kpol) - int(k) * a                           # sqrt(delta)
    toK = lambda m: pari.Mod(pari.subst(pari.lift(m), "w", pari.lift(a)), Kpol)
    PK = [[toK(P[r][c]) for c in range(3)] for r in range(3)]
    diag_pt = [root, pari.Mod(1, Kpol), pari.Mod(0, Kpol)]
    P0 = [sum(PK[r][c] * diag_pt[c] for c in range(3)) for r in range(3)]
    a2K = {key: toK(v) for key, v in a2.items()}; a3K = {key: toK(v) for key, v in a3.items()}
    LK, MK = Me.conic_cubic(a2K, a3K, T)
    assert pari.lift(pari.simplify(pari.substvec(LK, ["X1", "X2", "X3"], P0))) == 0
    q = Me.parametrise(LK, P0, T, pari("x"))
    S = pari.substvec(MK, ["X1", "X2", "X3"], q)
    ok = all(pari.lift(pari.simplify(x - y)) == 0 for x, y in zip(absinv_exact(S), [toK(z1), toK(z2), toK(z3)]))
    return {"over_F": False, "model": S, "exact_match": ok, "conic_diag": [str(pari.lift(x)) for x in d],
            "model_field": str(Kpol), "delta": str(delta)}


if __name__ == "__main__":
    d = int(sys.argv[1])
    N = json.load(open(f"invariants_F_{d}.json"))
    out = {}
    for name, r in sorted(N.items()):
        if r.get("u"):
            continue
        Fw = pari(r["F"]).subst("x", "w")
        res = build(Fw, r["J2^5/J10"], r["J2^3J4/J10"], r["J2^2J6/J10"])
        rec = {"F": r["F"], "over_F": res["over_F"], "conic_diagonal": res["conic_diag"]}
        if res.get("model") is not None:
            S = res["model"]
            rec["model_field"] = res.get("model_field", r["F"]); rec["delta"] = res.get("delta")
            rec["model_constant_first"] = [str(pari.lift(pari.polcoef(S, k, "x"))).replace("w", "theta") for k in range(7)]
            rec["exact_igusa_match"] = bool(res["exact_match"])
            rec["coefficient_digits"] = max(len(c) for c in rec["model_constant_first"])
        out[name] = rec
        print(name, "| model field degree:", int(pari.poldegree(pari(rec.get("model_field", r["F"])).subst("w", "x"))), "| point over F:", res["over_F"], "| exact Igusa match:", rec.get("exact_igusa_match"),
              "| largest coefficient (chars):", rec.get("coefficient_digits"), flush=True)
        json.dump(out, open(f"mestre_models_{d}.json", "w"), indent=1)
