"""Exact expressions of Mestre's coefficients A_ij, A_ijk (transvectant normalisation of mestre.py) and of
I2 = (f,f)_6 as polynomials in the Igusa invariants J2, J4, J6, J10 (even-degree invariants of the binary sextic).
Fitted by exact linear algebra on random integral sextics, then verified on further random sextics."""
import itertools, json, random, sys
sys.path.insert(0, ".")
import mestre as Me
pari = Me.pari
random.seed(7)


def raw(f):
    i = Me.tv(f, 6, f, 6, 4); y1 = Me.tv(f, 6, i, 4, 4); y2 = Me.tv(i, 4, y1, 2, 2); y3 = Me.tv(i, 4, y2, 2, 2)
    ys = {1: y1, 2: y2, 3: y3}
    out = {"I2": Me.tv(f, 6, f, 6, 6)}
    for p, q in itertools.combinations_with_replacement((1, 2, 3), 2):
        out[f"A{p}{q}"] = Me.tv(ys[p], 2, ys[q], 2, 2)
    for p, q, r in itertools.combinations_with_replacement((1, 2, 3), 3):
        out[f"A{p}{q}{r}"] = Me.tv(f, 6, ys[p] * ys[q] * ys[r], 6, 6)
    return {k: pari.simplify(v) for k, v in out.items()}


def degree(key):
    if key == "I2":
        return 2
    D = {1: 3, 2: 5, 3: 7}
    idx = [int(c) for c in key[1:]]
    return sum(D[k] for k in idx) + (1 if len(idx) == 3 else 0)


def monomials(dg):
    return [(a, b, c, e) for a in range(dg // 2 + 1) for b in range(dg // 4 + 1) for c in range(dg // 6 + 1)
            for e in range(dg // 10 + 1) if 2 * a + 4 * b + 6 * c + 10 * e == dg]


def sample():
    co = [random.randint(-9, 9) for _ in range(7)]
    co[6] = random.choice([1, 2, 3])
    P = sum((co[k] * pari("x") ** k for k in range(7)), pari(0))
    if pari.poldisc(P) == 0:
        return sample()
    f = sum((co[k] * pari("x") ** k * pari("y") ** (6 - k) for k in range(7)), pari(0))
    J = pari.genus2igusa(P)
    return raw(f), {"J2": J[0], "J4": J[1], "J6": J[2], "J10": J[4]}


keys = None
samples = [sample() for _ in range(60)]
keys = list(samples[0][0].keys())
formulas = {}
for key in keys:
    mons = monomials(degree(key))
    M = pari.matrix(40, len(mons), [s[1]["J2"] ** a * s[1]["J4"] ** b * s[1]["J6"] ** c * s[1]["J10"] ** e
                                    for s in samples[:40] for (a, b, c, e) in mons])
    rhs = pari.matrix(40, 1, [s[0][key] for s in samples[:40]])
    sol = pari.matinverseimage(M, rhs)
    assert int(pari.matsize(sol)[1]) == 1, key
    expr = sum((sol[k, 0] * pari("J2^%d*J4^%d*J6^%d*J10^%d" % mons[k]) for k in range(len(mons))), pari(0))
    for s in samples[40:]:
        assert pari.substvec(expr, ["J2", "J4", "J6", "J10"], [s[1]["J2"], s[1]["J4"], s[1]["J6"], s[1]["J10"]]) == s[0][key], key
    formulas[key] = str(expr)
    print(key, "degree", degree(key), "monomials", len(mons), "verified on 20 further sextics")
json.dump(formulas, open("mestre_formulas.json", "w"), indent=1)
