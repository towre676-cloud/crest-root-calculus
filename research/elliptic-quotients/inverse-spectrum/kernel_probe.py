"""Kernel of the square-coefficient matrix in rank-deficient genera: compute with m <= 300, test persistence
to m <= Mtest, and record the support/coefficients of the kernel vectors (are they ever a difference of two
forms, i.e. a square-equivalent pair? -- no: every pair already has an exact witness)."""
import glob, json, sys
import cypari2
from inverse import mat, level, sturm

pari = cypari2.Pari()
pari.allocatemem(12 * 10 ** 9)


def sq_rows(forms, B):
    rows = []
    for q in forms:
        r = pari.qfrep(mat(q), B * B)
        rows.append([int(r[m * m - 1]) for m in range(1, B + 1)])
    return rows


targets = {(-3, 17), (-7, 17), (-8, 13)}
Mtest = int(sys.argv[1]) if len(sys.argv) > 1 else 800
out = []
for f in sorted(glob.glob("inert_*.json")):
    for R in json.load(open(f)):
        if (R["d"], R["l"]) not in targets:
            continue
        for gi, g in enumerate(R["genera"]):
            forms = [x["q"] for x in g["forms"]]
            rows = sq_rows(forms, Mtest)
            M300 = pari.matrix(300, len(forms), [rows[j][i] for i in range(300) for j in range(len(forms))])
            K = pari.matker(M300)
            kdim = int(pari.matsize(K)[1])
            if kdim == 0:
                print(R["d"], R["l"], gi, "full rank at m<=300", flush=True)
                continue
            Mfull = pari.matrix(Mtest, len(forms), [rows[j][i] for i in range(Mtest) for j in range(len(forms))])
            Kt = pari.matker(Mfull)
            kdim_t = int(pari.matsize(Kt)[1])
            vecs = []
            for j in range(kdim_t):
                v = pari.matrix(len(forms), 1, [Kt[i, j] for i in range(len(forms))])
                v = v * pari.denominator(v)
                c = int(pari.content(v))
                vecs.append([int(v[i, 0]) // c for i in range(len(forms))])
            row = {"d": R["d"], "l": R["l"], "genus": gi, "forms": len(forms), "kernel_dim_m<=300": kdim,
                   f"kernel_dim_m<={Mtest}": kdim_t, "sturm": sturm(g["level_2q"]), "kernel": vecs,
                   "support_sizes": [sum(1 for x in v if x) for v in vecs]}
            out.append(row)
            print({k: v for k, v in row.items() if k != "kernel"}, flush=True)
            for v in vecs:
                print("   kernel vector:", v, flush=True)
json.dump(out, open("kernel_probe.json", "w"), indent=1)
