"""Independent check of the residual eigensystems' central behaviour (Cremona coefficients from ecdata)."""
import json
import cypari2
pari = cypari2.Pari()
pari.default("realprecision", 60)
CURVES = {"338c1": ([1, 0, 0, 81, 467], -8), "676b1": ([0, 1, 0, -4, -12], -8), "867e1": ([0, 1, 1, 6, -1], -3)}
out = {}
for lab, (a, D) in CURVES.items():
    E = pari.ellinit(a)
    Et = pari.ellinit(pari.elltwist(E, D))
    rk = pari.ellanalyticrank(Et)
    out[lab] = {"conductor": int(pari.ellglobalred(E)[0]), "rootno": int(pari.ellrootno(E)),
                "twist_by": D, "twist_conductor": int(pari.ellglobalred(Et)[0]),
                "twist_rootno": int(pari.ellrootno(Et)), "twist_analytic_rank": int(rk[0]),
                "twist_leading_coeff": str(rk[1])[:30],
                "untwisted_analytic_rank": int(pari.ellanalyticrank(E)[0])}
    print(lab, out[lab], flush=True)
json.dump(out, open("central_check.json", "w"), indent=1)
