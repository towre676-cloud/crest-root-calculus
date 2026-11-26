"""Prove that a kernel combination sum c_i theta_{q_i} (same genus, sum c_i = 0, hence a cusp form) has
vanishing coefficients at every square: check sum c_i r_{q_i}(m^2) = 0 for all m <= Sturm bound of the
weight-2 Shimura lift (level of 2q used, which bounds the lift's level from above)."""
import json, sys, time
import cypari2
from inverse import mat

pari = cypari2.Pari()
pari.allocatemem(16 * 10 ** 9)
d, l, gi = map(int, sys.argv[1:4])
K = next(r for r in json.load(open("kernel_probe.json")) if (r["d"], r["l"], r["genus"]) == (d, l, gi))
forms = None
import glob
for f in glob.glob("inert_*.json"):
    for R in json.load(open(f)):
        if (R["d"], R["l"]) == (d, l):
            forms = [x["q"] for x in R["genera"][gi]["forms"]]
B = K["sturm"]
t = time.time()
acc = [[0] * B for _ in K["kernel"]]
for i, q in enumerate(forms):
    if all(v[i] == 0 for v in K["kernel"]):
        continue
    r = pari.qfrep(mat(q), B * B)
    sq = [int(r[m * m - 1]) for m in range(1, B + 1)]
    for k, v in enumerate(K["kernel"]):
        if v[i]:
            acc[k] = [a + v[i] * s for a, s in zip(acc[k], sq)]
    print(f"form {i} done {round(time.time()-t)}s", flush=True)
ok = all(all(x == 0 for x in a) for a in acc)
res = {"d": d, "l": l, "genus": gi, "sturm": B, "kernel_vanishes_at_all_squares_m<=sturm": ok,
       "seconds": round(time.time() - t)}
print("RESULT", res, flush=True)
json.dump(res, open(f"kernel_proof_{d}_{l}_{gi}.json", "w"))
