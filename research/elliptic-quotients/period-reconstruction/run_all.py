"""Equation-only first-quotient reconstruction on the validation and target curves."""
import json
import sys
import time

from maps import first_quotient

CASES = {
    "C (torsion study, multiplicity one)": [-36, 0, 49, 0, -14, 0, 1],
    "C' (torsion study, multiplicity one, no involution)": [-3003, -14850, 11968, 6600, -3575, -990, 154],
    "T1: Richelot neighbour of C_2, pairing (1,-2)(-1,-4)(2,4)": [-1792, 7200, 5268, -900, -1605, -90, 154],
    "T2: Richelot neighbour of C_2, pairing (1,2)(-1,-4)(-2,4)": [-9856, 1440, 6420, 900, -1317, -450, 28],
}
if __name__ != "__main__":
    names = []
else:
    names = sys.argv[1:] or list(CASES)
out = {}
for n in names:
    t = time.time()
    o = first_quotient(CASES[n], 6)
    o["f"] = CASES[n]
    o["seconds"] = round(time.time() - t, 1)
    out[n] = o
    print(n, o.get("minimal_degree"), [(m["verified"], m["map_degree"], m["differential"]) for m in o.get("maps", [])], o["seconds"], flush=True)
if names:
    json.dump(out, open(f"first_quotients_{abs(hash(tuple(names))) % 10000}.json", "w"), indent=1)
