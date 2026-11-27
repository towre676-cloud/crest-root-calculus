"""Locate the two determinant-1440 full-spectrum targets as principal polarisations on C/a1 x C/a2 over Q(sqrt -10)."""
import json, sys
sys.path.insert(0, "../twin-curves"); sys.path.insert(0, "../inverse-spectrum")
from polarised import products, principal_polarisations_fast as pp
from gluing import isometric
c = json.load(open("../opposite-secant-spectra/delay/census.json"))["smallest_square_isospectral"]
targets = {"A": c["child_grams"][0], "B": c["child_grams"][1]}
box = int(sys.argv[1]) if len(sys.argv) > 1 else 25
found = {}
for key, Kf, Bs in products(-40, [1, 3]):
    for p in pp(Kf, Bs, box=box):
        for name, t in targets.items():
            if name not in found and isometric(p["q"], t):
                found[name] = {"product": [key[0], [[str(x) for x in e] for e in key[1]], key[2], [[str(x) for x in e] for e in key[3]]],
                               "H": [str(x) for x in p["H"]], "q": p["q"], "target": t}
                print(name, key[0], key[1], key[2], key[3], p["q"], flush=True)
    if len(found) == 2:
        break
json.dump(found, open("lattices_1440.json", "w"), indent=1)
print("found", sorted(found))
