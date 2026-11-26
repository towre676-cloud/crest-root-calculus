"""Locate every twin Humbert form of one discriminant as a principal polarisation on C/a1 x C/a2."""
import json, sys
sys.path.insert(0, "."); sys.path.insert(0, "../inverse-spectrum")
from polarised import products, principal_polarisations_fast as principal_polarisations
from gluing import isometric
d = int(sys.argv[1]); box = int(sys.argv[2]); fs = [int(x) for x in sys.argv[3].split(",")]
pairs = [c for c in json.load(open("../walk-to-products/collision_counts.json")) if c["d"] == d]
targets = {f"p{k}_{s}": c[s] for k, c in enumerate(pairs) for s in ("q1", "q2")}
found = {}
for key, Kf, B in products(d, fs):
    for pp in principal_polarisations(Kf, B, box=box):
        for name, t in targets.items():
            if name not in found and isometric(pp["q"], t):
                found[name] = {"product": [key[0], [[str(x) for x in e] for e in key[1]], key[2], [[str(x) for x in e] for e in key[3]]],
                               "H": [str(x) for x in pp["H"]], "q": pp["q"]}
    if len(found) == len(targets):
        break
print("found", len(found), "of", len(targets), sorted(found))
json.dump(found, open(f"twin_lattices_{d}.json", "w"), indent=1)
