"""Replay of the parallel branch's kernel formula dim V_split = h_split - h_nonsplit + dim V_nonsplit
(equal coprime observation rank in the two genera) on the 12 levels of uncertainty.json."""
import json
rows = []
for r in json.load(open("uncertainty.json")):
    (h1, k1), (h2, k2) = sorted([(g["size"], g["dim_V"]) for g in r["genera"]], reverse=True)
    rows.append({"d": r["d"], "l": r["l"], "h_split": h1, "dimV_split": k1, "h_nonsplit": h2, "dimV_nonsplit": k2,
                 "forced": h1 - h2, "residual": k2, "formula_holds": h1 - h2 + k2 == k1, "equal_rank": h1 - k1 == h2 - k2})
    print(rows[-1])
json.dump(rows, open("forced_check.json", "w"), indent=1)
