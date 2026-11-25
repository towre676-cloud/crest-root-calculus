"""Full replay: every curve, exact certificates, and the C_2 spectrum against PR #17 through degree 150."""
import json
import argparse
from collections import Counter
from math import gcd, isqrt

from endo import complete_spectrum
from maps import first_quotient
from run_all import CASES


def pr17_spectrum(X):
    cnt = Counter()
    R = isqrt(X) + 1
    for a in range(0, R + 1):
        for b in range(-R, R + 1):
            if (a == 0 and b != 1) or gcd(a, b) != 1:
                continue
            d = (a * a + b * b) if (a % 2 and b % 2) else 2 * (a * a + b * b)
            if d <= X:
                cnt[d] += 1
    return dict(sorted(cnt.items()))


parser = argparse.ArgumentParser(description=__doc__)
parser.add_argument("--output", default="period-replay.json")
args = parser.parse_args()

out = {"curves": {}}
r = complete_spectrum([-64, 0, 84, 0, -21, 0, 1], 150)
got = dict(sorted(Counter(r["subvarieties"].values()).items()))
out["C_2_spectrum_through_150"] = {"periods": got, "pr17_theorem": pr17_spectrum(150), "equal": got == pr17_spectrum(150)}
print("C_2 through 150 equals PR #17:", got == pr17_spectrum(150), flush=True)
for name, f in CASES.items():
    o = first_quotient(f, 6)
    o["f"] = f
    out["curves"][name] = o
    print(name, "minimal degree", o["minimal_degree"], "spectrum", o["spectrum_through_X"],
          "maps verified", [m["verified"] for m in o["maps"]], flush=True)
expected = [(2, 2), (4, 2), (4, 4), (4, 4)]
out["all_verified"] = out["C_2_spectrum_through_150"]["equal"] and all(
    o.get("minimal_degree") == degree and len(o.get("maps", [])) == number
    and all(m["verified"] and m["map_degree"] == degree for m in o["maps"])
    for o, (degree, number) in zip(out["curves"].values(), expected, strict=True)
)
with open(args.output, "w", encoding="utf-8") as handle:
    json.dump(out, handle, indent=1)
print("all verified:", out["all_verified"])
if not out["all_verified"]:
    raise SystemExit("The period replay did not reproduce its controls.")
