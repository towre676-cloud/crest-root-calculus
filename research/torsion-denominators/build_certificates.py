"""Build the stored certificates for C and C' (the Richelot pair)."""
import json
import sys
import time

from certify import certify_curve

C = [-36, 0, 49, 0, -14, 0, 1]
Cp = [-3003, -14850, 11968, 6600, -3575, -990, 154]
plan = [("C", C, 37, 2 ** 16), ("C'", Cp, 37, 2 ** 16), ("C", C, 61, 2 ** 16), ("C'", Cp, 61, 2 ** 16),
        ("C", C, 23, 2 ** 16)]
if "--with-23" in sys.argv:
    plan.append(("C'", Cp, 23, 2 ** 19))
out = []
for name, f, p, lim in plan:
    t = time.time()
    res, _ = certify_curve(f, p, rmax=6, sylow_limit=lim)
    res.update({"name": name, "curve": f, "seconds": round(time.time() - t, 1)})
    out.append(res)
    print(name, p, [B["prime_to_p_degree"] for B in res["blocks"]], res["seconds"], flush=True)
json.dump(out, open("certificates.json", "w"))
