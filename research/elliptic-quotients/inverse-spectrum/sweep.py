import json, sys, time
from analyse import analyse, cross_witnesses
from inverse import run
out = {}
for d in [int(x) for x in sys.argv[1].split(",")]:
    ls = [int(x) for x in sys.argv[2].split(",")]
    t = time.time()
    cl = run([d], ls)[d]
    rep = analyse(cl)
    cw = cross_witnesses(rep, int(sys.argv[3]) if len(sys.argv) > 3 else 200)
    out[d] = {"genera": rep, "cross": cw}
    for r in rep:
        sizes = [len(k["members"]) for k in r["square_classes"]]
        print(f"SUM d={d} l={r['levels_l']} det={r['det']} N={r['level_2q']} B={r['sturm_bound']} "
              f"iso={r['isometry_classes']} sq={len(sizes)} sizes={sizes}", flush=True)
    print(f"SUM d={d} cross-genus pairs separated {cw['pairs_separated']} unresolved {len(cw['unresolved_to_B'])}"
          f" {round(time.time()-t,1)}s", flush=True)
    json.dump({str(k): v for k, v in out.items()}, open(f"sweep_{sys.argv[1].replace(',','_')}.json", "w"), indent=1)
