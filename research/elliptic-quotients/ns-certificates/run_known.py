"""Certify the known curves: C, C', T1, T2 (period maps) and the six k = 3 bielliptic gluings."""
import json
import sys
from pathlib import Path

from biell import bielliptic
from load import load_curve
from pipeline import certify_curve
from providers import MapProvider, SumProvider, IsogenyPair
from setup import Reduction

results = []


def period_case(prefix, pair, isog=None, use=None, replay=None):
    name, f, recs = load_curve(prefix)
    if use is not None:
        recs = [recs[i] for i in use]
    need_i = any(r["gaussian"] for r in recs)

    def build(p, r):
        R = Reduction(f, recs, p, r)
        return R, [MapProvider(R, i) for i in range(len(recs))]
    return certify_curve(name, f, recs, build, need_i, X=150, pair=pair, isogenous_by=isog, replay=replay)


def biell_case(idx, replay=None):
    path = Path(__file__).resolve().parent.parent / "bielliptic-reconstruction/test_curves.json"
    cases = json.loads(path.read_text(encoding="utf-8"))
    case = cases[idx]
    f, recs = bielliptic(case["unscrambled_c"])

    def build(p, r):
        R = Reduction(f, recs, p, r)
        f1, f2 = MapProvider(R, 0), MapProvider(R, 1)
        rho = IsogenyPair(R.K, (recs[0]["g2"][0], recs[0]["g3"][0]), (recs[1]["g2"][0], recs[1]["g3"][0]), 3)
        rho.calibrate([P for P in f1.aux(30) if P is not None][:8])
        return R, [f1, f2, SumProvider(f1, f2, rho, None)]
    o = certify_curve(f"k=3 bielliptic, gluing {case['type']} beta={case['beta']}", f, recs, build, False, X=150,
                      pair=(0, 1), isogenous_by="explicit global 3-isogeny (PARI ellisomat, Velu), reduced exactly", replay=replay)
    o["expected_type"] = case["type"]
    return o


def known_case(which, replay=None):
    if which == "C":
        return period_case("C (", (0, 1), replay=replay)
    if which == "Cp":
        return period_case("C'", (0, 1), replay=replay)
    if which in ("T1", "T2"):
        return period_case(which, (2, 3), isog="deg mu(f_0, f_2) > 0 in the Gram matrix", replay=replay)
    return biell_case(int(which), replay=replay)


if __name__ == "__main__":
    which = sys.argv[1]
    o = known_case(which)
    json.dump(o, open(f"cert_{which}.json", "w"), indent=1, default=str)
    print(json.dumps({k: v for k, v in o.items() if k in
                     ("curve", "status", "rank_gram", "rank_certificate", "NS_discriminant",
                      "spectrum_through_X", "saturation_log", "seconds")}, default=str))
