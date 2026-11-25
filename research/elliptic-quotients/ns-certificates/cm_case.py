"""A CM test: glue the j = 1728 curve u(u-1)(u+1) to its twist u(u-2)(u+2) along each 2-torsion bijection."""
import itertools
import json
import sys
from fractions import Fraction

from biell import bielliptic
from pipeline import certify_curve
from providers import MapProvider, SumProvider, IsoPair
from setup import Reduction


def glue(a, b, beta):
    a = [Fraction(t) for t in a]
    bb = [Fraction(b[beta[i]]) for i in range(3)]
    R = (bb[1] - bb[0]) / (bb[2] - bb[0])
    den = (a[0] - a[1]) - R * (a[0] - a[2])
    s = (R * (a[0] - a[2]) * a[1] - (a[0] - a[1]) * a[2]) / den
    return [t + s for t in a]


def build_case(beta):
    cs = glue([0, 1, -1], [0, 2, -2], beta)
    if 0 in cs or len(set(cs)) < 3:
        return None
    f, recs = bielliptic(cs)
    g2a, g2b = recs[0]["g2"][0], recs[1]["g2"][0]
    assert recs[0]["g3"][0] == 0 and recs[1]["g3"][0] == 0, "not j = 1728 models"
    ratio = g2b / g2a                     # u^4
    from sympy import sqrt, nsimplify, Rational
    r = nsimplify(sqrt(Rational(ratio.numerator, ratio.denominator)))
    if not r.is_Rational:
        raise ValueError("u^2 irrational for this gluing")
    u2 = Fraction(int(r.p), int(r.q))

    def build(p, rr):
        R = Reduction(f, recs, p, rr, force_i=True)
        f1, f2 = MapProvider(R, 0), MapProvider(R, 1)
        src, dst = (recs[0]["g2"][0], 0), (recs[1]["g2"][0], 0)
        rho = IsoPair(R.K, src, dst, u2)
        rhoi = IsoPair(R.K, src, dst, u2, cm_i=True)
        return R, [f1, f2, SumProvider(f1, f2, rho, None), SumProvider(f1, f2, rhoi, None)]
    return f, recs, build, cs


def cm_case(which, replay=None):
    beta = tuple(map(int, which))
    f, recs, build, cs = build_case(beta)
    bound = replay["X"] if replay is not None else 150
    result = certify_curve(f"CM j=1728 gluing beta={beta}", f, recs, build, True, X=bound, pair=(0, 1),
                           isogenous_by="explicit isomorphism E1 -> E2 over Q(sqrt(u^2))", replay=replay)
    result["c"] = [str(c) for c in cs]
    return result


if __name__ == "__main__":
    which = sys.argv[1]
    result = cm_case(which)
    with open(f"cert_cm_{which}.json", "w", encoding="utf-8") as handle:
        json.dump(result, handle, indent=1, default=str)
    print(json.dumps({k: v for k, v in result.items() if k in
                     ("curve", "status", "rank_gram", "rank_certificate", "NS_discriminant",
                      "spectrum_through_X", "saturation_log", "seconds")}))
