"""End-to-end exact certificate: NS lattice, rank, saturation, complete spectrum."""
from __future__ import annotations

from math import isqrt
import time
from fractions import Fraction

import cypari2

from certify import gram_from_providers, saturate, theta_coords, provider_degree
from isotropic import spectrum_from_lattice
from ns import RATIONAL_CM_J
from setup import Reduction, good_primes

pari = cypari2.Pari()


def jinv(g2, g3):
    """j of Y^2 = 4X^3 - g2 X - g3 for Gaussian-rational (g2, g3) given as (re, im) pairs."""
    w = pari("Mod(w, w^2+1)")
    G2 = pari(str(g2[0])) + pari(str(g2[1])) * w
    G3 = pari(str(g3[0])) + pari(str(g3[1])) * w
    j = 1728 * G2 ** 3 / (G2 ** 3 - 27 * G3 ** 2)
    jl = pari.lift(j)
    if pari.poldegree(jl, "w") <= 0:
        return Fraction(str(pari.polcoef(jl, 0, "w")))
    return str(j)


def frob_data(R, pv):
    E = pv.E
    p = R.K.p
    assert R.K.r == 1
    a = p + 1 - (sum(1 for _ in E.points()) + 1)
    disc = a * a - 4 * p
    field = int(pari.coredisc(disc))
    return a, field, a % p != 0


def ell_trace(g2g3, p, i_val=None):
    """Trace of Frobenius of Y^2 = 4X^3 - g2 X - g3 over F_p (Gaussian coefficients use i -> i_val)."""
    def red(v):
        re_, im_ = Fraction(v[0]), Fraction(v[1])
        out = re_.numerator * pow(re_.denominator, -1, p)
        if im_:
            out += i_val * im_.numerator * pow(im_.denominator, -1, p)
        return out % p
    G2, G3 = red(g2g3[0]), red(g2g3[1])
    sq = [0] * p
    for t in range(p):
        sq[t * t % p] += 1
    n = 1 + sum(sq[(4 * X ** 3 - G2 * X - G3) % p] for X in range(p))
    return p + 1 - n


def rank_certificate(f, recs, build, need_i, pa, pb, isogenous_by, Gram):
    """Upper bound rho <= 2 + rank Hom(E_pa, E_pb), with exact witnesses."""
    out = {"pair": [pa, pb]}
    assert Gram[pa + 1][pb + 1] > 0, "chosen pair does not span J"
    ga = (recs[pa]["g2"], recs[pa]["g3"])
    gb = (recs[pb]["g2"], recs[pb]["g3"])
    ja = jinv(*ga)
    primes = []
    for p in good_primes(f, recs, need_i, start=11, stop=4000):
        primes.append(p)
        if len(primes) >= 60:
            break

    def ival(p):
        if p % 4 != 1:
            return None
        s = next(t for t in range(2, p) if pow(t, (p - 1) // 2, p) == p - 1)
        return pow(s, (p - 1) // 4, p)
    if isogenous_by:
        out["isogeny_witness"] = isogenous_by
        if isinstance(ja, Fraction) and ja.denominator == 1 and int(ja) in RATIONAL_CM_J:
            out["cm_witness"] = {"j": str(ja), "list": "the 13 rational CM j-invariants (Heegner-Stark)"}
            out["rank_Hom"] = 2
        else:
            fields = {}
            for p in primes:
                a = ell_trace(ga, p, ival(p))
                if a % p:
                    fields.setdefault(int(pari.coredisc(a * a - 4 * p)), p)
                if len(fields) >= 2:
                    break
            if len(fields) < 2:
                raise RuntimeError("could not certify non-CM")
            out["non_cm_witness"] = {str(k): v for k, v in fields.items()}
            out["rank_Hom"] = 1
    else:
        for p in primes:
            a1, a2 = ell_trace(ga, p, ival(p)), ell_trace(gb, p, ival(p))
            t0, t1, u0, u1 = 2, a1, 2, a2
            for _ in range(119):
                t0, t1 = t1, a1 * t1 - p * t0
                u0, u1 = u1, a2 * u1 - p * u0
            if t1 != u1:
                out["non_isogeny_witness"] = {"p": p, "traces": [a1, a2]}
                out["rank_Hom"] = 0
                break
        if "rank_Hom" not in out:
            raise RuntimeError("no non-isogeny witness found")
    out["rho_upper"] = 2 + out["rank_Hom"]
    return out


def _gram_at(build, p):
    R, pvs = build(p, 1)
    degs = []
    for pv in pvs:
        d = getattr(pv, "deg", None)
        if d is None:
            d1, d2, k = pv.f1.deg, pv.f2.deg, pv.rho.k
            bound = d1+k*d2+isqrt(4*d1*k*d2)+1
            d = provider_degree(pv, bound)
            pv.deg = d
        degs.append(d)
    return degs, gram_from_providers(pvs, degs)


def certify_curve(name, f, recs, build, need_i, X=100, pair=(0, 1), isogenous_by=None, gram_prime_start=1009, replay=None):
    t = time.time()
    ps = []
    if replay is not None:
        stored = replay["gram_prime"]
        ps = [stored] if isinstance(stored, int) else list(stored)
    if len(ps) < 2:
        for p in good_primes(f, recs, need_i, start=gram_prime_start, stop=20000):
            if p not in ps and _buildable(build, p):
                ps.append(p)
            if len(ps) == 2:
                break
    if len(ps) != 2 or ps[0] == ps[1]:
        raise ValueError("Two distinct good primes are required for the Gram replay.")
    degs, G = _gram_at(build, ps[0])
    degs2, G2 = _gram_at(build, ps[1])
    assert degs == degs2 and G == G2, "Gram matrix differs between primes"
    p0 = ps
    rank_gram = int(pari.matrank(pari.matrix(len(G), len(G), [x for r in G for x in r])))
    rc = rank_certificate(f, recs, build, need_i, pair[0], pair[1], isogenous_by, G)
    out = {"curve": name, "f": f, "gram_prime": p0, "provider_degrees": degs, "gram": G, "rank_gram": rank_gram,
           "rank_certificate": rc}
    if rank_gram != rc["rho_upper"]:
        out["status"] = "incomplete: need more independent providers"
        return out
    log = []
    B, GB, det = saturate(build, len(degs), G, log, need_i, f, recs,
                         replay=replay["saturation_log"] if replay is not None else None)
    tv = theta_coords(B)
    spec, classes = spectrum_from_lattice(GB, tv, X)
    out.update({"saturation_log": log, "NS_basis_in_generators": [[str(x) for x in v] for v in B],
                "NS_gram": GB, "NS_discriminant": det, "theta": tv,
                "spectrum_through_X": {str(k): v for k, v in spec.items()}, "X": X,
                "status": "certified", "seconds": round(time.time() - t, 1)})
    return out


def _buildable(build, p):
    try:
        build(p, 1)
        return True
    except (AssertionError, ValueError, NotImplementedError):
        return False
