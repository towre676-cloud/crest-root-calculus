"""Interpolate the Mumford functions of psi and recognise them over L = Q(zeta_8), then verify exactly."""
import json
import sys

import cypari2
import mpmath as mp

from qm import fit
from verify_morphism import verify

pari = cypari2.Pari()
pari.allocatemem(2 * 10 ** 9)
mp.mp.dps = 50
ZETA = mp.expj(mp.pi / 4)
NF = "t^4+1"


def load(which):
    raw = json.load(open(f"samples_{which}.json"))
    return [{k: mp.mpc(mp.mpf(v[0]), mp.mpf(v[1])) for k, v in p.items()} for p in raw]


def to_field(c, bound=10 ** 14):
    """c = a0 + a1 z + a2 z^2 + a3 z^3 with rational a_i (z = zeta_8), via an integer relation."""
    vec = [c, 1, ZETA, ZETA ** 2, ZETA ** 3]
    pv = [pari(f"{mp.nstr(mp.re(v), 48)} + I*({mp.nstr(mp.im(v), 48)})") for v in vec]
    rel = pari.lindep(pv)
    rel = [int(rel[i]) for i in range(5)]
    if rel[0] == 0 or max(abs(r) for r in rel) > bound:
        return None
    a = [pari(f"{-rel[k + 1]}/{rel[0]}") for k in range(4)]
    val = sum(mp.mpf(str(pari.numerator(a[k]))) / mp.mpf(str(pari.denominator(a[k]))) * ZETA ** k for k in range(4))
    if abs(val - c) > mp.mpf(10) ** -35:
        return None
    return "(" + " + ".join(f"({a[k]})*t^{k}" for k in range(4)) + ")"


def poly_str(cs):
    terms = []
    for i, c in enumerate(cs):
        s = to_field(c)
        if s is None:
            return None
        terms.append(f"{s}*x^{i}")
    return " + ".join(terms)


if __name__ == "__main__":
    which = sys.argv[1]
    pts = load(which)
    out = {}
    for key in ("s1", "s0", "g1", "g0"):
        found = None
        for deg in range(1, 16):
            if 2 * deg + 4 > len(pts):
                break
            N, D, sv = fit(pts, key, deg)
            if sv[0] < mp.mpf(10) ** -38 and sv[1] > mp.mpf(10) ** -10:
                Ns, Ds = poly_str(N), poly_str(D)
                found = {"degree": deg, "N": Ns, "D": Ds, "sv": [mp.nstr(sv[0], 3), mp.nstr(sv[1], 3)]}
                break
        out[key] = found
        print(key, found and {k: (v if k != "N" and k != "D" else (v[:120] + "..." if v else None)) for k, v in found.items()}, flush=True)
    json.dump(out, open(f"functions_{which}.json", "w"), indent=1)
    if all(out[k] and out[k]["N"] and out[k]["D"] for k in out):
        rf = lambda k: f"(({out[k]['N']})/({out[k]['D']}))"
        res = verify([-1, 5, -8, 4, -1, 1], rf("s1"), rf("s0"), rf("g1"), rf("g0"), NF)
        print("exact verification:", {k: (str(v) if k == "tangent" else v) for k, v in res.items()})
        json.dump({k: str(v) for k, v in res.items()}, open(f"verified_{which}.json", "w"), indent=1)
