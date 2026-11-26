"""Independent replay of two claims: (1) D = -163, l = 5: scalar parts 0 and sqrt(-163) give one genus,
agreeing through degree 41 and 8 vs 0 at degree 42 (signed); (2) Gaussian order Z[f i], l = 7, f = 1 mod 7:
the gluing forms include q0 = (7z)^2 + 4x^2 + 4 f^2 y^2 and q1 = (7z + 2y)^2 + 4x^2 + 4 f^2 y^2, first
separated at degree 2f."""
import itertools, json, sys
import cypari2
from gluing import ns_and_q, reduce_form, isometric
from inverse import mat, is_scalar, level, sturm
from discform import isomorphic as disc_iso
from scalar_part import parts, ops
from sweep_inert import prim_counts

pari = cypari2.Pari()


def classes(d, l, c=None):
    cls = []
    for a_, b_, c_, e_ in itertools.product(range(l), repeat=4):
        if (a_ * e_ - b_ * c_ + 1) % l:
            continue
        psi = [[a_, b_], [c_, e_]]
        if is_scalar(d, l, psi, c):
            continue
        r = ns_and_q(d, l, psi, c)
        if r is None:
            continue
        q = reduce_form(r["q"])
        a, _ = parts(d, l, psi, c)
        for cl in cls:
            if isometric(cl["q"], q):
                cl["a"].add(a); break
        else:
            cls.append({"q": q, "a": {a}})
    return cls


def first_diff(q1, q2, B):
    p1, p2 = prim_counts(q1, B), prim_counts(q2, B)
    m = next((k + 1 for k in range(B) if p1[k] != p2[k]), None)
    return m, (p1[m - 1], p2[m - 1]) if m else None


out = {}
# (1) D = -163, l = 5
cls = classes(-163, 5)
w_sqrt = None
mul, conj, sub, inv = ops(-163, 5)
sq = sub((0, 1), conj((0, 1)))                 # w - conj(w) = sqrt(-163)
c0 = next(c for c in cls if (0, 0) in c["a"])
c1 = next(c for c in cls if sq in c["a"])
m, cnt = first_diff(c0["q"], c1["q"], 80)
out["D-163_l5"] = {"classes": len(cls), "q_a=0": c0["q"], "q_a=sqrt-163": c1["q"], "same_genus": disc_iso(c0["q"], c1["q"]),
                   "first_degree": m, "counts_mod_sign": cnt}
print("D=-163 l=5:", out["D-163_l5"], flush=True)
# (2) conductor family
for f in (8, 15):
    cls = classes(0, 7, f * f)
    q0 = [[4, 0, 0], [0, 4 * f * f, 0], [0, 0, 49]]
    q1 = [[4, 0, 0], [0, 4 * f * f + 4, 14], [0, 14, 49]]
    h0 = [c for c in cls if isometric(c["q"], q0)]
    h1 = [c for c in cls if isometric(c["q"], q1)]
    m, cnt = first_diff(q0, q1, 4 * f)
    out[f"Z[{f}i]_l7"] = {"classes": len(cls), "q0_occurs": bool(h0), "q1_occurs": bool(h1),
                          "a_of_q0": sorted(h0[0]["a"]) if h0 else None, "a_of_q1": sorted(h1[0]["a"]) if h1 else None,
                          "same_genus": disc_iso(q0, q1), "first_degree": m, "counts_mod_sign": cnt}
    print(f"f={f}:", out[f"Z[{f}i]_l7"], flush=True)
json.dump(out, open("user_claims_replay.json", "w"), indent=1, default=str)
