"""Does the Humbert form of a nonscalar inert gluing depend only on the scalar part a, and does its genus
follow eps = ((1 + N(a)) / l)?  Exhaustive over all psi with det = -1, b != 0."""
import itertools, sys, json
import cypari2
from gluing import ns_and_q, reduce_form, isometric
from inverse import mat, is_scalar
from discform import isomorphic as disc_iso
from scalar_part import parts, ops

pari = cypari2.Pari()


def norm(d, l, a):
    mul, conj, _, _ = ops(d, l)
    n = mul(a, conj(a))
    assert n[1] == 0
    return n[0]


def run(d, l):
    cls = []
    for a_, b_, c_, e_ in itertools.product(range(l), repeat=4):
        if (a_ * e_ - b_ * c_ + 1) % l:
            continue
        psi = [[a_, b_], [c_, e_]]
        if is_scalar(d, l, psi):
            continue
        r = ns_and_q(d, l, psi)
        q = reduce_form(r["q"])
        a, b = parts(d, l, psi)
        Na = norm(d, l, a)
        eps = int(pari.kronecker(1 + Na, l))
        for c in cls:
            if isometric(c["q"], q):
                c["a"].add(a); c["Na"].add(Na); c["eps"].add(eps)
                break
        else:
            cls.append({"q": q, "a": {a}, "Na": {Na}, "eps": {eps}})
    # genus by discriminant form
    gen = []
    for c in cls:
        for i, g in enumerate(gen):
            if disc_iso(g, c["q"]):
                c["genus"] = i
                break
        else:
            gen.append(c["q"]); c["genus"] = len(gen) - 1
    a_sets_disjoint = all(not (x["a"] & y["a"]) for x, y in itertools.combinations(cls, 2))
    eps_by_genus = {}
    for c in cls:
        eps_by_genus.setdefault(c["genus"], set()).update(c["eps"])
    print(f"d={d} l={l} classes={len(cls)} a-sets disjoint (form is a function of a): {a_sets_disjoint}; "
          f"|a-set| per class {sorted(len(c['a']) for c in cls)}; eps per genus {eps_by_genus}", flush=True)
    return cls


if __name__ == "__main__":
    for spec in sys.argv[1:]:
        d, l = map(int, spec.split(":"))
        run(d, l)
