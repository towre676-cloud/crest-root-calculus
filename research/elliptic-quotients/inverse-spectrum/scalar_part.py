"""Scalar part a of a gluing: on F_{l^2} = O_K / l (l inert), psi(t) = a t + b conj(t)."""
import json, sys


def ops(d, l, c=None):
    c = (d * d - d) // 4 if c is None else c
    mul = lambda x, y: ((x[0] * y[0] - x[1] * y[1] * c) % l, (x[0] * y[1] + x[1] * y[0] + x[1] * y[1] * d) % l)
    conj = lambda x: ((x[0] + x[1] * d) % l, (-x[1]) % l)
    sub = lambda x, y: ((x[0] - y[0]) % l, (x[1] - y[1]) % l)
    def inv(x):
        for p in range(l):
            for q in range(l):
                if mul(x, (p, q)) == (1, 0):
                    return (p, q)
    return mul, conj, sub, inv


def parts(d, l, psi, c=None):
    mul, conj, sub, inv = ops(d, l, c)
    w = (0, 1)
    p1 = (psi[0][0] % l, psi[1][0] % l)
    pw = (psi[0][1] % l, psi[1][1] % l)
    s = inv(sub(w, conj(w)))
    a = mul(sub(pw, mul(p1, conj(w))), s)
    b = mul(sub(mul(p1, w), pw), s)
    return a, b


def gaussian(a, l):
    """write p + q w (w = -2 + i) as x + y i."""
    return ((a[0] - 2 * a[1]) % l, a[1] % l)


if __name__ == "__main__":
    d, l = -4, 19
    R = json.load(open("inert_-4_19.json"))[0]
    for gi, g in enumerate(R["genera"]):
        for w in g["witnesses"]:
            if w["degree"] >= 20:
                labs = []
                for qq in w["pair"]:
                    f = next(x for x in g["forms"] if x["q"] == qq)
                    a, b = parts(d, l, f["psi"])
                    labs.append(gaussian(a, l))
                print("genus", gi, "degree", w["degree"], "counts", w["counts"], "a (x + y i):", labs, w["pair"])
