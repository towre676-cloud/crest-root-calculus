"""Numerical Richelot graph from the certified bielliptic models (discovery only; exact witnesses are separate).

Starting from each bielliptic quartic-field model at every embedding of F, all 15 Richelot neighbours are formed
from the numerical branch points; each neighbour's absolute Igusa triple is compared with the recognised triples of
all forms at all embeddings.  Breadth-first to a given depth; with --twin-only the search only continues from twin forms.  Writes richelot_graph_<D>.json.

Usage: python richelot_explore.py D depth [--twin-only]
"""
import itertools, json, sys
import cypari2

pari = cypari2.Pari(); pari.allocatemem(2 * 10 ** 9); pari.set_real_precision(150)
KEYS = ("J2^5/J10", "J2^3J4/J10", "J2^2J6/J10")
x = pari("x")


def absinv(P):
    J2, J4, J6, J8, J10 = pari.genus2igusa(P)
    return [J2 ** 5 / J10, J2 ** 3 * J4 / J10, J2 ** 2 * J6 / J10]


def pairings():
    out = []
    def rec(rest, acc):
        if not rest:
            out.append(acc); return
        a = rest[0]
        for b in rest[1:]:
            rec([c for c in rest if c not in (a, b)], acc + [(a, b)])
    rec(list(range(6)), [])
    return out


def neighbours(P):
    r = pari.polroots(P, precision=600)
    lead = pari.polcoef(P, 6)
    res = []
    for pr in pairings():
        G = [(x - r[a]) * (x - r[b]) for a, b in pr]
        M = pari.matrix(3, 3, [pari.polcoef(g, i) for g in G for i in range(3)])
        if abs(pari.matdet(M)) < pari("1e-100"):
            continue
        br = [pari.deriv(G[i]) * G[(i + 1) % 3] - G[i] * pari.deriv(G[(i + 1) % 3]) for i in range(3)]
        Q = br[0] * br[1] * br[2]
        res.append((pr, Q))
    return res


def close(a, b, tol=pari("1e-60")):
    return all(abs(p - q) <= tol * (1 + abs(q)) for p, q in zip(a, b))


if __name__ == "__main__":
    D = int(sys.argv[1]); depth = int(sys.argv[2]); TWIN_ONLY = "--twin-only" in sys.argv
    tag = "_twin" if TWIN_ONLY else ""
    inv = json.load(open(f"invariants_F_{D}.json"))
    Fs = next(iter(inv.values()))["F"]
    thetas = pari.polroots(pari(Fs), precision=600)
    ev = lambda e, t: pari.subst(pari(e), "x", t)
    targets = {(n, k): [ev(r[key], t) for key in KEYS] for n, r in inv.items() for k, t in enumerate(thetas)}
    starts = {}
    for n, r in inv.items():
        if r.get("u") and r.get("v"):
            for k, t in enumerate(thetas):
                u, v = ev(r["u"], t), ev(r["v"], t)
                s1 = ((u + pari.sqrt(u * u - 4 * v ** 3)) / 2) ** (pari(1) / 3)
                starts[(n, k)] = x ** 6 - s1 * x ** 4 + (v / s1) * x ** 2 - 1
    edges, seen, nodes, hashes = [], {}, [], set()
    hkey = lambda I: tuple(str(pari.precision(v, 30))[:25] for v in (pari.real(I[0]), pari.imag(I[0])))
    frontier = []
    for key, P in starts.items():
        assert close(absinv(P), targets[key])
        seen[key] = 0; nodes.append((absinv(P), key)); hashes.add(hkey(absinv(P))); frontier.append((key, P, [list(key)]))
    for step in range(1, depth + 1):
        new = []
        for key, P, path in frontier:
            for pr, Q in neighbours(P):
                try:
                    I = absinv(Q)
                except Exception:
                    continue
                h = hkey(I)
                if h in hashes:
                    continue
                hashes.add(h)
                hit = [t for t, T in targets.items() if close(I, T)]
                label = hit[0] if hit else ("anon", len(nodes))
                nodes.append((I, label))
                p2 = path + [[list(p) for p in pr], list(label)]
                if hit:
                    edges.append({"path": p2, "step": step})
                    seen.setdefault(hit[0], step)
                if hit or not TWIN_ONLY:
                    new.append((label, Q, p2))
        frontier = new
        reached = sorted({n for (n, k) in seen})
        print("depth", step, "new nodes", len(new), "twin forms reached", len(reached), reached, flush=True)
    json.dump([[str(pari.precision(v, 60)) for v in I] + [str(lab)] for I, lab in nodes], open(f"richelot_nodes_{D}{tag}.json", "w"))
    json.dump({"reached": {f"{n}@{k}": s for (n, k), s in seen.items()}, "edges": edges},
              open(f"richelot_graph_{D}{tag}.json", "w"), indent=1)
