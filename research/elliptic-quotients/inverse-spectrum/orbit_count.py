"""Number of <units, conjugation>-orbits on {a in F_{l^2} : N(a) != -1}, compared with isometry classes."""
import json, glob
from scalar_part import ops
UNITS = {-3: 6, -4: 4}


def orbits(d, l):
    mul, conj, _, _ = ops(d, l)
    c = (d * d - d) // 4
    # unit group generators in O_K = Z[w]: -1 always; i = w + 2 for d = -4; zeta_6 = -w ... (w^2 = -3w - 3 for d=-3)
    units = [(1, 0), (l - 1, 0)]
    if d == -4:
        units = [(1, 0), (2 % l, 1), (l - 1, 0), ((-2) % l, l - 1)]
    if d == -3:
        z = (0, 1)                 # w = (-3 + sqrt -3)/2 has w^2 + 3w + 3 = 0; w + 1 is a 6th root of unity
        z = (2, 1)                 # zeta_6 = w + 2
        units, u = [], (1, 0)
        for _ in range(6):
            units.append(u); u = mul(u, z)
    seen, n = set(), 0
    for p in range(l):
        for q in range(l):
            a = (p, q)
            N = mul(a, conj(a))[0]
            if (N + 1) % l == 0 or a in seen:
                continue
            n += 1
            for u in units:
                for x in (mul(u, a), conj(mul(u, a))):
                    seen.add(x)
    return n


rows = []
for f in sorted(glob.glob("inert_*.json")):
    for r in json.load(open(f)):
        rows.append((r["d"], r["l"], r["isometry_classes"], orbits(r["d"], r["l"]),
                     [(g["size"], g["max_witness_degree"], len(g["square_equivalent_pairs"])) for g in r["genera"]]))
for row in sorted(set(map(lambda x: (x[0], x[1], x[2], x[3], str(x[4])), rows))):
    print(row)
