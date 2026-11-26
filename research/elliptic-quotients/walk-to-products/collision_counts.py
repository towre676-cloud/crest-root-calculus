"""For each collision pair: compare r(m^2) at coprime m and at all m (m <= 400), the full theta series, and
determinants; decide whether the coprime agreement is a theorem (Sturm bound of the t=1 Shimura lift)."""
import json, math, sys
sys.path.insert(0, "../inverse-spectrum")
import cypari2
from inverse import mat, level, sturm

pari = cypari2.Pari()
pari.allocatemem(8 * 10 ** 9)
C = json.load(open("collisions.json"))
out = []
for r in C:
    d, l = r["d"], r["l"]
    bad = 2 * l * abs(d)
    for grp in r["collision_groups"]:
        q1, q2 = grp[0]["q"], grp[1]["q"]
        N = level(q1)
        B = max(sturm(N), 400)
        B = min(B, 1200)
        r1 = pari.qfrep(mat(q1), B * B); r2 = pari.qfrep(mat(q2), B * B)
        cop = [m for m in range(1, B + 1) if math.gcd(m, bad) == 1]
        diff_cop = next((m for m in cop if r1[m * m - 1] != r2[m * m - 1]), None)
        diff_sq = next((m for m in range(1, B + 1) if r1[m * m - 1] != r2[m * m - 1]), None)
        diff_all = next((k for k in range(1, 2000) if r1[k - 1] != r2[k - 1]), None)
        out.append({"d": d, "l": l, "q1": q1, "q2": q2, "det": [int(pari.matdet(mat(q1))), int(pari.matdet(mat(q2)))],
                    "sturm_bound": sturm(N), "checked_to": B, "first_coprime_square_difference": diff_cop,
                    "first_square_difference": diff_sq, "first_theta_difference": diff_all,
                    "coprime_equal_proved": diff_cop is None and B >= sturm(N)})
        print(out[-1], flush=True)
json.dump(out, open("collision_counts.json", "w"), indent=1)
