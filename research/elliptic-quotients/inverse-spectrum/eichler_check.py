"""Empirical form of the theta/T(p^2) relation on the closure: is r(p^2 n) - eps_p(n) r(n) - p r(n/p^2)
= (A_p r)(n) + c r(n) for a constant c?  Tested on all n <= N for the closure at one level."""
import sys
import cypari2
from twins import closure_and_matrices
from inverse import mat

pari = cypari2.Pari()
d, l, p = map(int, sys.argv[1:4])
classes, ng, A = closure_and_matrices(d, l, [p])
n = len(classes)
N = 60
R = [[1] + [2 * int(x) for x in pari.qfrep(mat(q), p * p * N)] for q in classes]   # r(0..p^2 N), signed counts
Ap = [[int(A[p][i, j]) for j in range(n)] for i in range(n)]
# for each n, LHS_i(n) = r_i(p^2 n) + s * kron(-n, p) r_i(n) + p r_i(n/p^2); fit s, c by least squares on all (i, n)
import itertools
best = None
for s in (-1, 0, 1):
    rows = []
    for m in range(1, N + 1):
        for i in range(n):
            lhs = R[i][p * p * m] + s * int(pari.kronecker(-m, p)) * R[i][m] + (p * R[i][m // (p * p)] if m % (p * p) == 0 else 0)
            rhs = sum(Ap[i][j] * R[j][m] for j in range(n))
            rows.append((lhs - rhs, R[i][m]))
    # need lhs - rhs = c * r_i(m)
    cs = {Fr for Fr in []}
    ratios = set()
    ok = True
    for diff, r in rows:
        if r == 0:
            if diff != 0:
                ok = False; break
        else:
            ratios.add(diff / r)
    print("s", s, "consistent", ok, "c values", sorted(ratios)[:6])
