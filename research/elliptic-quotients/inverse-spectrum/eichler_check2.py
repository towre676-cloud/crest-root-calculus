"""Which middle term is right: -(-n/p) r(n)  or  +(Dn/p) r(n)?  Test the relation
r(p^2 n) + mid(n) r(n) + p r(n/p^2) = (A_p r)(n) exactly on the closure, all n <= N."""
import sys
import cypari2
from twins import closure_and_matrices
from inverse import mat

pari = cypari2.Pari()


def test(d, l, p, N=40):
    classes, ng, A = closure_and_matrices(d, l, [p])
    n = len(classes)
    R = [[1] + [2 * int(x) for x in pari.qfrep(mat(q), p * p * N)] for q in classes]
    Ap = [[int(A[p][i, j]) for j in range(n)] for i in range(n)]
    res = {}
    for name, mid in (("-(-n/p)", lambda m: -int(pari.kronecker(-m, p))), ("+(Dn/p)", lambda m: int(pari.kronecker(d * m, p)))):
        ok = True
        for m in range(1, N + 1):
            for i in range(n):
                lhs = R[i][p * p * m] + mid(m) * R[i][m] + (p * R[i][m // (p * p)] if m % (p * p) == 0 else 0)
                if lhs != sum(Ap[i][j] * R[j][m] for j in range(n)):
                    ok = False
        res[name] = ok
    return res


for spec in sys.argv[1:]:
    d, l, p = map(int, spec.split(":"))
    print(d, l, p, "(-D/p) =", int(pari.kronecker(-d, p)), test(d, l, p), flush=True)
