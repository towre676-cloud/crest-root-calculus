"""Multiplicative Hecke formula for coprime quotient counts.

From r(p^2 n) + (Dn/p) r(n) + p r(n/p^2) = (A_p r)(n) on the p-neighbour closure:
  for p not dividing m':  v(p m') = (A_p - (D/p)) v(m')         [(D p^0 m'^2 / p) = (D/p)]
  for k >= 1:             v(p^{k+1} m') = A_p v(p^k m') - p v(p^{k-1} m')
so for m coprime to the level,  v(m) := (r_x(m^2))_x = prod_{p^k || m} Q_{p,k}(A_p) v_0,
Q_{p,0} = 1, Q_{p,1} = A_p - (D/p), Q_{p,k+1} = A_p Q_{p,k} - p Q_{p,k-1}, v_0 = (r_x(1))_x.
r_x(1) counts the vectors of Humbert norm 1, i.e. degree-1 elliptic subvarieties: v_0 is supported on
the product (reducible) surfaces of the closure.  Verified here exactly against qfrep."""
import itertools, json, math, sys
sys.path.insert(0, "../inverse-spectrum")
import cypari2
from twins import closure_and_matrices
from inverse import mat

pari = cypari2.Pari()
pari.allocatemem(8 * 10 ** 9)


def Q(Ap, p, k, D, n):
    I = pari.matid(n)
    q0, q1 = I, Ap - I * int(pari.kronecker(D, p))
    if k == 0:
        return q0
    for _ in range(k - 1):
        q0, q1 = q1, Ap * q1 - q0 * p
    return q1


def check(d, l, primes, M=1100):
    classes, ng, A = closure_and_matrices(d, l, primes)
    n = len(classes)
    reps = [pari.qfrep(mat(q), M * M) for q in classes]
    r = lambda x, N: 2 * int(reps[x][N - 1])
    v0 = pari.matrix(n, 1, [r(x, 1) for x in range(n)])
    tested, ok = 0, True
    ms = sorted({math.prod(c) for k in range(1, 4) for c in itertools.combinations_with_replacement(primes, k)})
    for m in ms:
        if m > M:
            continue
        fac = pari.factor(m)
        P = pari.matid(n)
        for i in range(int(pari.matsize(fac)[0])):
            p, k = int(fac[i, 0]), int(fac[i, 1])
            P = Q(A[p], p, k, d, n) * P
        pred = P * v0
        for x in range(n):
            tested += 1
            if int(pred[x, 0]) != r(x, m * m):
                ok = False
    products = [x for x in range(n) if r(x, 1) > 0]
    return {"d": d, "l": l, "primes": primes, "closure": n, "gluing": ng, "m_tested": [m for m in ms if m <= M],
            "entries_tested": tested, "formula_exact": ok, "product_classes": len(products),
            "r1_on_products": [r(x, 1) for x in products],
            "product_classes_among_gluing": sum(1 for x in products if x < ng)}


if __name__ == "__main__":
    out = []
    for spec in sys.argv[1:]:
        d, l, ps = spec.split(":")
        res = check(int(d), int(l), [int(x) for x in ps.split(",")])
        out.append(res)
        print(res, flush=True)
    json.dump(out, open("../walk-to-products/product_formula.json", "w"), indent=1)
