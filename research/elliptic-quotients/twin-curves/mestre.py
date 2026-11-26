"""Mestre's reconstruction through transvectants (binary sextic f in x, y).

i = (f,f)_4, y1 = (f,i)_4, y2 = (i,y1)_2, y3 = (i,y2)_2 (quadratic covariants of degrees 3, 5, 7 in the coefficients),
A_ij = (y_i, y_j)_2, A_ijk = (f, y_i y_j y_k)_6, I2 = (f,f)_6.  Weight-zero normalisation
    a_ij = A_ij / I2^((d_i + d_j)/2),   a_ijk = A_ijk / I2^((1 + d_i + d_j + d_k)/2)
gives absolute invariants (independent of the model, hence Galois-equivariant).  Mestre: the curve is
    y^2 = M(q(t)),  L(q(t)) = 0,  L = sum a_ij X_i X_j,  M = sum a_ijk X_i X_j X_k,
for any parametrisation q of the conic L (up to a quadratic twist from the normalisation)."""
import itertools, math
import cypari2

pari = cypari2.Pari()
X, Y = pari("x"), pari("y")


def deriv(P, i, j):
    for _ in range(i):
        P = pari.deriv(P, "x")
    for _ in range(j):
        P = pari.deriv(P, "y")
    return P


def tv(Fp, m, G, n, k):
    s = pari(0)
    for i in range(k + 1):
        s += (-1) ** i * math.comb(k, i) * deriv(Fp, k - i, i) * deriv(G, i, k - i)
    return s * pari(math.factorial(m - k) * math.factorial(n - k)) / (math.factorial(m) * math.factorial(n))


DEG = {1: 3, 2: 5, 3: 7}


def invariants(f):
    """f homogeneous sextic in x, y.  Returns (a_ij dict, a_ijk dict) normalised."""
    i = tv(f, 6, f, 6, 4)                     # quartic
    y1 = tv(f, 6, i, 4, 4)                    # quadratic
    y2 = tv(i, 4, y1, 2, 2)
    y3 = tv(i, 4, y2, 2, 2)
    ys = {1: y1, 2: y2, 3: y3}
    I2 = tv(f, 6, f, 6, 6)
    a2 = {}
    for p, q in itertools.combinations_with_replacement((1, 2, 3), 2):
        a2[(p, q)] = tv(ys[p], 2, ys[q], 2, 2) / I2 ** ((DEG[p] + DEG[q]) // 2)
    a3 = {}
    for p, q, r in itertools.combinations_with_replacement((1, 2, 3), 3):
        a3[(p, q, r)] = tv(f, 6, ys[p] * ys[q] * ys[r], 6, 6) / I2 ** ((1 + DEG[p] + DEG[q] + DEG[r]) // 2)
    return a2, a3, I2


def conic_cubic(a2, a3, T):
    """Symmetric forms over the ring of the values (T: variables list [X1,X2,X3])."""
    L = pari(0)
    for (p, q), c in a2.items():
        L += (1 if p == q else 2) * c * T[p - 1] * T[q - 1]
    M = pari(0)
    for (p, q, r), c in a3.items():
        mult = math.factorial(3) // math.prod(math.factorial(list((p, q, r)).count(s)) for s in set((p, q, r)))
        M += mult * c * T[p - 1] * T[q - 1] * T[r - 1]
    return L, M


def parametrise(L, P0, T, t):
    """Quadratic parametrisation of the conic L through the point P0 (lines through P0 with direction (1, t, 0) or similar)."""
    # points P0 + s * D(t), D(t) = (1, t, 0) unless degenerate; L(P0 + s D) = s (2 B(P0, D) + s L(D)) -> s = -2 B(P0, D) / L(D)
    # direction family R0 + t R1 with (P0, R0, R1) linearly independent
    E = [[1, 0, 0], [0, 1, 0], [0, 0, 1]]
    for a, b in ((0, 1), (0, 2), (1, 2)):
        Mx = pari.matrix(3, 3, [P0[0], E[a][0], E[b][0], P0[1], E[a][1], E[b][1], P0[2], E[a][2], E[b][2]])
        if pari.matdet(Mx) != 0:
            break
    D = [E[a][k] + t * E[b][k] for k in range(3)]
    LD = pari.substvec(L, ["X1", "X2", "X3"], D)
    grad = [pari.deriv(L, v) for v in ("X1", "X2", "X3")]
    BPD = sum((pari.substvec(grad[k], ["X1", "X2", "X3"], P0) * D[k] for k in range(3)), pari(0)) / 2
    # homogeneous parametrisation: LD * P0 - 2 BPD * D
    return [LD * P0[k] - 2 * BPD * D[k] for k in range(3)]
