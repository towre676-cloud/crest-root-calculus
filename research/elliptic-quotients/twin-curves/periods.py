"""Period matrix, theta constants and Rosenhain branch points of a principally polarised (Lambda, H), Lambda in K^2."""
import itertools, sys
sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
import mpmath as mp
from gluing import K, E_of

mp.mp.dps = 280


def symplectic_basis(E):
    """Columns M (integer 4x4) with M^T E M = [[0, I], [-I, 0]] for unimodular alternating integer E."""
    n = 4
    vecs = [[1 if i == j else 0 for i in range(n)] for j in range(n)]
    form = lambda u, v: sum(u[i] * E[i][j] * v[j] for i in range(n) for j in range(n))
    es, fs = [], []
    while vecs:
        # pick e with some partner
        e = vecs[0]
        vals = [form(e, v) for v in vecs]
        # extended gcd to find f with form(e, f) = 1
        import math
        g, coeffs = 0, [0] * len(vecs)
        for k, val in enumerate(vals):
            if val == 0:
                continue
            if g == 0:
                g, coeffs = val, [0] * len(vecs); coeffs[k] = 1
                continue
            # combine g and val
            a, b = g, val
            x0, x1, y0, y1 = 1, 0, 0, 1
            while b:
                qq = a // b
                a, b = b, a - qq * b
                x0, x1 = x1, x0 - qq * x1
                y0, y1 = y1, y0 - qq * y1
            coeffs = [c * x0 for c in coeffs]; coeffs[k] += y0
            g = a
        if g < 0:
            coeffs = [-c for c in coeffs]; g = -g
        assert g == 1, "not unimodular"
        f = [sum(coeffs[k] * vecs[k][i] for k in range(len(vecs))) for i in range(n)]
        assert form(e, f) == 1
        es.append(e); fs.append(f)
        # project the rest onto the orthogonal complement of <e, f>
        rest = []
        for v in vecs[1:]:
            w = [v[i] - form(v, f) * e[i] + form(v, e) * f[i] for i in range(n)]
            rest.append(w)
        # basis of the span of rest (rank len(vecs) - 2) via HNF
        import cypari2
        pari = cypari2.Pari()
        if len(vecs) - 2 <= 0:
            vecs = []
            continue
        M = pari.matrix(n, len(rest), [rest[j][i] for i in range(n) for j in range(len(rest))])
        Hm = pari.mathnf(M)
        vecs = [[int(Hm[i, j]) for i in range(n)] for j in range(int(pari.matsize(Hm)[1]))]
    cols = es + fs
    return cols


def period_matrix(d, B, H):
    Kf = K(d)
    E = [[int(E_of(Kf, H, B[i], B[j])) for j in range(4)] for i in range(4)]
    cols = symplectic_basis(E)
    w = (d + mp.sqrt(mp.mpf(d))) / 2 if d < 0 else None
    w = (mp.mpf(d) + mp.mpc(0, mp.sqrt(-d))) / 2
    def z(col):
        v = [sum(Fr(col[k]) * B[k][t] for k in range(4)) for t in range(4)]
        return [mp.mpf(v[0].numerator) / v[0].denominator + (mp.mpf(v[1].numerator) / v[1].denominator) * w,
                mp.mpf(v[2].numerator) / v[2].denominator + (mp.mpf(v[3].numerator) / v[3].denominator) * w]
    Z = [z(c) for c in cols]
    Oa = mp.matrix([[Z[0][0], Z[1][0]], [Z[0][1], Z[1][1]]])
    Ob = mp.matrix([[Z[2][0], Z[3][0]], [Z[2][1], Z[3][1]]])
    for tau in (Oa ** -1 * Ob, Ob ** -1 * Oa, -(Oa ** -1 * Ob), -(Ob ** -1 * Oa)):
        sym = abs(tau[0, 1] - tau[1, 0]) < mp.mpf(10) ** (-200)
        im = mp.matrix([[tau[i, j].imag for j in range(2)] for i in range(2)])
        if sym and im[0, 0] > 0 and im[0, 0] * im[1, 1] - im[0, 1] * im[1, 0] > 0:
            return tau
    raise ValueError("no Siegel normalisation")


def theta(a, b, tau, N=None):
    """theta[a;b](0, tau), characteristic (a/2, b/2), a, b in {0,1}^2."""
    im = mp.matrix([[tau[i, j].imag for j in range(2)] for i in range(2)])
    lam = min(mp.eig(im)[0])
    if N is None:
        N = int(mp.sqrt(mp.mp.dps * 2.4 / (mp.pi * lam))) + 3
    s = mp.mpc(0)
    for n1 in range(-N, N + 1):
        for n2 in range(-N, N + 1):
            x = (n1 + mp.mpf(a[0]) / 2, n2 + mp.mpf(a[1]) / 2)
            q = x[0] * x[0] * tau[0, 0] + 2 * x[0] * x[1] * tau[0, 1] + x[1] * x[1] * tau[1, 1]
            s += mp.exp(mp.pi * 1j * q + mp.pi * 1j * (x[0] * b[0] + x[1] * b[1]))
    return s


def thetas(tau):
    out = {}
    for i in range(16):
        a = (i & 1, (i >> 1) & 1); b = ((i >> 2) & 1, (i >> 3) & 1)
        if (a[0] * b[0] + a[1] * b[1]) % 2 == 0:
            out[i] = theta(a, b, tau)
    return out


def rosenhain(th):
    t2 = {i: v ** 2 for i, v in th.items()}
    lam = t2[0] * t2[2] / (t2[1] * t2[3])
    mu = t2[2] * t2[8] / (t2[3] * t2[9])
    nu = t2[0] * t2[8] / (t2[1] * t2[9])
    return lam, mu, nu


def act(M, tau):
    """(A tau + B)(C tau + D)^{-1}, M = [[A, B], [C, D]] 4x4 integer."""
    A = mp.matrix([[M[0][0], M[0][1]], [M[1][0], M[1][1]]]); Bm = mp.matrix([[M[0][2], M[0][3]], [M[1][2], M[1][3]]])
    C = mp.matrix([[M[2][0], M[2][1]], [M[3][0], M[3][1]]]); D = mp.matrix([[M[2][2], M[2][3]], [M[3][2], M[3][3]]])
    return (A * tau + Bm) * (C * tau + D) ** -1


def siegel_reduce(tau, maxit=200):
    for _ in range(maxit):
        # Gauss-reduce Im tau with U, tau -> U tau U^T
        Y = [[tau[i, j].imag for j in range(2)] for i in range(2)]
        U = [[1, 0], [0, 1]]
        a, b, c = Y[0][0], Y[0][1], Y[1][1]
        changed = True
        while changed:
            changed = False
            if c < a:
                a, c = c, a; U = [U[1], U[0]]; changed = True
            k = int(mp.nint(b / a))
            if k:
                b = b - k * a; c = c - 2 * k * (b + k * a) + k * k * a if False else None
                U = [U[0], [U[1][0] - k * U[0][0], U[1][1] - k * U[0][1]]]
                Ym = mp.matrix(Y)
                Um = mp.matrix(U)
                Yn = Um * Ym * Um.T
                a, b, c = Yn[0, 0], Yn[0, 1], Yn[1, 1]
                changed = True
        Um = mp.matrix(U)
        tau = Um * tau * Um.T
        # translate real parts
        for i in range(2):
            for j in range(i, 2):
                r = mp.nint(tau[i, j].real)
                tau[i, j] -= r
                if i != j:
                    tau[j, i] -= r
        if abs(tau[0, 0]) >= 1 - mp.mpf(10) ** -30:
            return tau
        M = [[0, 0, -1, 0], [0, 1, 0, 0], [1, 0, 0, 0], [0, 0, 0, 1]]
        tau = act(M, tau)
    raise ValueError("reduction did not converge")
