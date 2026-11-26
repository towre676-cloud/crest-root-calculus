"""Validate the theta -> Rosenhain -> Igusa pipeline on curves with known invariants (rigorous periods)."""
import sys
sys.path.insert(0, "."); sys.path.insert(0, "../rigorous-search")
import mpmath as mp, cypari2
from periods import thetas, rosenhain, siegel_reduce, symplectic_basis
from periods_arb import period_matrix_arb
pari = cypari2.Pari(); pari.set_real_precision(60)
mp.mp.dps = 60


def tau_of_curve(f):
    S = period_matrix_arb(f, bits=220)
    E = [[int(x) for x in r] for r in S["E"]]
    cols = symplectic_basis(E)
    cv = lambda a: mp.mpf(a.mid().str(60, radius=False))
    Pi = [[mp.mpc(cv(z.real), cv(z.imag)) for z in row] for row in S["Pi"]]
    Z = [[sum(cols[c][k] * Pi[r][k] for k in range(4)) for r in range(2)] for c in range(4)]
    Oa = mp.matrix([[Z[0][0], Z[1][0]], [Z[0][1], Z[1][1]]]); Ob = mp.matrix([[Z[2][0], Z[3][0]], [Z[2][1], Z[3][1]]])
    for tau in (Oa ** -1 * Ob, Ob ** -1 * Oa):
        im = [[tau[i, j].imag for j in range(2)] for i in range(2)]
        if abs(tau[0, 1] - tau[1, 0]) < 1e-40 and im[0][0] > 0 and im[0][0] * im[1][1] - im[0][1] ** 2 > 0:
            return siegel_reduce(tau)
    raise ValueError


def absinv(J):
    J2, J4, J6, J8, J10 = J
    return [J2 ** 5 / J10, J2 ** 3 * J4 / J10, J2 ** 2 * J6 / J10]


def compare(Pstr, rosen=rosenhain):
    P = pari(Pstr)
    f = [int(pari.polcoef(P, k)) for k in range(7)]
    while f and f[-1] == 0:
        f.pop()
    tau = tau_of_curve(f)
    lam, mu, nu = rosen(thetas(tau))
    s = lambda z: "(%s + I*(%s))" % (mp.nstr(z.real, 55), mp.nstr(z.imag, 55))
    Q = pari("x*(x-1)*(x-%s)*(x-%s)*(x-%s)" % (s(lam), s(mu), s(nu)))
    a, b = absinv(pari.genus2igusa(P)), absinv(pari.genus2igusa(Q))
    err = max(float(abs(x - y) / (abs(x) + 1)) for x, y in zip(a, b))
    return a, b, err


if __name__ == "__main__":
    for Pstr in ("x*(x-1)*(x-2)*(x-3)*(x-5)", "x^6 + 3*x^4 - 2*x + 7", "x^5 - x + 1"):
        a, b, err = compare(Pstr)
        print(Pstr, "relative error", err)
        print("   true ", [str(pari.precision(x, 12)) for x in a])
        print("   theta", [str(pari.precision(x, 12)) for x in b])
