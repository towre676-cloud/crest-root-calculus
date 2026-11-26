"""Galois orbits of CM polarised surfaces by Shimura reciprocity, and exact invariant polynomials.

For (Lambda = a1 + a2, H) with O_3 acting diagonally, the conjugates are (b Lambda, H / N b) for b over the proper
ideal classes of O_3, together with complex conjugates.  The orbit of each absolute Igusa invariant gives
prod (X - z_i); its coefficients are recognised as rationals (bestappr) and accepted only when the rounding error is
at least 10^-(DPS/2) relative AND every orbit value is a root to working precision.  The other two invariants are
expressed as Q-polynomials in z by exact interpolation over the orbit (Galois-equivariant, so rational)."""
import itertools, json, math, os, sys
sys.path.insert(0, "."); sys.path.insert(0, "../inverse-spectrum")
from fractions import Fraction as Fr
import mpmath as mp
import cypari2
from gluing import K, E_of
from polarised import ideals
from periods import period_matrix, thetas, rosenhain, siegel_reduce

DPS = int(os.environ.get("TWIN_DPS", "500"))
pari = cypari2.Pari(); pari.allocatemem(8 * 10 ** 9); pari.set_real_precision(DPS - 20)
mp.mp.dps = DPS


def lattice_mult(Kf, b, a):
    gens = [Kf.mul(x, y) for x in b for y in a]
    den = 1
    for g in gens:
        for c in g:
            den = den * c.denominator // math.gcd(den, c.denominator)
    M = pari.matrix(2, len(gens), [int(gens[j][i] * den) for i in range(2) for j in range(len(gens))])
    Hm = pari.mathnf(M)
    return [(Fr(int(Hm[0, j]), den), Fr(int(Hm[1, j]), den)) for j in range(2)]


def covol(Kf, a):
    return abs(Kf.mul(a[0], Kf.conj(a[1]))[1])


def invariants(d, B, H):
    tau = siegel_reduce(period_matrix(d, B, H))
    lam, mu, nu = rosenhain(thetas(tau))
    s = lambda z: "(%s + I*(%s))" % (mp.nstr(z.real, DPS - 15), mp.nstr(z.imag, DPS - 15))
    J2, J4, J6, J8, J10 = pari.genus2igusa(pari("x*(x-1)*(x-%s)*(x-%s)*(x-%s)" % (s(lam), s(mu), s(nu))))
    return [J2 ** 5 / J10, J2 ** 3 * J4 / J10, J2 ** 2 * J6 / J10]


def orbit(d, rec, f=3):
    Kf = K(d)
    f1, a1, f2, a2 = rec["product"]
    a1 = [tuple(Fr(x) for x in e) for e in a1]; a2 = [tuple(Fr(x) for x in e) for e in a2]
    H = tuple(Fr(x) for x in rec["H"])
    vals = []
    for b in ideals(d, f):
        b = [tuple(e) for e in b]
        Nb = covol(Kf, b) / covol(Kf, [(Fr(1), Fr(0)), (Fr(0), Fr(f))])      # index of b in O_f
        A1, A2 = lattice_mult(Kf, b, a1), lattice_mult(Kf, b, a2)
        B = [(A1[0][0], A1[0][1], Fr(0), Fr(0)), (A1[1][0], A1[1][1], Fr(0), Fr(0)),
             (Fr(0), Fr(0), A2[0][0], A2[0][1]), (Fr(0), Fr(0), A2[1][0], A2[1][1])]
        Hb = tuple(x / Nb for x in H)
        E = [[E_of(Kf, Hb, B[i], B[j]) for j in range(4)] for i in range(4)]
        assert all(x.denominator == 1 for r in E for x in r)
        assert abs(int(pari.matdet(pari.matrix(4, 4, [int(x) for r in E for x in r])))) == 1
        vals.append(invariants(d, B, Hb))
    return vals


def distinct(zs, tol):
    out = []
    for z in zs:
        if all(abs(z - w) > tol * (1 + abs(z)) for w in out):
            out.append(z)
    return out


def rational_poly(roots):
    P = pari(1)
    for r in roots:
        P = P * (pari("x") - r)
    n = len(roots)
    coeffs, worst = [], pari(0)
    for k in range(n + 1):
        c = pari.polcoef(P, k)
        re = pari.real(c)
        q = pari.bestappr(re, pari("10^%d" % (DPS // 3)))
        err = abs(c - q) / (1 + abs(c))
        worst = max(worst, err)
        coeffs.append(q)
    Q = sum((coeffs[k] * pari("x") ** k for k in range(n + 1)), pari(0))
    return Q, worst


if __name__ == "__main__":
    d = int(sys.argv[1])
    L = json.load(open(f"twin_lattices_{d}.json"))
    names = sys.argv[2].split(",") if len(sys.argv) > 2 else sorted(L)
    out = {}
    for name in names:
        vals = orbit(d, L[name])
        z = [v[0] for v in vals]
        zs = distinct(z + [pari.conj(x) for x in z], pari("10^-%d" % (DPS // 2)))
        Q, worst = rational_poly(zs)
        red = pari.polredbest(Q) if int(pari.poldegree(Q)) <= 24 else None
        irr = bool(pari.polisirreducible(Q))
        out[name] = {"orbit_size": len(zs), "minpoly_z": str(Q), "irreducible": irr, "rounding_error": float(pari.log(worst + pari("10^-%d" % (DPS * 2)))) / 2.302585,
                     "field_polredbest": str(red), "field_disc": str(pari.factor(pari.nfdisc(red))) if red is not None else None}
        print(name, "orbit", len(zs), "| irreducible", irr, "| log10 rounding error %.0f" % out[name]["rounding_error"],
              "| max coeff digits", max(len(str(abs(int(pari.numerator(pari.polcoef(Q, k)))))) for k in range(int(pari.poldegree(Q)) + 1)),
              "| field", red, "| disc", out[name]["field_disc"], flush=True)
    json.dump(out, open(f"orbit_{d}.json", "w"), indent=1)
