"""Algebraic inputs and rank witnesses for the nine period-search cases."""
import importlib.util
import json
from fractions import Fraction
from math import isqrt
from pathlib import Path
import sys

import sympy as sp
from flint import acb, acb_mat, arb

from exact_lattices import gaussian, pari, rational


PARENT = Path(__file__).resolve().parent.parent
NS = PARENT / "ns-certificates"
QM = PARENT / "qm-certificate"
sys.path.insert(0, str(NS))

from biell import bielliptic
from load import load_curve
from pipeline import ell_trace, jinv


CASES = {
    "C": (150, 8),
    "Cp": (150, 16),
    "T1": (150, 34),
    "T2": (150, 34),
    "0": (150, 16),
    "1": (150, 32),
    "3": (150, 33),
    "cm_102": (40, 80),
    "QM": (60, 24),
}


def isogeny_multiplier(first, second):
    """Differential multiplier of the dual of a verified degree-three isogeny."""
    As, Bs = (-Fraction(first[key][0]) / 4 for key in ("g2", "g3"))
    Ad, Bd = (-Fraction(second[key][0]) / 4 for key in ("g2", "g3"))
    source = pari.ellinit([0, 0, 0, pari(str(As)), pari(str(Bs))])
    target = pari.ellinit([0, 0, 0, pari(str(Ad)), pari(str(Bd))])
    entries, degrees = pari.ellisomat(source)[:2]
    idx = next(j for j in range(len(entries)) if int(degrees[0, j]) == 3
               and pari.ellinit(entries[j][0]).j() == target.j())
    model, forward = entries[idx][:2]
    A, B = map(lambda t: Fraction(str(t)), model)
    u2 = Bd * A / (B * Ad)
    if u2 ** 2 * A != Ad or u2 ** 3 * B != Bd:
        raise ArithmeticError("Invalid target isomorphism.")
    x = sp.Symbol("x")

    def poly(value):
        value = pari.subst(value, "y", 1)
        return sum(sp.Rational(str(pari.polcoef(value, j, "x"))) * x ** j
                   for j in range(int(pari.poldegree(value, "x")) + 1))

    N, M, D = (poly(value) for value in forward)
    U, V = sp.cancel(N / D ** 2), sp.cancel(M / D ** 3)
    n, d = sp.fraction(U)
    source_cubic = x ** 3 + As * x + Bs
    if (max(sp.degree(n, x), sp.degree(d, x)) != 3 or sp.degree(n, x) <= sp.degree(d, x)
            or sp.cancel(V ** 2 * source_cubic - U ** 3 - A * U - B) != 0):
        raise ArithmeticError("Invalid explicit degree-three isogeny.")
    c = sp.cancel(U.diff(x) / V)
    if not c.is_Rational or not c:
        raise ArithmeticError("Invalid isogeny differential.")
    return acb(rational(u2)).sqrt() * rational(3 / c)


def rank_bound(recs, evidence):
    """Verify the split endomorphism rank bound from elliptic reductions."""
    left, right = [recs[i] for i in evidence["pair"]]

    def invariants(m):
        return m["g2"], m["g3"]

    def trace(m, p):
        if not sp.isprime(p) or p < 5:
            raise ArithmeticError("Invalid reduction prime.")
        g2, g3 = invariants(m)
        i = next((t for t in range(1, p) if t*t % p == p-1), None)

        def mod(z):
            a, b = map(Fraction, z)
            if b and i is None:
                raise ArithmeticError("The coefficient field does not split.")
            return (a.numerator * pow(a.denominator, -1, p)
                    + (i or 0) * b.numerator * pow(b.denominator, -1, p)) % p

        if (mod(g2) ** 3 - 27 * mod(g3) ** 2) % p == 0:
            raise ArithmeticError("Bad elliptic reduction.")
        a = ell_trace((g2, g3), p, i)
        if a % p == 0:
            raise ArithmeticError("The witness requires ordinary reduction.")
        return a, int(pari.coredisc(a*a - 4*p))

    if "non_isogeny_witness" in evidence:
        p = evidence["non_isogeny_witness"]["p"]
        if trace(left, p)[1] == trace(right, p)[1]:
            raise ArithmeticError("The two geometric reduction fields coincide.")
        for m in (left, right):
            j = jinv(*invariants(m))
            if not isinstance(j, Fraction) or j.denominator == 1:
                raise ArithmeticError("The witness does not exclude CM.")
        return 2
    if "non_cm_witness" in evidence:
        fields = {trace(left, p)[1] for p in evidence["non_cm_witness"].values()}
        if len(fields) < 2:
            raise ArithmeticError("The two ordinary reductions do not exclude CM.")
        return 4
    if jinv(*invariants(left)) != 1728 or jinv(*invariants(right)) != 1728:
        raise ArithmeticError("Unexpected CM case.")
    return 8


def split_case(name):
    receipt = json.loads((NS / f"cert_{name}.json").read_text(encoding="utf-8"))
    factor = None
    if name.isdecimal():
        cases = json.loads((PARENT / "bielliptic-reconstruction/test_curves.json").read_text(encoding="utf-8"))
        f, recs = bielliptic(cases[int(name)]["unscrambled_c"])
        factor = isogeny_multiplier(*recs)
    elif name == "cm_102":
        f, recs = bielliptic(receipt["c"])
        ratio = recs[1]["g2"][0] / recs[0]["g2"][0]
        u2 = Fraction(isqrt(ratio.numerator), isqrt(ratio.denominator))
        if u2 ** 2 != ratio or any(m["g3"] != (0, 0) for m in recs):
            raise ArithmeticError("Invalid exact CM isomorphism.")
        factor = acb(rational(u2)).sqrt()
    else:
        _, f, recs = load_curve({"C": "C (", "Cp": "C'"}.get(name, name))
    if f != receipt["f"]:
        raise ArithmeticError("The equation differs from the rank-witness model.")
    p, q = (recs[i] for i in receipt["rank_certificate"]["pair"])
    a, b = (sp.Rational(row[0]) + sp.I * sp.Rational(row[1]) for row in (p["a"], p["b"]))
    c, d = (sp.Rational(row[0]) + sp.I * sp.Rational(row[1]) for row in (q["a"], q["b"]))
    if sp.expand(a*d - b*c) == 0:
        raise ArithmeticError("The maps do not give an isogeny with a product.")
    rank = rank_bound(recs, receipt["rank_certificate"])
    diffs = [(gaussian(r["a"]), gaussian(r["b"]), r["degree"]) for r in recs]
    if factor is not None:
        a, b, _ = diffs[0]
        c, d, _ = diffs[1]
        mults = (factor, -acb(0, 1) * factor) if name == "cm_102" else (factor,)
        degree = 4 if name == "cm_102" else 8
        # The two complementary involutions give orthogonal targets. Hence
        # deg(f1 + dual(rho) f2) = 2 + 2 deg(rho), also for the CM automorphism.
        diffs.extend((a + m*c, b + m*d, degree) for m in mults)
    return f, diffs, rank, receipt


def qm_tangents():
    sys.path.insert(0, str(QM))
    spec = importlib.util.spec_from_file_location("qm_certificate", QM / "certificate.py")
    module = importlib.util.module_from_spec(spec)
    spec.loader.exec_module(module)
    ma, mg = module.tangent("alpha"), module.tangent("gamma")
    identity = pari.matid(2)
    beta = ma * mg - mg * ma
    if (pari.lift(ma * ma) != 2 * identity or pari.lift(mg * mg + mg - identity) != 0 * identity
            or pari.lift(beta * beta) != -6 * identity or pari.lift(ma * beta + beta * ma) != 0 * identity
            or [int(pari.hilbert(2, -6, p)) for p in (2, 3)] != [-1, -1]):
        raise ArithmeticError("The exact correspondences do not generate the claimed division algebra.")
    left, right = (module.geometric_frobenius_field(p) for p in (7, 11))
    if not left or not right or left[0] == right[0]:
        raise ArithmeticError("The ordinary reductions do not rule out a CM square.")

    def embed(M, k):
        t = acb((arb(k) / 4).cos_pi(), (arb(k) / 4).sin_pi())
        return acb_mat([[sum((rational(str(pari.polcoef(pari.lift(M[i, j]), e, "t"))) * t ** e
                             for e in range(4)), acb(0)) for j in range(2)] for i in range(2)])

    return module.F, [[embed(M, k) for M in (ma, mg)] for k in (1, 3, 5, 7)]
