"""Verify the stored algebraic correspondences and certify QM absence."""
import json
from pathlib import Path
import sys

import cypari2

from verify_morphism import verify

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "exact-completeness"))
from absence import good_model, power_polynomial, weil_polynomial  # noqa: E402

pari = cypari2.Pari()
F = [-1, 5, -8, 4, -1, 1]


def tangent(which):
    data = json.loads((HERE / f"functions_{which}.json").read_text(encoding="utf-8"))
    functions = [f"({data[key]['N']})/({data[key]['D']})" for key in ("s1", "s0", "g1", "g0")]
    verified = verify(F, *functions, "t^4+1")
    if not (verified["member"] and verified["tangent_linear"]):
        raise ArithmeticError(f"Invalid correspondence for {which}.")
    entries = [pari(entry) for row in verified["tangent"] for entry in row]
    return pari.matrix(2, 2, entries)


def geometric_frobenius_field(p, rmax=12):
    if not good_model(F, p):
        return None
    P = weil_polynomial(F, p)
    for r in range(1, rmax + 1):
        Q = pari.Pol([int(c) for c in power_polynomial(P, p, r).all_coeffs()])
        fac = pari.factor(Q)
        if len(fac[0]) == 1 and int(fac[1][0]) == 2 and pari.poldegree(fac[0][0]) == 2:
            q2 = fac[0][0]
            a = -int(pari.polcoef(q2, 1))
            if a % p == 0:
                return None   # supersingular
            return int(pari.coredisc(a * a - 4 * p ** r)), r
    return None


if __name__ == "__main__":
    Ma, Mg = tangent("alpha"), tangent("gamma")
    I = pari.matid(2)
    Mb = Ma * Mg - Mg * Ma            # anticommutes with alpha whenever alpha^2 is scalar
    Mb2 = pari.lift(Mb * Mb)
    b = Mb2[0, 0]
    rel = {"alpha^2 = 2": pari.lift(Ma * Ma) == 2 * I,
           "gamma^2 + gamma - 1 = 0": pari.lift(Mg * Mg + Mg - I) == 0 * I,
           "gamma not in Q(alpha)": Mb2 != 0 * I,
           "beta = alpha gamma - gamma alpha has rational square": Mb2 == b * I and str(pari.type(b)) in ("t_INT", "t_FRAC"),
           "alpha beta = - beta alpha": pari.lift(Ma * Mb + Mb * Ma) == 0 * I}
    b = int(b) if str(pari.type(b)) == "t_INT" else b
    hilbert_all = {int(p): int(pari.hilbert(2, b, p)) for p in [2, 3, 5, 7] + [int(q) for q in pari.factor(abs(pari.numerator(b)))[0]]}
    hilbert = -1 if -1 in hilbert_all.values() else 1
    fields = {}
    for p in (7, 11):
        gf = geometric_frobenius_field(p)
        if gf is not None:
            fields[gf[0]] = (p, gf[1])
    out = {"relations": {k: bool(v) for k, v in rel.items()}, "beta_squared": str(b),
           "hilbert_symbols_(2,b)_p": hilbert_all,
           "division_algebra": hilbert == -1, "frobenius_fields": {str(k): v for k, v in fields.items()},
           "not_CM_square": len(fields) >= 2}
    out["no_elliptic_quotient"] = all(out["relations"].values()) and out["division_algebra"] and out["not_CM_square"]
    print(json.dumps(out, indent=1))
    receipt = json.loads((HERE / "certificate.json").read_text(encoding="utf-8"))
    if not out["no_elliptic_quotient"] or json.loads(json.dumps(out)) != receipt:
        raise SystemExit("The exact QM replay did not match its certificate.")
