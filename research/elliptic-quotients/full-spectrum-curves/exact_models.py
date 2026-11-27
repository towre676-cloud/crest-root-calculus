"""Exact genus-two equations for the determinant-1440 full-spectrum pair, with two exact degree-three covers each.

Input: recognised_<A|B>.json (discovery, gluing33.py).  Every assertion below is exact arithmetic in the ring class
field H = Q(theta, sqrt(-10)) of the order O_3 = Z[3 sqrt(-10)], theta a root of the polredabs model F of H_{-360},
j1 = phi(theta).  For each curve:

  E1 : Y^2 = g1(X) = X^3 + A X + B,  A = 3 j1 (1728 - j1),  B = 2 j1 (1728 - j1)^2,  H_{-360}(j1) = 0;
  Q  = x^3 + x^2 + q1 x + q0,  f1 = v + kappa x / Q;
  Q^3 g1(f1) = (v Q + kappa x)^3 + A (v Q + kappa x) Q^2 + B Q^3 = c1 R S1^2   with R cubic, S1 cubic;
  C  : y^2 = Q R  (squarefree of degree 6, so genus two), phi1 = (f1, sqrt(c1) S1 y / Q^2) : C -> E1 of degree 3;
  f2 = P2 / R with P2 = x^2 + p1 x + p0, g2(X) = Res_x(Q, X R - P2) / Res_x(Q, R),
  R^3 g2(f2) = c2 Q S2^2,  phi2 : C -> E2 : Y^2 = g2(X) of degree 3,  H_{-360}(j(E2)) = 0.
Hence both quotients have CM by exactly O_3, the two degree-3 maps are distinct (different triple fibres), and
chi = phi1^* + phi2^* : E1 x E2 -> J(C) is an isogeny.  With psi = (phi1_*, phi2_*), psi chi = [[3, b], [b', 3]], whose
adjugate product is (9 - deg b) I, so ker chi is killed by some e in {1, ..., 9}, e | 2520.  The scalar 3 sqrt(-10)
in O_3 therefore gives the integral homology matrix 7560 sqrt(10) J on J(C) (up to sign), and its integral
centralizer is the full geometric endomorphism ring (mestre_link.certify_cm_order).

Usage: python exact_models.py --verify [A B] [--embeddings=0,5]
Without --verify the program regenerates the two receipts.
"""
import argparse
import json
import sys
from pathlib import Path

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent / "twin-curves"))
sys.path.insert(0, str(HERE.parent / "inverse-spectrum"))
from mestre_link import normalised_quintic, period_matrix_acb, certify_cm_order, humbert_from_order, pari
from gluing import isometric
from igusa_exact import absolute_invariants

pari.allocatemem(2 * 10 ** 9)


def require(c, msg):
    if not c:
        raise AssertionError(msg)


def field():
    Hp = pari.polclass(-360)
    F, a = pari.polredabs(Hp, 1)
    Fy = F.subst("x", "y")
    nf = pari.nfinit(Fy)
    Pabs, alpha, k = pari.rnfequation(nf, pari("x^2 + 10"), 1)       # x = s + k*theta, s^2 = -10
    Pw = Pabs.subst("x", "y")
    require(Pw.poldegree() == 16 and pari.polisirreducible(Pw), "degree-16 coefficient field")
    th = pari.Mod(pari.lift(alpha).subst("x", "y"), Pw)                # theta in the absolute field
    s = pari.Mod(pari("y"), Pw) - int(k) * th
    require(s ** 2 == -10, "sqrt(-10)")
    phi = pari.lift(a)
    j1 = pari.subst(phi, "x", th)
    require(pari.subst(Hp, "x", j1) == 0, "j1 is a root of H_{-360}")
    return {"F": F, "Pabs": Pw, "theta": th, "s": s, "j1": j1, "H360": Hp}


def element(Fd, rec):
    ev = lambda e: pari.subst(pari(e), "x", Fd["theta"])
    return ev(rec["alpha"]) + ev(rec["beta"]) * Fd["s"]


XV = pari.varhigher("X")                                              # above the field variable y


def jinv(g):
    a2, a1, a0 = (pari.polcoef(g, k, XV) for k in (2, 1, 0))
    A = a1 - a2 ** 2 / 3
    B = a0 - a1 * a2 / 3 + 2 * a2 ** 3 / 27
    return 1728 * 4 * A ** 3 / (4 * A ** 3 + 27 * B ** 2)


def build(name):
    Fd = field()
    rec = json.load(open(HERE / f"recognised_{name}.json"))
    require(rec["F"] == str(Fd["F"]), "field polynomial")
    val = {k: element(Fd, r) for k, r in rec["values"].items()}
    x = pari("x")
    X = XV
    j1 = Fd["j1"]
    A = 3 * j1 * (1728 - j1); B = 2 * j1 * (1728 - j1) ** 2
    Q = x ** 3 + x ** 2 + val["q1"] * x + val["q0"]
    require(val["kappa"] != 0 and val["q0"] != 0, "f1 = v + kappa x / Q has degree three")
    L = val["v"] * Q + val["kappa"] * x
    H1 = L ** 3 + A * L * Q ** 2 + B * Q ** 3
    require(pari.poldegree(H1, "x") == 9, "degree of Q^3 g1(f1)")
    S1 = pari.gcd(H1, pari.deriv(H1, "x"))
    S1 = S1 / pari.polcoef(S1, pari.poldegree(S1, "x"), "x")
    require(pari.poldegree(S1, "x") == 3, "S1 is a cubic")
    Rq, rem = pari.divrem(H1, S1 ** 2, "x")
    require(rem == 0, "S1^2 divides Q^3 g1(f1)")
    c1 = pari.polcoef(Rq, 3, "x")
    R = Rq / c1
    sextic = Q * R
    require(pari.poldegree(pari.gcd(sextic, pari.deriv(sextic, "x")), "x") == 0, "Q R is squarefree")
    # second map
    require(val["p1"] == pari(2) / 3, "normalised linear coefficient of P2")
    P2 = x ** 2 + val["p1"] * x + val["p0"]
    require(pari.polresultant(P2, R, "x") != 0, "f2 = P2 / R has degree three")
    resQR = pari.polresultant(Q, R, "x")
    require(resQR != 0, "Q and R coprime")
    g2 = pari.polresultant(Q, X * R - P2, "x") / resQR
    require(pari.poldegree(g2, X) == 3 and pari.polcoef(g2, 3, X) == 1, "g2 monic cubic")
    H2 = sum(pari.polcoef(g2, k, X) * P2 ** k * R ** (3 - k) for k in range(4))
    S2 = pari.gcd(H2, pari.deriv(H2, "x"))
    S2 = S2 / pari.polcoef(S2, pari.poldegree(S2, "x"), "x")
    require(pari.poldegree(S2, "x") == 3, "S2 is a cubic")
    Qq, rem2 = pari.divrem(H2, S2 ** 2, "x")
    require(rem2 == 0, "S2^2 divides R^3 g2(f2)")
    c2 = pari.polcoef(Qq, 3, "x")
    require(Qq == c2 * Q, "the singleton fibre of phi2 is Q")
    j2 = jinv(g2)
    require(pari.subst(Fd["H360"], "x", j2) == 0, "j(E2) is a root of H_{-360}")
    require(j2 != j1, "distinct quotient j-invariants")
    tostr = lambda e: str(pari.lift(e)).replace("y", "theta")
    return {
        "field_polynomial_w": str(Fd["Pabs"]).replace("y", "w"),
        "theta_generator_note": "theta is a root of field_polynomial_w; H = Q(theta) is the ring class field of Z[3 sqrt(-10)]",
        "sextic_constant_first": [tostr(pari.polcoef(sextic, k, "x")) for k in range(7)],
        "Q": [tostr(pari.polcoef(Q, k, "x")) for k in range(4)],
        "R": [tostr(pari.polcoef(R, k, "x")) for k in range(4)],
        "phi1": {"X": "v + kappa*x/Q(x)", "v": tostr(val["v"]), "kappa": tostr(val["kappa"]),
                 "j": tostr(j1), "E": "Y^2 = X^3 + 3 j (1728 - j) X + 2 j (1728 - j)^2", "c": tostr(c1),
                 "S": [tostr(pari.polcoef(S1, k, "x")) for k in range(4)]},
        "phi2": {"X": "P2(x)/R(x)", "P2": [tostr(pari.polcoef(P2, k, "x")) for k in range(3)],
                 "E": [tostr(pari.polcoef(g2, k, X)) for k in range(4)], "j": tostr(j2), "c": tostr(c2),
                 "S": [tostr(pari.polcoef(S2, k, "x")) for k in range(4)]},
    }


MODULI_POLYNOMIAL = "x^8 - 5*x^6 + 12*x^4 - 5*x^2 + 1"


def verify_moduli_field(model):
    """Find the field generated by the absolute invariants in the J2 != 0 chart."""
    modulus = pari(model["field_polynomial_w"].replace("w", "y"))
    coefficients = [pari.Mod(pari(c.replace("theta", "y")), modulus)
                    for c in model["sextic_constant_first"]]
    invariants = absolute_invariants(pari, coefficients)
    polynomial = pari.minpoly(invariants[0], "x")
    require(polynomial.poldegree() == 8, "Expected a degree-eight first absolute invariant.")

    def vector(value):
        return pari.Col([pari.polcoef(pari.lift(value), j, "y") for j in range(16)])

    basis = pari.Mat([vector(invariants[0] ** j) for j in range(8)])
    require(all(len(basis.matinverseimage(vector(value))) == 8 for value in invariants[1:]),
            "The remaining absolute invariants do not lie in the octic subfield.")
    require(pari.polredabs(polynomial) == pari(MODULI_POLYNOMIAL), "Unexpected octic field of moduli.")


def embedding_form(model, embedding):
    # Each precision is accepted only after the enclosures certify the exact matrix.
    for root_bits, period_bits in ((2200, 400), (6000, 1200)):
        try:
            f = normalised_quintic(model["sextic_constant_first"], model["field_polynomial_w"],
                                    embedding, root_bits)
            data = period_matrix_acb(f, period_bits)
            _, order = certify_cm_order(data, 7560, discriminant=-10)
            return humbert_from_order(data, order), len(order)
        except (ArithmeticError, AssertionError, ValueError):
            if root_bits == 6000:
                raise


def certify(name, model, embeddings, target, other):
    out = {}
    for e in embeddings:
        q, rank = embedding_form(model, e)
        iso = bool(isometric(q, target))
        require(iso and not isometric(q, other), f"{name} embedding {e}: Humbert form {q} is not exactly the target")
        out[str(e)] = {"humbert_form": q, "isometric_to_target": iso, "rank_End": rank}
        print(name, "embedding", e, "Humbert form", q, "isometric to target", iso, flush=True)
    return out


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("names", nargs="*", metavar="MODEL")
    parser.add_argument("--embeddings", default=",".join(map(str, range(16))))
    parser.add_argument("--verify", action="store_true", help="compare with the archived receipts without writing")
    args = parser.parse_args()
    names = args.names or ["A", "B"]
    embs = [int(t) for t in args.embeddings.split(",")]
    if any(n not in ("A", "B") for n in names) or any(not 0 <= e < 16 for e in embs):
        parser.error("models must be A or B and embeddings must be in 0–15")
    census_path = HERE.parent / "opposite-secant-spectra" / "delay" / "census.json"
    census = json.loads(census_path.read_text())["smallest_square_isospectral"]
    targets = dict(zip(("A", "B"), census["child_grams"], strict=True))
    models_path = HERE / "models_1440.json"
    certificate_path = HERE / "certificate_1440.json"
    models = json.loads(models_path.read_text()) if models_path.exists() else {}
    certificates = json.loads(certificate_path.read_text()) if certificate_path.exists() else {}
    for name in names:
        model = build(name)
        verify_moduli_field(model)
        if args.verify:
            require(model == models[name], f"{name}: exact model differs from its receipt.")
        target = pari.matrix(3, 3, [entry for row in targets[name] for entry in row])
        counts = target.qfrep(9)
        require(counts[0] == counts[3] == 0 and counts[8] == 1,
                f"{name}: unexpected degree-two or degree-three count.")
        print(f"{name}: exact covers and octic moduli field verified", flush=True)
        other = targets["B" if name == "A" else "A"]
        results = certify(name, model, embs, targets[name], other)
        if args.verify:
            archived = certificates[name]
            require(archived["target"] == targets[name], f"{name}: archived target differs from the census.")
            for e, result in results.items():
                old = archived["embeddings"][e]
                require(old["isometric_to_target"] and result["rank_End"] == old["rank_End"]
                        and isometric(result["humbert_form"], old["humbert_form"]),
                        f"{name} embedding {e}: differs from the archived certificate.")
        else:
            models[name] = model
            certificate = certificates.setdefault(name, {"target": targets[name], "embeddings": {}})
            certificate["embeddings"].update(results)
            models_path.write_text(json.dumps(models, indent=1))
            certificate_path.write_text(json.dumps(certificates, indent=1))


if __name__ == "__main__":
    main()
