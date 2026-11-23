"""Independent characteristic-zero audit of the standard atlas root-ratio witnesses.

Run: python original/audit_powers.py [receipt.json]. Requires SymPy.
This checks primed Frobenius polynomials recorded in a receipt, independently of
the finite-field routines. Point counts, smoothness and specialization are
separate parts of the atlas certificate.
"""
import argparse
import json
from math import gcd
from pathlib import Path
import sympy as sp


def root_power_sums(coefficients, limit):
    """Newton sums s_k = sum(r^k) of all roots r, for a descending monic list."""
    degree = len(coefficients) - 1
    if degree != 6 or coefficients[0] != 1:
        raise ValueError("Expected a monic sextic Frobenius polynomial")
    sums = [degree]
    for k in range(1, limit + 1):
        if k <= degree:
            previous = sum(coefficients[j] * sums[k-j] for j in range(1,k))
            sums.append(-previous-k*coefficients[k])
        else:
            sums.append(-sum(coefficients[j]*sums[k-j] for j in range(1,degree+1)))
    return sums


def powered_polynomial(power_sums, m, x):
    """Form product (x-r^m) using the power sums s_m through s_6m."""
    coefficients = [1]
    for k in range(1, 7):
        numerator = -sum(power_sums[m*j]*coefficients[k-j] for j in range(1,k+1))
        coefficient, remainder = divmod(numerator, k)
        if remainder:
            raise ArithmeticError("Newton identity produced a nonintegral coefficient")
        coefficients.append(coefficient)
    return sp.Poly.from_list(coefficients, x, domain=sp.ZZ)


def audit(receipt):
    # A pair of different roots has at most 6*5=30 conjugates. If their
    # quotient is primitive of order m, phi(m)<=30. Since phi(m)>=sqrt(m/2),
    # all possibilities occur by 1800. Enumerate by gcd directly.
    orders = [m for m in range(2,1801) if sum(gcd(j,m)==1 for j in range(1,m+1))<=30]
    x = sp.Symbol("x")
    output = []
    for row in receipt["curves"]:
        if row["status"] != "no_elliptic_map_over_Qbar":
            continue
        witnesses = [r for r in row["reductions"] if r.get("simplicity_certificate")]
        if not witnesses:
            raise ValueError("Missing root-ratio witness in "+row["id"])
        for witness in witnesses:
            coefficients = list(reversed(witness["frobenius_ascending"]))
            poly = sp.Poly.from_list(coefficients, x, domain=sp.ZZ)
            if not poly.is_irreducible:
                raise ValueError("Reducible source polynomial for "+row["id"])
            sums = root_power_sums(coefficients, 6*max(orders))
            for m in orders:
                if not powered_polynomial(sums, m, x).is_sqf:
                    raise ValueError(f"Root collision at exponent {m} for {row['id']}")
            output.append({"id":row["id"], "good_prime":witness["p"], "irreducible":True,
                           "squarefree_power_polynomials":len(orders)})
    return {"tested_orders":orders, "rows":output,
            "total_squarefree_checks":sum(r["squarefree_power_polynomials"] for r in output)}


if __name__ == "__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("receipt",type=Path,nargs="?",default=Path(__file__).resolve().parents[1]/"receipt.json")
    args=parser.parse_args()
    try:
        print(json.dumps(audit(json.loads(args.receipt.read_text(encoding="utf-8"))),indent=2))
    except (KeyError, ValueError, ArithmeticError, OSError) as exc:
        parser.exit(2, "Audit error: "+str(exc)+"\n")
