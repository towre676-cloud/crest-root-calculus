"""Finite-order bounds from good reduction on the fifteen character Jacobians."""
from functools import reduce
from math import gcd, prod
from pathlib import Path
import json
import subprocess
import sympy as sp

from rank import polynomial

HERE = Path(__file__).resolve().parent


def torsion_bound(rows):
    first, order = rows[0]
    possible = set(sp.factorint(order)) | {first}
    bound = 1
    for ell in possible:
        others = [n for p, n in rows if p != ell]
        if not others:
            raise ValueError("Missing prime-to-residue-characteristic comparison")
        common = reduce(gcd, others)
        while common % ell == 0:
            bound *= ell
            common //= ell
    return bound


def verify(binary=None):
    d = json.loads((HERE / "torsion-orders.json").read_text())
    common_prime = d["common_good_prime"]
    groups = {}
    code = ["default(parisizemax,500000000);"]
    expected = []
    for key, rows in d["sector_orders"].items():
        indices = tuple(map(int, key))
        f = polynomial(indices)
        sector_rank = len(indices) - (1 if len(indices) % 2 else 2)
        if torsion_bound(rows) != 2**sector_rank:
            raise ArithmeticError("The reduction bound does not match rational branch two-torsion")
        groups[key] = sector_rank
        for p, n in dict([*rows, [common_prime, d["common_orders"][key]]]).items():
            if p == 2 or not sp.isprime(p) or f.discriminant() % p == 0 or n <= 0:
                raise ValueError("Invalid good-prime data")
            if binary:
                coeffs = list(map(int, f.all_coeffs()))
                code.append(f"print(subst(hyperellcharpoly(Mod(1,{p})*Pol({coeffs})),x,1));")
                expected.append(str(n))
    if binary:
        result = subprocess.run([binary, "-fq"], input="\n".join(code), text=True, capture_output=True, check=True)
        errors = [l for l in result.stderr.splitlines() if "Warning: " not in l]
        if errors or result.stdout.splitlines() != expected:
            raise ArithmeticError("Regenerated finite Jacobian orders differ")
    n = prod(d["common_orders"].values())
    return {"elementary_two_ranks": groups, "whole_group_exponent_divides": 8,
            "whole_group_order_divides": n & -n, "common_good_prime": common_prime}
