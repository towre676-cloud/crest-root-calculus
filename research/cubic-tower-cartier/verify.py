"""Reconstruct every geometric simplicity witness and the genus-55 congruence."""
import argparse
from collections import Counter, defaultdict
from itertools import combinations
import json
from pathlib import Path
import sys

import numpy as np

from base_change import powered_congruence
from cartier import Cartier, numerators
from geometric_dimensions import subdegrees

HERE = Path(__file__).resolve().parent
TOWER = HERE.parent / "cubic-tower-factors"
sys.path.insert(0, str(TOWER))
from tower import check_stage


def check_witnesses():
    rows = json.loads((HERE / "triple-witnesses.json").read_text())
    expected = {",".join(map(str, support)) + "/" + sign
                for support in combinations(range(10), 3) for sign in "+-"}
    if set(rows) != expected:
        raise ArithmeticError("Incomplete signed support coverage")
    groups = defaultdict(lambda: defaultdict(dict))
    counts = Counter()
    for key, witnesses in rows.items():
        support, sign = key.split("/")
        support = tuple(map(int, support.split(",")))
        possible = set(range(1, 7))
        for witness in witnesses:
            prime = witness["p"]
            coeffs = witness["unit_polynomial"]
            options = subdegrees(prime, tuple(coeffs))
            if options is None or list(options) != witness["permitted_half_dimensions"]:
                raise ArithmeticError("Dimension restrictions failed on " + key)
            possible.intersection_update(options)
            groups[prime][support][sign] = coeffs
        if possible:
            raise ArithmeticError("Proper dimension survived on " + key)
        counts[len(witnesses)] += 1
    for prime, supports in sorted(groups.items()):
        engine = Cartier(prime, 3)
        for support, saved in supports.items():
            check_stage(prime, support)
            polys = numerators(engine.matrices(support), prime)
            if any(poly != saved[sign] for sign, poly in zip("+-", polys) if sign in saved):
                raise ArithmeticError(f"Cartier witness differs for {support} at {prime}")
    return {"classes": len(expected), "matrix_pairs": sum(map(len, groups.values())),
            "largest_prime": max(groups), "witness_counts": {str(n): counts[n] for n in sorted(counts)}}


def crt3p(a, b, p):
    """Return the canonical coefficient satisfying x=a mod 3 and x=b mod p."""
    return b + p * (((a-b)*pow(p, -1, 3)) % 3)


def initial_coefficients(traces):
    result = [1]
    for k in range(1, len(traces)+1):
        n = -sum(result[k-i]*traces[i-1] for i in range(1, k+1))
        c, remainder = divmod(n, k)
        if remainder:
            raise ArithmeticError("Failed exact Newton identity")
        result.append(c)
    return result


def genus55():
    p = 11
    series = json.loads((HERE / "extension-counts.json").read_text())
    prev = json.loads((HERE.parent / "cubic-tower-higher-supports/mod-three.json").read_text())["11"]
    for d, saved in enumerate(prev["direct_counts"], 1):
        row = series[d-1]
        if row != [d, p**d, *saved["ordinary_and_twisted_counts"], *saved["signed_triple_traces"]]:
            raise ArithmeticError("Past extension counts disagree")
    if any(r[:2] != [n, p**n] for n,r in enumerate(series, 1)):
        raise ArithmeticError("Extension orders disagree")
    # All degrees above the dimension in a Weil numerator vanish modulo p.
    base = json.loads((TOWER / "receipt.json").read_text())["stage_two"]["complete_factors_low_to_high"]["base"]
    full = np.array([v % p for v in base], dtype=np.int64)
    triple = {}
    for k in range(1,4):
        engine = Cartier(p,k)
        for support in combinations(range(3),k):
            check_stage(p,support)
            polys = numerators(engine.matrices(support),p)
            for sign,poly in zip("+-",polys):
                full = np.convolve(full,poly) % p
                if k==3:
                    triple[sign] = poly + [0]*12
    full = full.tolist()
    full += [0]*(111-len(full))
    items = {"whole":(full, prev["full_numerator_low_to_high_mod_3"],
                      [p**d+1-r[2] for d,r in enumerate(series,1)]),
             "+":(triple["+"], prev["signed_triples_low_to_high_mod_3"]["0,1,2/+"], [r[4] for r in series]),
             "-":(triple["-"], prev["signed_triples_low_to_high_mod_3"]["0,1,2/-"], [r[5] for r in series])}
    result = {}
    for label,(modp,mod3,traces) in items.items():
        mod33 = [crt3p(a,b,p) for a,b in zip(mod3,modp,strict=True)]
        first = initial_coefficients(traces)
        if [x % 33 for x in first] != mod33[:len(first)]:
            raise ArithmeticError("Exact coefficients disagree with the two independent congruences")
        lift = powered_congruence(mod33,33,33)
        if lift["coefficient_modulus"] != 1089:
            raise ArithmeticError("Wrong promised base-change precision")
        result[label] = {"coefficients_0_through_8":first,"full_numerator_mod_33":mod33,
                         "degree_33_base_change":lift}
    return result


def small_controls():
    # Synthetic roots give exact extension polynomials independently of Newton.
    for p,m,extension,roots in [(3,3,9,[1,4,7]),(2,8,4,[2,3,18]),(11,33,33,[2,5,7,8])]:
        h=[1]; exact=[1]
        for a in roots:
            h=[x-y*a for x,y in zip(h+[0],[0]+h)]
            exact=[x-y*a**extension for x,y in zip(exact+[0],[0]+exact)]
        powered=powered_congruence([c%m for c in h],m,extension)
        if powered["numerator"] != [v%powered["coefficient_modulus"] for v in exact]:
            raise ArithmeticError("Independent base-change example failed")
    # Independent finite-fibre traces and affine residues use different bases.
    sys.path.insert(0, str(HERE.parent / "cubic-tower-geometric-separation"))
    from affine_frobenius import signed_numerators
    triple = Cartier(41, 3)
    for support in combinations(range(10), 3):
        check_stage(41, support)
        if numerators(triple.matrices(support), 41) != signed_numerators(41, support):
            raise ArithmeticError("The independent triple matrices disagree")
    # Every coefficient on the low supports was already determined integrally.
    source = json.loads((TOWER / "receipt.json").read_text())["stage_ten"]["complete_factors_low_to_high"]
    for k in (1,2):
        engine=Cartier(41,k)
        for s in combinations(range(10),k):
            check_stage(41,s)
            polys=numerators(engine.matrices(s),41)
            for sign,poly in zip("+-",polys):
                expected=source[",".join(map(str,s))+"/"+sign]
                if poly+[0]*(len(poly)-1) != [v%41 for v in expected]:
                    raise ArithmeticError("Integral low-support control failed")


if __name__=="__main__":
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--compute",action="store_true")
    args=parser.parse_args()
    small_controls()
    report={"simplicity":check_witnesses(),"genus55_at_11":genus55()}
    if args.compute:
        print(json.dumps(report, indent=2))
    elif report!=json.loads((HERE/"receipt.json").read_text()):
        raise ArithmeticError("The fresh computation differs from the receipt")
    else:
        print("Verified absolute simplicity of all 240 triple factors and genus-55 coefficients modulo 33.")
