"""Replay spectra and optionally reconstruct the Gram matrices and saturation."""
import argparse
import json
from collections import Counter
from math import gcd, isqrt
from pathlib import Path

from isotropic import spectrum_from_lattice

HERE = Path(__file__).resolve().parent
CASES = ("C", "Cp", "T1", "T2", "0", "1", "2", "3", "4", "5", "cm_102", "cm_120")

FIX = {"identity": 4, "transposition": 2, "3-cycle": 1}


def bielliptic_law(k, t, X):
    cnt = Counter()
    Ra, Rb = isqrt(2 * k * k * X) + 1, isqrt(2 * k * X) + 1
    for a in range(0, Ra + 1):
        for b in range(-Rb, Rb + 1):
            if (a == 0 and b != 1) or gcd(a, b) != 1:
                continue
            two = FIX[t] if (a % 2 and b % 2) else 1
            kk = k if (a % k == 0) else 1
            d = 2 * (a * a + k * b * b) // (two * kk)
            if d <= X:
                cnt[d] += 1
    return dict(sorted(cnt.items()))


def theorem41_q0(X):
    """T: degrees 4m with multiplicity 2^(omega(m)+2)... via the doubled C_2 law: d = 2(a^2+b^2)
    (both odd) else 4(a^2+b^2) is the doubled spectrum of C_2, i.e. 2 x {2m} with the same counts."""
    cnt = Counter()
    R = isqrt(X) + 2
    for a in range(0, R + 1):
        for b in range(-R, R + 1):
            if (a == 0 and b != 1) or gcd(a, b) != 1:
                continue
            d = (a * a + b * b) if (a % 2 and b % 2) else 2 * (a * a + b * b)
            if 2 * d <= X:
                cnt[2 * d] += 1
    return dict(sorted(cnt.items()))


def full_replay(name, receipt):
    from cm_case import cm_case
    from run_known import known_case

    result = cm_case(name[3:], receipt) if name.startswith("cm_") else known_case(name, receipt)
    critical = ("f", "gram_prime", "provider_degrees", "gram", "rank_gram",
                "saturation_log", "NS_basis_in_generators", "NS_gram", "NS_discriminant",
                "theta", "spectrum_through_X", "X", "status")
    canonical = json.loads(json.dumps(result, default=str))
    if isinstance(receipt["gram_prime"], int):
        canonical["gram_prime"] = canonical["gram_prime"][0]
    for key in critical:
        if canonical.get(key) != receipt[key]:
            raise ArithmeticError(f"{name} full replay differs at {key}.")
    if result["rank_certificate"]["rho_upper"] != receipt["rank_certificate"]["rho_upper"]:
        raise ArithmeticError(f"{name} rank bound differs.")


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", metavar="CASE")
    parser.add_argument("--full", action="store_true", help="Recompute maps, Gram and saturation checks")
    parser.add_argument("--cm-bound", type=int, default=40)
    args = parser.parse_args()
    if any(name not in CASES for name in args.cases) or args.cm_bound < 0:
        parser.error("Use a known case and a nonnegative CM bound.")
    references = json.loads((HERE / "period_reference.json").read_text(encoding="utf-8"))
    cm_checks = []
    for name in args.cases or CASES:
        c = json.loads((HERE / f"cert_{name}.json").read_text(encoding="utf-8"))
        if c.get("status") != "certified":
            raise ArithmeticError(f"{name} receipt has no certificate.")
        spec, classes = spectrum_from_lattice(c["NS_gram"], c["theta"], c["X"])
        stored = {int(k): v for k, v in c["spectrum_through_X"].items()}
        if len({tuple(v) for d, v in classes}) != len(classes) or spec != stored:
            raise ArithmeticError(f"{name} spectrum differs from its receipt.")
        if name in ("0", "1", "2", "3", "4", "5"):
            independent = bielliptic_law(3, c["expected_type"], c["X"])
        elif name in ("T1", "T2"):
            independent = theorem41_q0(c["X"])
        elif name.startswith("cm_"):
            independent = {int(k): v for k, v in references["cm_102"].items()}
        else:
            independent = {2 if name == "C" else 4: 2}
        if independent != spec:
            raise ArithmeticError(f"{name} spectrum differs from its reference.")
        print(f"{name}: spectrum through {c['X']} matches the independent reference.", flush=True)
        if args.full:
            full_replay(name, c)
            print(f"{name}: exact maps, Gram matrix and saturation replayed.", flush=True)
        if name.startswith("cm_") and args.cm_bound > c["X"]:
            extended, all_classes = spectrum_from_lattice(c["NS_gram"], c["theta"], args.cm_bound)
            if {k: v for k, v in extended.items() if k <= c["X"]} != spec:
                raise ArithmeticError(f"{name} extended enumeration changes low degrees.")
            cm_checks.append(extended)
            print(f"{name}: {len(all_classes)} classes in {len(extended)} degrees through {args.cm_bound}.", flush=True)
    if len(cm_checks) == 2 and cm_checks[0] != cm_checks[1]:
        raise ArithmeticError("The extended CM spectra disagree.")
    print("All requested checks passed.")


if __name__ == "__main__":
    main()
