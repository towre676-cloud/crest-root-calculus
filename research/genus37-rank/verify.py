"""Verify lower bounds from rational divisors and upper bounds from two-descent."""
import argparse
from itertools import combinations
import json
import os
from pathlib import Path
import shutil

from torsion import verify as torsion_groups
from rank import verify_sector
from upper import verify as upper_bounds

HERE = Path(__file__).resolve().parent


def independent_reductions(binary, records):
    from fractions import Fraction
    from upper import gp_value
    import subprocess
    from rank import polynomial

    code = [(HERE / "reductions.gp").read_text()]
    for key, sector in records.items():
        coeffs = [int(c) for c in polynomial(tuple(map(int, key))).all_coeffs()]
        divs = [[] if row.get("infinity") else
                [[str(Fraction(z)) if Fraction(z).denominator != 1 else int(Fraction(z))
                  for z in row[k]] for k in ("u", "v")] for row in sector["divisors"]]
        args = [divs, sector["ell"], sector["controls"], sector["witnesses"]]
        code.append("print(check_reductions(Pol(" + gp_value(coeffs) + ")," +
                    ",".join(map(gp_value, args)) + "));")
    proc = subprocess.run([binary, "-fq"], input="\n".join(code), text=True, capture_output=True, check=True)
    if proc.stderr or proc.stdout.splitlines() != ["1"] * len(records):
        raise ArithmeticError("Independent PARI divisor check failed:\n" + proc.stderr[:1800])


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp", nargs="?", const="gp", default=os.environ.get("GP") or shutil.which("gp"),
                        help="Regenerate all local images, number fields and global descent matrices with PARI/GP")
    args = parser.parse_args()
    records = json.loads((HERE / "certificate.json").read_text())["sectors"]
    if args.gp:
        independent_reductions(args.gp, records)
    upper = upper_bounds(args.gp)
    sectors = ["".join(c) for k in range(1, 5) for c in combinations("0123", k)]
    if set(upper) != set(sectors) or set(records) - set(sectors):
        raise ArithmeticError("Unexpected support in the character decomposition")
    result = {}
    for s in sectors:
        low = verify_sector(tuple(map(int, s)), records[s]) if s in records else 0
        if low > upper[s]:
            raise ArithmeticError("The independent rank bounds contradict each other")
        result[s] = [low, upper[s]]
    print(json.dumps({"sectors": result, "genus_37_rank": [sum(v[i] for v in result.values()) for i in (0,1)],
                      "torsion": torsion_groups(args.gp),
                      "exact_sectors": [s for s, v in result.items() if v[0] == v[1]]}, indent=2))


if __name__ == "__main__":
    main()
