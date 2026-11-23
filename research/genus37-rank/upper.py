"""Replay the 15 two-descent upper bounds and, optionally, regenerate them with PARI/GP."""
import json
import os
from pathlib import Path
import re
import shutil
import subprocess

HERE = Path(__file__).resolve().parent


def gp_value(value):
    if isinstance(value, list):
        return "[" + ",".join(map(gp_value, value)) + "]"
    if isinstance(value, int):
        return str(value)
    if isinstance(value, str) and re.fullmatch(r"-?\d+/[1-9]\d*", value):
        return value
    if value in ("branch", "quadratic", "extension"):
        return json.dumps(value)
    raise ValueError("Invalid local divisor witness")


def recompute(binary):
    inputs = json.loads((HERE / "local-witnesses.json").read_text())
    code = ['default(parisizemax,1000000000);', (HERE / "selmer.gp").read_text()]
    for key, case in inputs.items():
        places = [0, *case["primes"]]
        if set(map(int, case["places"])) != set(places):
            raise ValueError("Missing local witness")
        fields = [int(i) + 1 for i in key]
        witnesses = [case["places"][str(p)]["generators"] for p in places]
        code.append("sector(" + ",".join(map(gp_value, [fields, case["primes"], witnesses])) + ");")
    proc = subprocess.run([binary, "-fq"], input="\n".join(code), text=True,
                          capture_output=True, check=True)
    diagnostics = [line for line in proc.stderr.splitlines() if "Warning: new maximum stack" not in line]
    if diagnostics:
        raise RuntimeError("\n".join(diagnostics))
    records = {}
    for line in proc.stdout.splitlines():
        fields, dim, diagonal, torsion, upper, bits = json.loads(line)
        key = "".join(str(i - 1) for i in fields)
        records[key] = dict(dimension=dim, diagonal=diagonal, two_torsion=torsion,
                            rank_upper_bound=upper, rows=[hex(b) for b in bits])
    if set(records) != set(inputs):
        raise ArithmeticError("Incomplete PARI computation")
    return records


def matrix_rank(rows, dimension):
    pivots = {}
    for row in rows:
        r = int(row, 16)
        if r < 0 or r.bit_length() > dimension:
            raise ValueError("Invalid binary restriction row")
        while r:
            i = r.bit_length() - 1
            if i not in pivots:
                pivots[i] = r
                break
            r ^= pivots[i]
    return len(pivots)


def bounds(records):
    result = {}
    for key, data in records.items():
        dimension = data["dimension"]
        rank = matrix_rank(data["rows"], dimension)
        upper = dimension - rank - data["diagonal"] - data["two_torsion"]
        if upper != data["rank_upper_bound"] or upper < 0:
            raise ArithmeticError("Invalid descent matrix or rank bound")
        result[key] = upper
    return result


def verify(binary=None):
    stored = json.loads((HERE / "descent-matrices.json").read_text())
    result = bounds(stored)
    if binary:
        rebuilt = bounds(recompute(binary))
        if rebuilt != result:
            raise ArithmeticError("The independently regenerated upper bounds differ")
    return result


if __name__ == "__main__":
    import argparse
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--gp", nargs="?", const="gp", default=os.environ.get("GP") or shutil.which("gp"))
    args = parser.parse_args()
    print(json.dumps(verify(args.gp), indent=2))
