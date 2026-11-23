"""Replay the signed Jacobian factor certificates of the S3 cubic tower."""
import argparse
from functools import cache
from itertools import combinations
import json
from math import prod
from pathlib import Path

from field import Field, coefficients, multiply
from simplicity import absolutely_simple, disjoint, power_bound
from tower import check_stage, trace_from_factor, trace_sums, weil_polynomial


def polynomial_product(factors):
    result = [1]
    for factor in factors:
        merged = [0]*(len(result)+len(factor)-1)
        for i, a in enumerate(result):
            for j, b in enumerate(factor):
                merged[i+j] += a*b
        result = merged
    return result


def check_by_enumeration(p, extension, indices):
    """Directly count both forms above every base point; no cubic histograms."""
    field = Field(p, extension)
    q = field.q
    digits = [coefficients(a, p, extension) for a in range(q)]
    mul = [[multiply(a, b, p, field.modulus) for b in range(q)] for a in range(q)]
    add = [[sum(((x+y) % p)*p**j for j, (x, y) in enumerate(zip(a, b)))
            for b in digits] for a in digits]
    cubes = [mul[x][mul[x][x]] for x in range(q)]
    squares = [mul[x][x] for x in range(q)]
    nonsquare = next(a for a in range(1, q) if a not in squares)
    result = []
    for s in (1, nonsquare):
        count = 0
        for z in range(q):
            rhs = add[add[mul[-4 % p][cubes[z]]][squares[z]]][add[z][1]]
            for w in range(q):
                if mul[s][squares[w]] != rhs:
                    continue
                fibre = 1
                for i in indices:
                    a = mul[3][add[z][i % p]]
                    fibre *= sum(add[mul[s][cubes[t]]][mul[a][t]] == w for t in range(q))
                count += fibre
        infinity = []
        for i in indices:
            radicand = ((-12*i-1)*pow(12, -1, p)) % p
            infinity.append(1+sum(mul[s][b2] == radicand for b2 in squares))
        result.append(count+prod(infinity))
    return result


def run_stage(p, indices, orders):
    traces = {}
    base = []
    counts = []
    for extension, order in enumerate(orders, 1):
        tr, data, full = trace_sums(p, indices, extension, order)
        base.append(tr)
        counts.append(full)
        for support, pair in data.items():
            seq = traces.setdefault(support, [[], []])
            for sign in (0, 1):
                seq[sign].append(pair[sign])
    base_factor = weil_polynomial(base, p, 1)
    complete = {"base": base_factor}
    for support, pair in traces.items():
        dimension = len(support)*2**(len(support)-1)
        for sign, seq in zip(("+", "-"), pair):
            if len(seq) < dimension:
                continue
            factor = weil_polynomial(seq, p, dimension)
            if any(trace_from_factor(factor, j) != v for j, v in enumerate(seq, 1)):
                raise ArithmeticError("A sign-factor Frobenius power failed replay")
            complete[",".join(str(i) for i in support)+"/"+sign] = factor
    if any(trace_from_factor(base_factor, j) != v for j, v in enumerate(base, 1)):
        raise ArithmeticError("Elliptic base failed replay")
    genus = 1 + 2*len(indices)*3**(len(indices)-1)
    degree = sum(len(f)-1 for f in complete.values())
    result = {
        "indices": list(indices), "p": p, "genus": genus,
        "ordinary_and_involution_twisted_counts": counts,
        "certified_factor_degree": degree, "total_numerator_degree": 2*genus,
        "complete_factors_low_to_high": complete,
    }
    if degree == 2*genus:
        numerator = polynomial_product(complete.values())
        if any(p**j+1-trace_from_factor(numerator, j) != c[0] for j, c in enumerate(counts, 1)):
            raise ArithmeticError("Complete numerator failed point-count replay")
        result["full_numerator_low_to_high"] = numerator
    return result


def classify_factors(stage):
    """Certify absolute simplicity and pairwise nonisogeny across all computed factors."""
    stage_prime = stage["p"]
    factors = stage["complete_factors_low_to_high"]
    alternatives = [p for p in (5, 7, 11, 13, 17, 19) if p != stage_prime]

    @cache
    def smaller_stage(prime, indices):
        try:
            check_stage(prime, indices)
        except ValueError:
            return {}
        depth = len(indices)*2**(len(indices)-1)
        return run_stage(prime, indices, [len(indices)]*depth)["complete_factors_low_to_high"]

    def auxiliary_factor(prime, label):
        if label == "base":
            # Use any one good cover to compute the common elliptic base factor.
            for i in stage["indices"]:
                factor = smaller_stage(prime, (i,)).get("base")
                if factor is not None:
                    return factor
            return None
        indices = tuple(map(int, label.split("/")[0].split(",")))
        return smaller_stage(prime, indices).get(label)

    exceptions_simple = {}
    for name, f in factors.items():
        if len(f) == 3 or absolutely_simple(f):
            continue
        for p in alternatives:
            other = auxiliary_factor(p, name)
            if other is not None and absolutely_simple(other):
                exceptions_simple[name] = p
                break
        else:
            raise ArithmeticError(f"No absolute-simplicity certificate for {name}")

    exceptions_disjoint = {}
    pairs = 0
    for a, b in combinations(factors, 2):
        if len(factors[a]) != len(factors[b]):
            continue
        pairs += 1
        if disjoint(factors[a], factors[b]):
            continue
        for p in alternatives:
            f, g = auxiliary_factor(p, a), auxiliary_factor(p, b)
            if f is not None and g is not None and disjoint(f, g):
                exceptions_disjoint[f"{a} | {b}"] = p
                break
        else:
            raise ArithmeticError(f"No geometric nonisogeny certificate for {a} and {b}")
    return {
        "rational_frobenius_power_multiple": power_bound(8),
        "auxiliary_coefficient_prime": 1000003,
        "absolutely_simple_factors": len(factors),
        "equal_dimension_pairs_separated": pairs,
        "simplicity_witness_prime_when_different": exceptions_simple,
        "separation_witness_prime_when_different": exceptions_disjoint,
    }


def benchmark():
    first = run_stage(11, [0, 1], [2]*5)
    tenth = run_stage(41, list(range(10)), [10, 10, 2, 2])
    if first["certified_factor_degree"] != 26 or tenth["certified_factor_degree"] != 762:
        raise ArithmeticError("Factor coverage does not match the dimension theorem")
    for k in (1, 2):
        if check_by_enumeration(11, k, [0, 1]) != first["ordinary_and_involution_twisted_counts"][k-1]:
            raise ArithmeticError("Small-extension direct enumeration disagrees")
    prior = {}
    for p, expected in ((41, 81), (53, 2268), (59, 2952), (71, 243)):
        count = check_by_enumeration(p, 1, list(range(10)))
        if count[0] != expected:
            raise ArithmeticError("Stage-ten historical count disagrees")
        prior[str(p)] = count
    if prior["41"] != tenth["ordinary_and_involution_twisted_counts"][0]:
        raise ArithmeticError("Two stage-ten counters disagree")
    try:
        check_stage(5, [0, 1])
    except ValueError:
        pass
    else:
        raise ArithmeticError("The bad reduction at 5 passed the stage check")
    first["geometric_factor_checks"] = classify_factors(first)
    tenth["geometric_factor_checks"] = classify_factors(tenth)
    return {"stage_two": first, "stage_ten": tenth, "stage_ten_prime_field_checks": prior}


def main():
    cli = argparse.ArgumentParser(description=__doc__)
    cli.add_argument("--compute", action="store_true", help="print the recomputed certificate")
    args = cli.parse_args()
    result = benchmark()
    if args.compute:
        print(json.dumps(result, indent=2))
    else:
        saved = json.loads(Path(__file__).with_name("receipt.json").read_text(encoding="utf-8"))
        if result != saved:
            raise ArithmeticError("Saved certificate differs")
        print("Verified the stage-two numerator, 111 simple nonisogenous stage-ten factors, and point counts.")


if __name__ == "__main__":
    main()
