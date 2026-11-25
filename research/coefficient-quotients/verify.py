"""Replay quotient certificates, completeness controls and the frozen experiment."""

from collections import Counter
from pathlib import Path
import statistics
import sys

import sympy as sp

local_import_path = sys.path[:]
from panel import generate
sys.path[:] = local_import_path
from benchmark import (HERE, REPETITIONS, TIMEOUTS, assert_identities, digest,
                       encode, expected_and_schedule, map_list, read, write)
from discover import (X, coefficients, constant_degrees, discover, nonconstant_degrees, normalization, polynomial,
                      variable_incidence)
from finite import rational_points, verify as verify_finite
from freeze import source_paths


def check(condition, reason):
    if not condition:
        raise ArithmeticError(reason)


def independent_constant(f):
    """Use leading cancellation and direct powers, separate from both timing paths."""
    result = []
    for d, m in constant_degrees(f):
        work = sp.Poly(X**d, X, domain=sp.QQ)
        leading = f.LC()
        for k in range(1, d):
            i = m * d - k
            missing = (f.nth(i) / leading - (work**m).nth(i)) / m
            work += sp.Poly(missing * X**(d-k), X, domain=sp.QQ)
        residual = f
        out = {}
        for j in range(m, -1, -1):
            out[j] = residual.nth(d*j)
            if out[j]:
                residual -= (work**j).mul_ground(out[j])
        if not residual:
            result.append({"degree": d, "target_degree": m, "inner": coefficients(work),
                           "outer": coefficients(sp.Poly.from_dict(out, X, domain=sp.QQ))})
    return result


def verify_variable_row(source, row):
    d, m, e = (row[key] for key in ("degree", "target_degree", "vertical_degree"))
    variables, equations, inner, vertical = variable_incidence(source, d, m, e)
    cert = row["certificate"]
    if variables:
        check(cert["variables"] == list(map(str, variables))
              and cert["equations"] == list(map(str, equations)), "Certificate changed its incidence.")
        verify_finite(cert)
        points, count = cert["rational_points"], cert["geometric_points"]
    else:
        check(cert is None, "Unexpected finite-algebra certificate in a zero-variable stratum.")
        points = [] if equations else [[]]
        count = len(points)
    check(row["geometric_points"] == count, "Geometric support count changed.")
    certified = []
    for point in points:
        assignment = dict(zip(variables, map(sp.Rational, point), strict=True))
        u = sp.Poly(inner.as_expr().subs(assignment), X, domain=sp.QQ)
        v = sp.Poly(vertical.as_expr().subs(assignment), X, domain=sp.QQ)
        certified.append((coefficients(u), coefficients(v)))
    observed = []
    for mapping in row["maps"]:
        u, v, g = map(polynomial, (mapping["inner"], mapping["vertical"], mapping["outer"]))
        check(u.degree() == d and v.degree() == e and g.degree() == m and u.LC() == v.LC() == 1
              and u.nth(0) == 0 and g.gcd(g.diff()).degree() == 0 and u.diff().rem(v).is_zero
              and g.compose(u) == source * v**2, "Invalid exported nonconstant map.")
        observed.append((mapping["inner"], mapping["vertical"]))
    check(sorted(observed) == sorted(certified), "Rational map list is incomplete.")


def raw_cubic(source, e):
    """Eliminate raw coefficients by exact leading cancellation and Gröbner reduction."""
    horizontal = sp.symbols("a1:3")
    vertical = sp.symbols(f"b0:{e}")
    outer = sp.symbols("c0:3")
    u = X**3 + sum(a * X**i for i, a in enumerate(horizontal, 1))
    v = X**e + sum(b * X**i for i, b in enumerate(vertical))
    g = source.LC() * X**3 + sum(c * X**i for i, c in enumerate(outer))
    identity = sp.Poly(g.subs(X, u) - source.as_expr() * v**2, X)
    substitutions = []
    # In each chosen coefficient the pivot occurs linearly with a nonzero
    # constant coefficient. Substitution therefore gives exact elimination.
    for power, pivot in zip((8, 7, 6, 3, 0), (*reversed(horizontal), *reversed(outer)), strict=True):
        equation = sp.Poly(identity.nth(power), pivot)
        check(equation.degree() == 1 and equation.LC().is_Rational,
              "Raw coefficient elimination lost its nonzero constant pivot.")
        value = -equation.nth(0) / equation.LC()
        substitutions.append((pivot, value))
        identity = sp.Poly(identity.as_expr().subs(pivot, value), X)
    result = rational_points(identity.coeffs(), vertical)
    maps = []
    for point in result["rational_points"]:
        assignment = list(zip(vertical, map(sp.Rational, point), strict=True))
        maps.append({"inner": coefficients(sp.Poly(u.subs(substitutions + assignment), X)),
                     "outer": coefficients(sp.Poly(g.subs(substitutions + assignment), X)),
                     "vertical": coefficients(sp.Poly(v.subs(assignment), X))})
    return {"geometric_points": result["geometric_points"], "maps": sorted(maps, key=encode)}


def variable_output(row):
    return {"geometric_points": row["geometric_points"], "maps": sorted(row["maps"], key=encode)}


def independent_provenance(source, provenance, maps):
    if "diagonal" in provenance:
        diagonal, bonds, repeats = (provenance[key] for key in ("diagonal", "bonds", "repetitions"))
        check(all(a + bonds[i - 1] + bonds[i] == p >= 0
                  for i, (a, p) in enumerate(zip(diagonal, provenance["pinning"], strict=True))),
              "Jacobi metadata do not give the reported pinning model.")
        # Full transfer products verify the expanded polynomial without the
        # repeated-cell recurrence used when the panel was created.
        matrix = sp.eye(2)
        for i in range(len(bonds) * repeats):
            j = i % len(bonds)
            matrix = (sp.Matrix([[X - diagonal[j], -bonds[j - 1]], [bonds[j], 0]]) * matrix).applyfunc(sp.expand)
        expected = sp.Poly(sp.trace(matrix)**2 - 4 * sp.prod(bonds)**(2*repeats), X, domain=sp.QQ)
        check(expected == source, "A Jacobi-derived radicand disagrees with its independent transfer product.")
    elif "inner" in provenance:
        u, g = map(polynomial, (provenance["inner"], provenance["outer"]))
        if provenance["study"] == "nonconstant":
            v = polynomial(provenance["vertical"])
            shift = provenance["shift"]
            expected = g.compose(u).exquo(v**2).shift(shift)
            if "perturbation" not in provenance:
                u, v = u.shift(shift), v.shift(shift)
                g = g.shift(u.nth(0))
                u -= u.nth(0)
                check({"inner": coefficients(u), "vertical": coefficients(v), "outer": coefficients(g)} in maps,
                      "Missing a planted ramification map.")
        else:
            expected = g.compose(u)
            if provenance.get("affine"):
                scale, shift, unit = map(sp.Rational, provenance["affine"])
                expected = expected.compose(sp.Poly(scale*X+shift, X, domain=sp.QQ)).mul_ground(unit)
            if "inflator" in provenance:
                expected *= polynomial(provenance["inflator"])**2
        if "perturbation" in provenance:
            k, amount = provenance["perturbation"]
            expected += sp.Poly(sp.Rational(amount)*X**k, X, domain=sp.QQ)
        check(source == expected, "Expanded input differs from its source provenance.")


def compare(rows, first, second, case_names, metric):
    def field(row):
        if metric == "wall":
            return row["wall"]
        if metric == "cpu":
            return row["process_cpu"]
        return sum(value for name, value in row["phases"].items() if "gp_process" not in name)
    ratios = []
    for case in case_names:
        sample = [[field(row) for row in rows if row["case"] == case and row["method"] == method
                   and row["status"] == "completed"] for method in (first, second)]
        if all(len(times) == REPETITIONS for times in sample):
            ratios.append(statistics.median(sample[1]) / statistics.median(sample[0]))
    return {"paired_inputs": len(ratios), "first_wins": sum(r > 1 for r in ratios),
            "second_wins": sum(r < 1 for r in ratios),
            "median_ratio": round(statistics.median(ratios), 6) if ratios else None}


def audit(directory):
    names = ("freeze", "inputs", "provenance", "constant-reference", "measurements")
    freeze, inputs, sidecar, reference, measured = (read(directory/f"{name}.json") for name in names)
    check(set(freeze["sha256"]) == set(source_paths()), "The freeze omits or adds sources.")
    for relative, value in freeze["sha256"].items():
        check(digest((HERE.parent.parent / relative).read_bytes()) == value, f"Frozen source changed: {relative}")
    check(generate(freeze["panel_seed"]) == (inputs, sidecar), "The seeded held-out panel cannot be regenerated.")
    measured_protocol = measured["protocol"]
    check(measured_protocol["limits_seconds"] == TIMEOUTS and measured_protocol["repetitions"] == REPETITIONS
          and measured_protocol["shuffle_seed"] == freeze["schedule_seed"], "Protocol mismatch.")
    check(set(measured_protocol["file_sha256"]) == {f"{n}.json" for n in names[:-1]}, "Incomplete file hashes.")
    for filename, value in measured_protocol["file_sha256"].items():
        check(digest((directory/filename).read_bytes()) == value, f"The timed file changed: {filename}")
    schedule, _ = expected_and_schedule(inputs, sidecar, freeze["schedule_seed"])
    rows = measured["rows"]
    actual = [tuple(row[n] for n in ("case", "study", "method", "repetition", "stratum")) for row in rows]
    check(actual == schedule, "The measured schedule is missing, reordered or duplicated.")

    expected_variable, raw_cases = {}, set()
    stats = Counter()
    for name, entry in inputs.items():
        check(set(entry) == {"coefficients"}, "A worker input exposes provenance.")
        values = entry["coefficients"]
        p = polynomial(values)
        a, f = normalization(p)
        info = sidecar["cases"][name]
        result = discover(values, nonconstant=info["study"]=="nonconstant")
        standard = independent_constant(f)
        check(map_list(result["constant"]) == standard, "Independent leading-power search differed.")
        positive = {"square": coefficients(a), "source": coefficients(f), "maps": standard}
        assert_identities(positive, values)
        stats["constant_degree_checks"] += len(constant_degrees(f))
        if info["study"] == "constant":
            check(reference[name] == positive, "Incorrect frozen constant reference.")
            stats["constant_maps"] += len(standard)
            independent_provenance(p, info, standard)
            continue
        check([(r["degree"], r["target_degree"], r["vertical_degree"]) for r in result["nonconstant"]]
              == nonconstant_degrees(f), "An admissible stratum was skipped.")
        d, m, e = info["degree_stratum"]
        for row in result["nonconstant"]:
            verify_variable_row(f, row)
            stats["nonconstant_strata"] += 1
            stats["nonconstant_maps_all_strata"] += len(row["maps"])
            if [row["degree"], row["target_degree"], row["vertical_degree"]] == [d, m, e]:
                expected_variable[name] = variable_output(row)
                raw_key = (e, "perturbation" in info)
                if d == m == 3 and raw_key not in raw_cases:
                    check(raw_cubic(f, e) == expected_variable[name], "Raw cubic elimination disagreed.")
                    raw_cases.add(raw_key)
        independent_provenance(p, info, expected_variable[name]["maps"])
    check(len(raw_cases) == 4, "The raw cubic check lacks a positive or perturbed stratum.")

    for row in rows:
        status = row["status"]
        check(row["wall"] > 0 and row["process_cpu"] >= 0 and status in
              {"completed", "wrong_answer", "timeout", "process_error"}, "Invalid timing/status data.")
        if "phases" in row:
            check(all(v >= 0 for v in row["phases"].values()), "Negative phase time.")
        name = row["case"]
        if row["study"] == "constant" and status == "completed":
            check(row["output_sha256"] == digest(encode(reference[name])), "Wrong output hash.")
        elif row["study"] == "constant" and status == "wrong_answer":
            check(row["output"] != reference[name] and row["output_sha256"] == digest(encode(row["output"])),
                  "Misclassified wrong answer.")
        elif row["study"] == "nonconstant" and status == "completed":
            output = row["output"]
            check(row["output_sha256"] == digest(encode(output))
                  and {**output, "maps": sorted(output["maps"], key=encode)} == expected_variable[name],
                  "A timed nonconstant result differs from its verified certificate.")

    by_method = {}
    for row in rows:
        by_method.setdefault(row["method"], Counter())[row["status"]] += 1
    const = [n for n, v in sidecar["cases"].items() if v["study"] == "constant"]
    var = list(expected_variable)
    paired = {m: {metric: compare(rows, "crest", m, const, metric) for metric in ("cpu", "wall", "solver")}
              for m in ("crest-no-sieve", "corrected-dense", "singular", "gp", "stock-sympy")}
    cohorts = {}
    for cohort in sorted({info["cohort"] for info in sidecar["cases"].values()}):
        named = [n for n in const if sidecar["cases"][n]["cohort"] == cohort]
        if named:
            cohorts[cohort] = {
                "cases": len(named), "with_quotient": sum(bool(reference[n]["maps"]) for n in named),
                "quotients": sum(len(reference[n]["maps"]) for n in named),
                "solver_comparisons": {m: compare(rows, "crest", m, named, "solver")
                                       for m in ("corrected-dense", "singular", "gp")}}
    return {"attempts": len(rows), "constant_inputs": len(const), "nonconstant_inputs": len(var),
            "raw_cubic_cross_checks": len(raw_cases), "algebra": dict(stats),
            "by_method": {m: dict(count) for m, count in sorted(by_method.items())},
            "constant_comparisons": paired, "cohorts": cohorts,
            "nonconstant_comparisons": {m: compare(rows, "derivative", "without-derivative", var, m)
                                        for m in ("cpu", "wall", "solver")}}


def higher_genus_control():
    """Replay the full public search on an inflated quintic target fixture."""
    u = sp.Poly(X**3-3*X, X, domain=sp.QQ)
    v = sp.Poly(X**2-1, X, domain=sp.QQ)
    g = sp.Poly(-sp.Rational(7,5)*(X**2-4)*(X**3+X+3), X, domain=sp.QQ)
    f = g.compose(u).exquo(v**2)
    inflator = sp.Poly(X**2+X+1, X, domain=sp.QQ)
    result = discover(coefficients(f*inflator**2), nonconstant=True)
    check(result["source"] == coefficients(f) and result["square"] == coefficients(inflator),
          "Higher-genus control lost the nonsquare rational unit or normalization factor.")
    expected = {(5, 3, 2): 0, (7, 3, 5): 0, (9, 3, 8): 0, (3, 5, 2): 1}
    actual = {(row["degree"], row["target_degree"], row["vertical_degree"]): row["geometric_points"]
              for row in result["nonconstant"]}
    check(actual == expected, "The higher-genus target control lost a stratum.")
    for row in result["nonconstant"]:
        verify_variable_row(f, row)
    target = next(row for row in result["nonconstant"] if row["target_degree"] == 5)
    check(target["maps"] == [{"inner": coefficients(u), "outer": coefficients(g),
                              "vertical": coefficients(v)}], "Wrong map to the genus-two target.")


if __name__ == "__main__":
    directory = HERE / "study"
    results = audit(directory)
    higher_genus_control()
    if "--write-summary" in sys.argv:
        write(directory / "summary.json", results)
    else:
        check(results == read(directory / "summary.json"), "The summary differs from the frozen replay.")
    print(f"Verified {results['attempts']} scheduled attempts, {results['algebra']['nonconstant_strata']} held-out strata,"
          f" {results['raw_cubic_cross_checks']} raw elimination controls, and an inflated genus-two target.")
