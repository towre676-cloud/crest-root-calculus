"""Enumerate level-three opposite secants with bounded Gram determinant."""

import argparse
from collections import Counter, defaultdict
from fractions import Fraction
from itertools import combinations, product
import json
from math import gcd, lcm
from pathlib import Path
import subprocess
import sys
from tempfile import TemporaryDirectory

HERE = Path(__file__).resolve().parent
sys.path.insert(0, str(HERE.parent.parent / "forced-square-kernel"))
from star import (act, cross, gram_matrix, matrix_entries, normal_partition, orbits,
                  pari, residue_field, sublattice)
from neighbours import neighbours


def local_obstruction(parent, primes):
    """Test representation of 1 over Z_p at all primes of the child determinant."""
    det, a, b, c, u, v, w = parent
    if not any((a*x*x+b*y*y+c*z*z+2*(u*x*y+v*x*z+w*y*z)) % 8 == 1
               for x, y, z in product(range(4), repeat=3)):
        return 2
    for p in primes:
        if p == 2:
            continue
        minors = (a*b-u*u, a*c-v*v, b*c-w*w, a*w-u*v, b*v-u*w, c*u-v*w)
        if any(m % p for m in minors):
            continue
        value = next((d % p for d in (a, b, c) if d % p), 0)
        if not value or pow(value, (p-1)//2, p) != 1:
            return p
    # At 3 the two planes contain the same line of norm one.
    return 0


def secant_pairs(parent, g, field):
    vertices, _ = normal_partition(parent, field)
    edges = list(combinations(sorted(vertices), 2))
    actions = [{x: act(matrix_entries(k), x, 3) for x in vertices} for k in pari.qfauto(g)[1]]
    labels, _ = orbits(edges, [
        lambda edge, table=table: tuple(sorted((table[edge[0]], table[edge[1]])))
        for table in actions
    ])
    representatives = {}
    for edge in edges:
        other = tuple(sorted(set(vertices)-set(edge)))
        i, j = labels[edge], labels[other]
        if i != j:
            representatives.setdefault(tuple(sorted((i, j))), (edge, other))
    return list(representatives.values())


def sturm_bound(left, right):
    level = 1
    for form in (left, right):
        inverse = form**-1
        for j in range(3):
            for i in range(j+1):
                entry = inverse[i, j] / (4 if i == j else 2)
                level = lcm(level, int(entry.denominator()))
    if level % 4:
        raise ArithmeticError("Expected a ternary theta level divisible by four.")
    lift_level = level//2
    index = lift_level
    for p in map(int, pari.factor(lift_level)[0]):
        index = index//p*(p+1)
    return index//6


def shell_batch(binary, cap, cases):
    if not cases:
        return []
    inputs = [" ".join(map(str, (*parent, *(n for edge in edges for point in edge
                                           for n in point))))
              for parent, edges, _ in cases]
    result = subprocess.run([binary, str(cap)], input="\n".join(inputs)+"\n",
                            capture_output=True, check=True, text=True)
    rows = [tuple(map(int, row.split())) for row in result.stdout.splitlines()]
    if len(rows) != len(cases):
        raise ArithmeticError("Shell enumeration returned an unexpected number of cases.")
    return rows


def primitive_square_counts(ordinary):
    primitive = ordinary.copy()
    for degree in range(1,len(primitive)):
        for multiple in range(2*degree,len(primitive),degree):
            primitive[multiple] -= primitive[degree]
    return primitive


def neighbour_comparison(binary, pending, field):
    certificates = {}
    theta_jobs = {}
    ordinary_witnesses = {}
    for index, (parent, edges, _) in pending:
        det, a, b, c, u, v, w = parent
        gram = ((a,u,v),(u,b,w),(v,w,c))
        g = gram_matrix(gram)
        points, _ = normal_partition(gram,field)
        local = dict(zip(points,neighbours(g,3),strict=True))
        left, right = ([local[x] for x in edge] for edge in edges)
        if any(left[0].qfisom(right[j]) and left[1].qfisom(right[1-j]) for j in (0,1)):
            certificates[index] = ("neighbour_isometry", 0, 0)
            continue
        # Exact ordinary coefficients detect equality failure before restriction to squares.
        ordinary_counts = [q.qfrep(128) for q in (*left,*right)]
        ordinary_witnesses[index] = next(
            (n for n,(a,b,c,d) in enumerate(zip(*ordinary_counts,strict=True),1) if a+b != c+d),
            None,
        )
        bound = sturm_bound(g,g)
        jobs = []
        for q in (*left,*right):
            if int(q.matdet()) != det or sturm_bound(g,q) != bound:
                raise ArithmeticError("A level-three neighbour changed its genus invariants.")
            transform = q.qflllgram()
            reduced = transform.mattranspose()*q*transform
            a,b,c = (int(reduced[j,j]) for j in range(3))
            u,v,w = (int(reduced[i,j]) for i,j in ((0,1),(0,2),(1,2)))
            jobs.append((det,a,b,c,u,v,w,bound))
        theta_jobs[index] = jobs
    jobs = dict.fromkeys(job for group in theta_jobs.values() for job in group)
    if jobs:
        cap = max(max(job[-1] for job in jobs),3)
        if cap > 8192:
            raise ArithmeticError("The parent Sturm bound exceeds the compiled shell limit.")
        inputs = "\n".join(" ".join(map(str,job)) for job in jobs)+"\n"
        computed = subprocess.run([binary,str(cap),"--theta"], input=inputs,
                                  capture_output=True, check=True, text=True)
        counts = {job: list(map(int,row.split()))
                  for job,row in zip(jobs,computed.stdout.splitlines(),strict=True)}
        for index, group in theta_jobs.items():
            differences = (a+b-c-d for a,b,c,d in zip(*(counts[job] for job in group),strict=True))
            first = next(((3*m,d) for m,d in enumerate(differences,1) if d),None)
            certificates[index] = ("neighbour_sturm", *(first or (0,0)))
    return certificates, ordinary_witnesses


def certify_case(case, row):
    parent, edges, obstruction = case
    det, a, b, c, u, v, w = parent
    delta, excess, witness, _, checked = row
    gram = [[a,u,v],[u,b,w],[v,w,c]]
    forms = [gram_matrix(sublattice(gram, cross(*edge), 3)) for edge in edges]
    if any(int(g.matdet()) != 9*det for g in forms) or forms[0].qfisom(forms[1]):
        raise ArithmeticError("Unexpected determinant or isometry in the pair.")
    if any(int(g.qfrep(1)[0]) for g in forms):
        raise ArithmeticError("The pair contains a product.")
    replay_bound = max(delta,81)
    counts = [g.qfrep(replay_bound**2) for g in forms]
    differences = [0] + [2*int(counts[0][m*m-1]-counts[1][m*m-1])
                         for m in range(1,replay_bound+1)]
    if any(differences[:delta]) or differences[delta] != excess:
        raise ArithmeticError("Independent PARI shell counts disagree.")
    if witness and (gcd(witness,6*det) != 1 or
                    any(not g.qfrep(witness**2)[witness**2-1] for g in forms)):
        raise ArithmeticError("Independent realization witness check failed.")
    points, _ = normal_partition(gram,residue_field(-4,3))
    local = dict(zip(points,neighbours(gram_matrix(gram),3),strict=True))
    left, right = [[local[point].qfrep((replay_bound//3)**2) for point in edge] for edge in edges]
    differences_at_parent = [0] + [
        2*int(left[0][m*m-1]+left[1][m*m-1]-right[0][m*m-1]-right[1][m*m-1])
        for m in range(1,replay_bound//3+1)
    ]
    if any(differences_at_parent[:delta//3]) or differences_at_parent[delta//3] != excess:
        raise ArithmeticError("Exact neighbour theta counts disagree with the first difference.")
    primitive_child = primitive_square_counts(differences)
    primitive_parent = primitive_square_counts(differences_at_parent)
    for m in range(1,replay_bound+1):
        predicted = primitive_parent[m//3]-primitive_child[m//3] if m%3 == 0 else 0
        if primitive_child[m] != predicted:
            raise ArithmeticError("Independent primitive counts disagree with the transfer.")
    return {"parent_gram": gram, "child_grams": [matrix_entries(g) for g in forms],
            "gram_determinant": 9*det, "parent_squared_minima": [a,b,c],
            "first_difference": delta, "signed_excess": excess,
            "realization_witness": witness or None, "local_obstruction": obstruction or None,
            "sturm_bound": sturm_bound(*forms)}


def check_equal_family(field):
    edges = (((0,1,2),(1,1,0)),((1,0,2),(1,1,2)))
    results = []
    for t in (0,1,2,4,5):
        a = 5+36*t
        parent = ((a,1-36*t,2),(1-36*t,a,2),(2,2,8))
        g = gram_matrix(parent)
        if int(g.matdet()) != 160*(1+18*t):
            raise ArithmeticError("The full-square family determinant changed.")
        points, _ = normal_partition(parent,field)
        local = dict(zip(points,neighbours(g,3),strict=True))
        if any(not local[p].qfisom(local[q]) for p,q in zip(*edges)):
            raise ArithmeticError("An isometric neighbour pair is missing.")
        children = [gram_matrix(sublattice(parent,cross(*edge),3)) for edge in edges]
        expected = 5+36*t
        first = [2*int(a-b) for a,b in zip(*(q.qfrep(expected) for q in children))]
        if any(first[:-1]) or first[-1] != -2:
            raise ArithmeticError("The full-square family theta witness changed.")
        results.append({"t": t, "child_determinant": 1440*(1+18*t),
                        "first_ordinary_theta_difference": expected})
    return results


def census(limit, initial_cap):
    compiler = ["g++", "-O3", "-std=c++17", "-Wall", "-Wextra", "-Werror"]
    with TemporaryDirectory() as directory:
        parents_bin, shells_bin = (str(Path(directory)/name) for name in ("parents", "shells"))
        for program, binary in (("parents", parents_bin), ("shells", shells_bin)):
            subprocess.run([*compiler, str(HERE/f"{program}.cpp"), "-o", binary], check=True)
        enumerated = subprocess.check_output([parents_bin,str(limit)], text=True).splitlines()
        by_minima = defaultdict(list)
        prime_cache = {}
        field = residue_field(-4,3)
        cases = []
        for line in enumerated:
            parent = tuple(map(int,line.split()))
            det, a, b, c, u, v, w = parent
            gram = ((a,u,v),(u,b,w),(v,w,c))
            g = gram_matrix(gram)
            key = (det,a,b,c)
            if any(g.qfisom(other) for other in by_minima[key]):
                continue
            by_minima[key].append(g)
            if det not in prime_cache:
                prime_cache[det] = tuple(map(int,pari.factor(det)[0]))
            obstruction = local_obstruction(parent, prime_cache[det])
            cases.extend((parent, edges, obstruction) for edges in secant_pairs(gram,g,field))
        rows = shell_batch(shells_bin, initial_cap, cases)
        retry = [(i,case) for i,(case,row) in enumerate(zip(cases,rows))
                 if not case[2] and not row[0] and not row[3]]
        certificates, ordinary_witnesses = neighbour_comparison(shells_bin,retry,field)
        for i, (_,delta,excess) in certificates.items():
            _, _, witness, is_product, checked = rows[i]
            if delta and delta <= checked:
                raise ArithmeticError("A neighbour witness contradicts the initial shell enumeration.")
            rows[i] = delta,excess,witness,is_product,checked

    archive = json.loads((HERE.parent/"certificate.json").read_text())["level_three"]
    bins = defaultdict(Counter)
    ordinary_witness_histogram = Counter()
    ranked = defaultdict(list)
    checks = []
    for i,(case,row) in enumerate(zip(cases,rows)):
        parent, edges, obstruction = case
        det, _, _, c, *_ = parent
        delta, excess, witness, is_product, checked = row
        if witness and obstruction:
            raise ArithmeticError("An actual square contradicts a local obstruction.")
        if is_product and obstruction:
            raise ArithmeticError("A product contradicts a local obstruction.")
        band = bins[9*det//10000*10000]
        band["all_pairs"] += 1
        if is_product:
            band["product_pairs"] += 1
        elif obstruction:
            band["locally_obstructed_pairs"] += 1
            if delta:
                ranked["obstructed_by_minimum"].append((Fraction(delta**2,c), case, row))
        else:
            band["jacobian_pairs"] += 1
            if not delta:
                band["jacobian_square_isospectral"] += 1
                method = certificates[i][0]
                band[method] += 1
                if method == "neighbour_sturm":
                    witness_norm = ordinary_witnesses[i]
                    ordinary_witness_histogram[str(witness_norm) if witness_norm else "unseen"] += 1
            else:
                band["jacobian_separated"] += 1
                band["max_first_difference"] = max(band["max_first_difference"],delta)
                ranked["by_determinant"].append((Fraction(delta**2,9*det), case, row))
                ranked["by_minimum"].append((Fraction(delta**2,c), case, row))
                ranked["by_degree"].append((delta, case, row))
            if 9*det in {36*abs(level["D"]) for level in archive}:
                checks.append((case,row))

    # Each archived pair belongs to the larger search. Exact isometries match it.
    archived_pairs = 0
    for level in archive:
        for twin in level["twin_pairs"]:
            left, right = (gram_matrix(twin[side]) for side in ("left","right"))
            det = int(left.matdet())
            if det > limit:
                continue
            for case, row in checks:
                if det != 9*case[0][0] or row[0] != twin["first_degree"]:
                    continue
                _, a, b, c, u, v, w = case[0]
                gram = ((a,u,v),(u,b,w),(v,w,c))
                l, r = [gram_matrix(sublattice(gram,cross(*edge),3)) for edge in case[1]]
                if ((l.qfisom(left) and r.qfisom(right)) or (l.qfisom(right) and r.qfisom(left))):
                    archived_pairs += 1
                    break
            else:
                raise ArithmeticError(f"An archived pair at determinant {det} is missing.")

    maxima = {key: certify_case(case,row)
              for key, values in ranked.items() if values
              for _, case, row in [max(values)]}
    smallest_equal = min(((cases[i],kind) for i,(kind,delta,_) in certificates.items()
                          if not delta), default=None)
    sample = None
    if smallest_equal:
        (parent, edges, _), kind = smallest_equal
        det, a,b,c,u,v,w = parent
        gram = ((a,u,v),(u,b,w),(v,w,c))
        g = gram_matrix(gram)
        children = [gram_matrix(sublattice(gram,cross(*edge),3)) for edge in edges]
        if children[0].qfisom(children[1]):
            raise ArithmeticError("The full-square sample is isometric.")
        theta = [child.qfrep(64) for child in children]
        square_counts = children[0].qfrep(initial_cap**2)
        sample = {"gram_determinant": 9*det, "parent_gram": gram,
                  "child_grams": [matrix_entries(q) for q in children],
                  "secants": edges, "proof": kind, "parent_sturm_bound": sturm_bound(g,g),
                  "coprime_square_witness": next((m for m in range(1,initial_cap+1)
                                                  if gcd(m,6*det)==1 and square_counts[m*m-1]),None),
                  "first_ordinary_theta_difference": next((n for n in range(1,65)
                                                           if theta[0][n-1] != theta[1][n-1]),None)}
    return {"determinant_limit": limit, "initial_degree_cap": initial_cap,
            "reduced_matrices": len(enumerated),
            "parent_isometry_classes": sum(map(len,by_minima.values())),
            "nonisometric_pairs": len(cases), "archived_pairs_matched": archived_pairs,
            "determinant_bands": dict(sorted(bins.items())), "maxima": maxima,
            "ordinary_neighbour_witness_norms": dict(sorted(ordinary_witness_histogram.items())),
            "smallest_square_isospectral": sample, "equal_family_replay": check_equal_family(field)}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--determinant-limit", type=int, default=100000)
    parser.add_argument("--initial-cap", type=int, default=128)
    parser.add_argument("--verify", type=Path)
    options = parser.parse_args()
    if not 21 <= options.initial_cap <= 4096:
        parser.error("require 21 <= initial cap <= 4096")
    result = census(options.determinant_limit,options.initial_cap)
    data = json.dumps(result,indent=2)
    if options.verify and json.loads(data) != json.loads(options.verify.read_text()):
        raise ArithmeticError("The census differs from the stored record.")
    print(data)


if __name__ == "__main__":
    main()
