"""Exact certificates for opposite secants with unbounded normalized finite delay."""

import argparse
from itertools import product
import json
from math import gcd, isqrt
from pathlib import Path


def det(gram):
    (a,u,v), (_,b,w), (_,_,c) = gram
    return a*b*c + 2*u*v*w - a*w*w - b*v*v - c*u*u


def norm(gram, point):
    (a,u,v), (_,b,w), (_,_,c) = gram
    x,y,z = point
    return a*x*x + b*y*y + c*z*z + 2*(u*x*y + v*x*z + w*y*z)


def short_norms_are_separated(n, s, error):
    # Test the two squares surrounding n²R, for each nonsquare R.
    n_squared = n*n
    for r in range(1,s*s+1):
        if isqrt(r)**2 == r:
            continue
        value = n_squared*r
        root = isqrt(value)
        if min(value-root*root, (root+1)**2-value) <= error*r:
            return False
    return True


def direct_square_counts(gram, n, s):
    """Independent shell check using the spatial cutoff in the proof."""
    square_counts = [{}, {}]
    limit = s*s
    for x in range(-s,s+1):
        y_limit = isqrt(limit-x*x)
        for y in range(-y_limit,y_limit+1):
            z_limit = isqrt(limit-x*x-y*y)
            for z in range(-z_limit,z_limit+1):
                if not (x or y or z):
                    continue
                value = norm(gram,(x,y,z))
                degree = isqrt(value)
                if degree > n*s or degree*degree != value:
                    continue
                for counts, residual in zip(square_counts,(x-y,x+y),strict=True):
                    if residual % 3 == 0:
                        counts[degree] = counts.get(degree,0)+1
    common = {n*j: 2 for j in range(1,s+1)}
    excess = {**common, n*s: 4}
    if square_counts != [excess,common]:
        raise ArithmeticError("The independent lattice-shell comparison failed.")


def construct(k, replay_shells):
    if k < 2 or k%3 != 2:
        raise ValueError("Require k >= 2 and k congruent to 2 modulo 3.")
    y,z = 2*k,2*k*k
    s = z+1
    power = 1
    while power <= 11*s*s:
        power *= 3
    coefficient = 6*(power+1)
    error = coefficient*(1+y)**2 + 6*(1+y+s+(1+y)*s)
    perturbation = (
        (coefficient*y*y+6*(y*z-y-z), -coefficient*y+3*(1-z), 3*(1-y)),
        (-coefficient*y+3*(1-z), coefficient, 3),
        (3*(1-y), 3, 0),
    )
    start = 2*error*s*(12*s**3+1)+1
    start += (1-start)%6
    n = next(
        (j for j in range(start,start+6*(s*s+1),6)
         if short_norms_are_separated(j,s,error)),
        None,
    )
    if n is None:
        raise ArithmeticError("No witness in the proved search interval.")
    gram = tuple(tuple(value+n*n*(i==j) for j,value in enumerate(row))
                 for i,row in enumerate(perturbation))
    parent_det = det(gram)
    first_difference = n*s
    if (det(perturbation) != -54*power or gcd(n,9*parent_det) != 1
            or (n*n-error)*(s*s+1) <= first_difference**2 or n*n-error <= 1):
        raise ArithmeticError("The determinant, coprimality or outer shell test failed.")
    if any(gram[i][j]%3 != int(i==j) for i,j in product(range(3),repeat=2)):
        raise ArithmeticError("Expected the identity modulo three.")
    if any(norm(gram,w)%4 != sum(w)**2%4 for w in product(range(4),repeat=3)):
        raise ArithmeticError("The modulo-four norm condition failed.")
    if y%3 != 1 or norm(gram,(1,y,z)) != first_difference**2 or norm(gram,(0,0,1)) != n*n:
        raise ArithmeticError("The two prescribed square norms or the secant label failed.")
    if replay_shells:
        direct_square_counts(gram,n,s)
    return {
        "k": k,
        "squared_minima_lower_bound": n*n-error,
        "squared_minima_upper_bound": n*n+error,
        "parent_gram": gram,
        "child_gram_determinant": 9*parent_det,
        "common_first_degree": n,
        "first_distinguishing_degree": first_difference,
        "primitive_counts_at_first_difference": [2,0],
        "normalized_delay_target": s,
        "distinguishing_vector": [1,y,z],
        "checked_by_direct_lattice_shells": replay_shells,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--verify", type=Path)
    options = parser.parse_args()
    cases = [construct(k, k <= 8) for k in (2,5,8,20)]
    receipt = json.dumps(cases,indent=2)
    if options.verify and json.loads(receipt) != json.loads(options.verify.read_text()):
        raise ArithmeticError("The exact certificates differ from the stored receipt.")
    print(receipt)


if __name__ == "__main__":
    main()
