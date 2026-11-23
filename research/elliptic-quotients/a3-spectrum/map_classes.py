"""Count all geometric elliptic-map classes of a given degree on the a=3 quartic.

Classes are considered up to isomorphism of the target and translation.
The optional --directions switch prints primitive vectors, their optimal
degrees, and the number of target kernels contributing at the requested degree.
"""
import argparse
from itertools import product
import json
from math import gcd, isqrt


def optimal_degree(a, b, c):
    if gcd(a, b, c) != 1:
        raise ValueError("The direction must be a primitive integer triple.")
    kernel_order = 2 if b % 2 and c % 2 else 1
    return 2*(a*a+b*b+c*c)//kernel_order


def sum_divisors(n):
    total = 0
    for i in range(1, isqrt(n)+1):
        if n % i == 0:
            total += i
            if i*i != n:
                total += n//i
    return total


def count_classes(degree, list_directions=False):
    if type(degree) is not int or degree < 1:
        raise ValueError("The degree must be a positive integer.")
    # An optimal direction of degree <=d has squared norm <=d.
    bound = isqrt(degree)
    total, optimal, rational_optimal = 0, 0, 0
    directions = []
    for vector in product(range(-bound,bound+1), repeat=3):
        if gcd(*vector) != 1 or next(q for q in vector if q) < 0:
            continue
        base = optimal_degree(*vector)
        if degree % base:
            continue
        weight = sum_divisors(degree//base)
        total += weight
        a,b,c = vector
        rational = c==0 or (a,b)==(0,0)
        if base == degree:
            optimal += 1
            rational_optimal += int(rational)
        if list_directions:
            directions.append({"primitive_vector": list(vector), "optimal_degree": base,
                               "optimal_class_field": "Q" if rational else "Q(sqrt(-2))",
                               "isogeny_kernels_at_requested_degree": weight})
    answer = {"degree":degree, "all_geometric_classes":total,
              "geometrically_primitive_classes":optimal,
              "primitive_classes_defined_over_Q":rational_optimal}
    if list_directions:
        answer["directions"] = directions
    return answer


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("degree",type=int)
    parser.add_argument("--directions", action="store_true")
    args = parser.parse_args()
    try:
        print(json.dumps(count_classes(args.degree, args.directions), indent=2))
    except ValueError as error:
        parser.exit(2, f"Input error: {error}\n")


if __name__ == "__main__":
    main()
