"""Prove density bounds by filtering the good Shimura coefficient indices."""
import argparse
from collections import defaultdict
from fractions import Fraction
import json
from math import gcd
from pathlib import Path
import sys

from coprime import CLASS_NUMBER_ONE, compare_receipt, counts_on_shell, primitive_shells, residue_field, scalar_orbits


def previously_excluded_two_adic():
    path = Path(__file__).resolve().parents[1] / "prime-density"
    small = json.loads((path / "receipt.json").read_text())["levels"]
    large = json.loads((path / "two-adic.json").read_text())
    return {r["ell"] for r in small if r["D"] == -8 and len(r["excluded_twists"]) == 2
            and all(t["character"] % 2 == 0 for t in r["excluded_twists"])} | {
        r["ell"] for r in large if set(r["gcd_degrees"].values()) == {0}}


TWO_ADIC_EXCLUDED = previously_excluded_two_adic()


def odd_prime_factors(n):
    """Primes of odd valuation."""
    factors = []
    p = 2
    while p*p <= n:
        odd = False
        while n%p == 0:
            n //= p
            odd = not odd
        if odd:
            factors.append(p)
        p += 1
    if n > 1:
        factors.append(n)
    return factors


def divisor_terms(D, bound):
    """The convolution of the constant-one arithmetic function with chi_D."""
    char = [0]*(bound+1)
    for n in range(1,bound+1,2):
        if D == -4:
            char[n] = 1 if n%4==1 else -1
        elif D == -8:
            char[n] = 1 if n%8 in (1,3) else -1
        else:
            v = pow(n%(-D),(-D-1)//2,-D)
            char[n] = -1 if v==-D-1 else v
    sums = [0]*(bound+1)
    for d in range(1,bound+1,2):
        if char[d]:
            for m in range(d,bound+1,2*d):
                sums[m] += char[d]
    terms = [[] for _ in range(bound+1)]
    for n in range(1,bound+1,2):
        for d in range(1,bound//n+1,2):
            if sums[d]:
                terms[n*d].append((n,sums[d]))
    return terms


def value_bit(field, m):
    return int(field.legendre(m)==-1) | (int(m%4==3)<<1)


def negative(character, bits):
    return bool((character & bits).bit_count()%2)


def classify_indices(D, field, factorizations, group, sectors):
    """Map each allowed coefficient index to the sign sectors it can certify."""
    # Odd D forces Steinberg at q=|D|. Imaginary characters ramified at 2
    # would require a CM newform of 2-exponent at least five.
    imaginary = [1] if D%2==0 and field.p%4==3 else []
    prime_values = {}
    for m in range(1,len(factorizations),2):
        if gcd(m, D*field.p)!=1:
            continue
        factor_signs = []
        for p in factorizations[m]:
            if p not in prime_values:
                prime_values[p] = value_bit(field,p)
            factor_signs.append(prime_values[p])
        if any(not any(negative(c, s) for s in factor_signs) for c in imaginary):
            continue
        values = value_bit(field,m)
        valid = [s for s in sectors if all(negative(c,s) or not negative(c,values) for c in group)]
        if valid:
            yield m,valid


def refine(groups, counts, degree):
    unresolved = []
    last = None
    for group in groups:
        bins = defaultdict(list)
        for a in group:
            bins[counts[a]].append(a)
        if len(bins)>1:
            lower,higher = sorted(bins)[:2]
            last = {"index":degree,"centers":[bins[lower][0],bins[higher][0]],"lift_counts":[lower,higher]}
        unresolved.extend(g for g in bins.values() if len(g)>1)
    return unresolved,last


def certify(D, field, shells, terms, factorizations):
    centers = list(scalar_orbits(field))
    traces = {a:(field.trace(a),field.trace(field.mul(a,(0,1)))) for a in centers}
    has_extra_characters = D == -8 and field.p not in TWO_ADIC_EXCLUDED
    number_sectors = 4 if has_extra_characters else 2
    characters = range(1,4) if has_extra_characters else (1,)
    partitions = [[centers] if len(centers)>1 else [] for _ in range(number_sectors)]
    last = [None]*number_sectors
    for m,valid in classify_indices(D,field,factorizations,characters,range(number_sectors)):
        if not any(partitions):
            break
        targets = {a for s in valid for group in partitions[s] for a in group}
        if not targets:
            continue
        active = sorted(targets)
        values = [0]*len(active)
        for n,weight in terms[m]:
            counts = counts_on_shell(active,traces,shells[n-1],field.p)
            for i,count in enumerate(counts):
                values[i] += weight*count
        counts = dict(zip(active,values))
        for s in valid:
            partitions[s], witness = refine(partitions[s],counts,m)
            last[s] = witness or last[s]
    success = sum(not u for u in partitions)
    return {"D":D,"ell":field.p,"forms":len(centers),"density_lower_bound":str(Fraction(success,number_sectors)),
            "sectors":[{"mask":s,"last_witness":last[s],"unresolved_pairs":sum(len(g)*(len(g)-1)//2 for g in groups),
                        "example_unresolved_pair":groups[0][:2] if groups else None}
                       for s,groups in enumerate(partitions)]}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--prime-limit",type=int,default=199)
    parser.add_argument("--degree-limit",type=int,default=1024)
    parser.add_argument("--discriminant",type=int,choices=CLASS_NUMBER_ONE,action="append")
    parser.add_argument("--level",type=int,help="Check only this gluing prime.")
    parser.add_argument("--verify",type=Path,help="Compare with the saved receipt.")
    args=parser.parse_args()
    if args.prime_limit<3 or args.degree_limit<1:
        parser.error("Use a gluing prime bound of at least 3 and a positive degree limit.")
    factorizations = [odd_prime_factors(n) for n in range(args.degree_limit+1)]
    output=[]
    for D in args.discriminant or CLASS_NUMBER_ONE:
        shells=None
        terms=divisor_terms(D,args.degree_limit)
        for p in ([args.level] if args.level else range(3,args.prime_limit+1,2)):
            try:
                field=residue_field(D,p)
            except ValueError:
                continue
            if shells is None:
                shells=primitive_shells(field,args.degree_limit)
            result=certify(D,field,shells,terms,factorizations)
            output.append(result)
            print(D,p,result["density_lower_bound"],
                  [s["unresolved_pairs"] for s in result["sectors"]],file=sys.stderr,flush=True)
    compare_receipt(output,args.verify)


if __name__ == "__main__":
    main()
