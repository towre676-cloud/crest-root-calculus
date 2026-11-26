"""Exact local-conductor and neighbour tests for density of separating primes."""
import argparse
import json
from fractions import Fraction
from math import gcd, isqrt
from pathlib import Path
import sys

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hecke-observation"))
from neighbours import Closure, gram_matrix, pari
from observe import form_gram, residue_field, run, scalar_orbits


def primes_up_to(limit):
    return [p for p in range(2,limit+1) if all(p % q for q in range(2,isqrt(p)+1))]


def fundamental_discriminant(squarefree):
    return squarefree if squarefree % 4 == 1 else 4*squarefree


def allowed_characters(D, ell):
    """Discriminants of quadratic twists whose conductor square divides M/2."""
    level = 2*abs(D)*ell*ell
    candidates = [-1, ell, -ell]
    return sorted(d for a in candidates if level % (d := fundamental_discriminant(a))**2 == 0)


def legendre(D, p):
    return int(pari.kronecker(D,p))


def hecke_matrices(D, ell, primes):
    field = residue_field(D,ell)
    closure = Closure([form_gram(field,a) for a in scalar_orbits(field)], primes)
    n = len(closure.forms)
    result = []
    for p in primes:
        A = pari.matrix(n,n,[r.count(i) for r in closure.rows[p] for i in range(n)])
        result.append(A)
    return closure, result


def audit_local(D, ell, primes, matrices):
    level = 2*abs(D)*ell*ell
    n = len(matrices[0])
    for i,A in enumerate(matrices):
        for B in matrices[:i]:
            if A*B != B*A:
                raise ArithmeticError("The good-prime operators must commute.")
    forbidden = []
    possible = []
    for d in allowed_characters(D,ell):
        signs = [legendre(d,p) for p in primes]
        for power in (1,2,3):
            w = [i**power for i in range(1,len(primes)+1)]
            X = sum(c*A for c,A in zip(w,matrices))
            Y = sum(c*s*A for c,s,A in zip(w,signs,matrices))
            g = pari.gcd(X.charpoly(),Y.charpoly())
            if g.poldegree() == 0:
                forbidden.append({"character":d,"weights":w})
                break
        else:
            possible.append(d)
    cm = []
    cm_excluded = []
    for d in possible:
        inert = [(p,A) for p,A in zip(primes,matrices) if legendre(d,p)==-1]
        if d<0:
            if inert and sum(A*A for _,A in inert).matdet() != 0:
                cm_excluded.append(d)
            else:
                cm.append(d)
    return {"D":D,"ell":ell,"level_bound":level,"primes":primes,"closure":n,
            "excluded_twists":forbidden,"remaining_twists":possible,
            "excluded_cm_characters":cm_excluded,"remaining_cm_characters":cm}


def krylov_partition(matrix, initial, num_forms):
    """Use exact integer observations to separate the gluing coordinates."""
    remaining = [list(range(num_forms))]
    vector = initial
    for exponent in range(len(initial)):
        next_groups = []
        for group in remaining:
            buckets = {}
            for i in group:
                buckets.setdefault(int(vector[i]), []).append(i)
            next_groups.extend(g for g in buckets.values() if len(g) > 1)
        if not next_groups:
            return exponent, []
        remaining = next_groups
        vector = matrix * vector
    return exponent, remaining


def certify_density(local, max_prime, closure, matrices):
    """Use exact finite witnesses in every remaining character class."""
    D,ell = local["D"],local["ell"]
    characters = local["remaining_twists"]
    modulus = 4*ell
    signatures = {tuple(legendre(d,r) for d in characters) for r in range(1,modulus) if gcd(r,modulus)==1}
    cm = local["remaining_cm_characters"]
    operators = dict(zip(local["primes"],matrices))
    initial = pari.Col([2*int(q.qfrep(1)[0]) for q in closure.forms])
    for d in cm:
        inert = [A for p,A in operators.items() if legendre(d,p) == -1]
        if not inert:
            return {"witnesses":[],"density_lower_bound":"0"}
        initial = sum(A*A for A in inert) * initial
    found = {}
    for p in (operators if cm else primes_up_to(max_prime)):
        if gcd(p,2*D*ell)!=1:
            continue
        key = tuple(legendre(d,p) for d in characters)
        if key in found:
            continue
        if cm:
            exponent, unresolved = krylov_partition(operators[p],initial,closure.nseeds)
            if unresolved:
                exponent = None
        else:
            obs = run(D,ell,p)
            exponent = None if obs["joint_unresolved"] else obs["tested_through_exponent"]
        if exponent is not None:
            found[key] = {"p":p,"signature":key,"exponent":exponent,"CM_filtered":bool(cm)}
        if signatures==set(found):
            break
    return {"witnesses":list(found.values()),"density_lower_bound":str(Fraction(len(found),len(signatures)))}


def exact_failed_primes():
    results = []
    for D,ell,p,expected in [(-3,5,17,2),(-8,7,17,11),(-8,7,29,3)]:
        closure, (matrix,) = hecke_matrices(D,ell,(p,))
        initial = pari.Col([2*int(q.qfrep(1)[0]) for q in closure.forms])
        exponent, unresolved = krylov_partition(matrix,initial,closure.nseeds)
        pairs = sum(len(g)*(len(g)-1)//2 for g in unresolved)
        if pairs != expected:
            raise ArithmeticError("The number of all-powers ties changed.")
        results.append({"D":D,"ell":ell,"p":p,"tested_through_exponent":exponent,"tied_pairs":pairs})
    return results


def split_prime_collision():
    """Prove that the same pair is tied at every split auxiliary prime of Q(sqrt(-11))."""
    D, ell, character = -3, 11, -11
    field = residue_field(D, ell)
    centers = [(0,1),(2,2)]
    if len({field.legendre(1+field.norm(a)) for a in centers}) != 1:
        raise ArithmeticError("Expected two forms in the same genus.")
    if not all(a in scalar_orbits(field) for a in centers):
        raise ArithmeticError("Expected distinct reduced scalar orbits.")
    # M=1452, target level 726 (weight-two Sturm bound 264).
    bound = 264
    r0,r1 = [gram_matrix(form_gram(field,a)).qfrep(bound**2) for a in centers]
    squares = [0] + [2*int(r0[m*m-1] - r1[m*m-1]) for m in range(1,bound+1)]
    lift = [0]*(bound+1)
    for d in range(1,bound+1):
        if gcd(d,2*abs(D)*ell) == 1:
            for m in range(1,bound//d+1):
                lift[d*m] += legendre(D,d)*squares[m]
    # Recover F from F(z)+F(2z) using the coefficient at each lower index.
    coeffs = lift.copy()
    for m in range(2,bound+1,2):
        coeffs[m] -= coeffs[m//2]
    # The Sturm bound at level 363 is 88. Use all available coefficients
    # to verify the full identity at level 726 as well as the twist sign.
    space = pari.mfinit([363,2],1)
    F = pari.mflinear(space,pari.mftobasis(space,pari(coeffs)))
    if pari.mfcoefs(F,bound) != pari(coeffs):
        raise ArithmeticError("The recovered cusp form does not give the square lift.")
    if any(a and legendre(character,m)!=-1 for m,a in enumerate(coeffs)):
        raise ArithmeticError("The odd-level cusp form is not anti-invariant under the twist.")
    return {"D":D,"ell":ell,"centers":centers,"twisting_character":character,
            "level_of_cusp_form":363,"square_lift_sturm_bound":bound,
            "cusp_coefficients_through_64":coeffs[:65],"nonseparating_character_sign":1}


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("levels", nargs="*", help="D:ell")
    parser.add_argument("--verify", type=Path, help="Recompute and check an exact receipt.")
    parser.add_argument("--operator-prime-limit", type=int, default=23)
    parser.add_argument("--witness-prime-limit", type=int, default=43)
    args=parser.parse_args()
    if bool(args.levels) == bool(args.verify):
        parser.error("Provide levels or --verify RECEIPT.")
    recorded = json.loads(args.verify.read_text()) if args.verify else None
    levels = [f"{r['D']}:{r['ell']}" for r in recorded["levels"]] if recorded else args.levels
    pari.allocatemem(32_000_000, 1_000_000_000, silent=True)
    results=[]
    for spec in levels:
        D,ell=map(int,spec.split(":"))
        ps=[p for p in primes_up_to(args.operator_prime_limit) if gcd(p,2*D*ell)==1]
        closure, matrices = hecke_matrices(D,ell,ps)
        local=audit_local(D,ell,ps,matrices)
        result={**local,**certify_density(local,args.witness_prime_limit,closure,matrices)}
        results.append(result)
        print(result,file=sys.stderr,flush=True)
    result = {"levels": results, "split_prime_collision": split_prime_collision(),
              "failed_primes": exact_failed_primes()}
    rendered = json.dumps(result, indent=2)
    if recorded is not None and json.loads(rendered) != recorded:
        raise ArithmeticError("The recomputed receipt differs.")
    print(rendered)


if __name__=="__main__":
    main()
