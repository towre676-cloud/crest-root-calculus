"""Identify the residual Hecke eigensystems and verify their twisted central values."""
import argparse
import json
from fractions import Fraction
from math import gcd, lcm
import sys

from exact import ambient_genus, depleted_sturm, fiber, kernel, pari, residue_field, sturm_bound
from neighbours import Closure


# Weierstrass models from the indicated Cremona classes.
CURVES = {
    (-8,13): (("676b1",[0,1,0,-4,-12]), ("338c1",[1,0,0,81,467])),
    (-3,17): (("867e1",[0,1,1,6,-1]),),
}


def lift_coefficients(counts, relation, D, ell, t, bound):
    coefficients = [0]*(bound+1)
    square_class = [0]+[sum(2*c*a for c,a in zip(relation,column)) for column in zip(*counts)]
    for d in range(1,bound+1):
        if gcd(d,2*D*ell)==1:
            character = int(pari.kronecker(D*t,d))
            for m in range(1,bound//d+1):
                coefficients[d*m] += character*square_class[m]
    return coefficients


def oldform_weights(coefficients, curve, target, conductor):
    if target%conductor:
        raise ArithmeticError("The newform level does not divide the lift level.")
    bound = len(coefficients)-1
    dilations = list(map(int,pari.divisors(target//conductor)))
    a = [0]+list(map(int,curve.ellan(bound)))
    remainder = coefficients.copy()
    weights = {}
    for d in dilations:
        if d>bound:
            raise ArithmeticError("The comparison cutoff does not see every dilation.")
        c = remainder[d]
        if c:
            weights[d] = c
            for m in range(1,bound//d+1):
                remainder[d*m] -= c*a[m]
    if any(remainder) or not weights:
        raise ArithmeticError("The Shimura lift differs from the claimed elliptic oldforms.")
    return weights


def run(D,ell):
    records = {"split":[],"nonsplit":[]}
    field = residue_field(D,ell)
    for Q,group in ambient_genus(D):
        fiber(Q,group,field,records)
    prime = 3 if D==-8 else 5
    closure = Closure([Q for Q,_ in records["nonsplit"]],(prime,))
    n = len(closure.forms)
    if n != len(records["nonsplit"]):
        raise ArithmeticError("The exhaustive genus list must be closed under neighbours.")
    _,depleted_bound = depleted_sturm(D,ell)
    K = kernel(closure.forms,D,ell,depleted_bound)
    A = pari.matrix(n,n,[row.count(i) for row in closure.rows[prime] for i in range(n)]).mattranspose()
    restriction = K.matinverseimage(A*K)
    if len(K) != len(CURVES[D,ell]) or K*restriction != A*K:
        raise ArithmeticError("The claimed residual Hecke module changed.")
    target = 2*abs(D)*ell*ell
    bound = sturm_bound(target,D,ell)
    t = 5
    counts = []
    for q in closure.forms:
        reps = q.qfrep(t*bound**2)
        counts.append([int(reps[t*m*m-1]) for m in range(1,bound+1)])
    del reps
    found = set()
    results = []
    for label,model in CURVES[D,ell]:
        curve = pari.ellinit(model)
        conductor = int(curve.ellglobalred()[0])
        trace = int(curve.ellap(prime))
        eigenline = (restriction-trace*pari.matid(len(K))).matker()
        if len(eigenline)!=1 or trace in found:
            raise ArithmeticError("The exact Hecke eigenspace is not the claimed line.")
        found.add(trace)
        rational = [Fraction(str(c)) for c in K*eigenline[0]]
        denominator = lcm(*(c.denominator for c in rational))
        v = [int(c*denominator) for c in rational]
        if sum(v)!=0:
            raise ArithmeticError("Expected a cuspidal theta combination.")
        lift = lift_coefficients(counts,v,D,ell,t,bound)
        weights = oldform_weights(lift,curve,target,conductor)
        twisted = pari.ellinit(curve.elltwist(D))
        space,symbol = twisted.msfromell(1)
        central_ratio = pari.mseval(space,symbol,[pari("oo"),0])
        root_number = int(twisted.ellrootno())
        if central_ratio!=0 or root_number!=1:
            raise ArithmeticError("The expected central zero and sign did not certify.")
        results.append({"label":label,"ainvs":model,"conductor":conductor,"auxiliary_prime":prime,
                        "neighbor_eigenvalue":trace,"shimura_class":t,"comparison_bound":bound,
                        "oldform_dilations":weights,"twist":D,"twisted_conductor":int(twisted.ellglobalred()[0]),
                        "twisted_root_number":root_number,"central_L_over_real_period":str(central_ratio)})
        print(label,weights,root_number,central_ratio,file=sys.stderr,flush=True)
    return {"D":D,"ell":ell,"nonsplit_kernel_dimension":len(K),"packets":results}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--discriminant",type=int,choices=(-8,-3))
    args = parser.parse_args()
    pari.allocatemem(128_000_000,2_500_000_000,silent=True)
    print(json.dumps([run(D,ell) for D,ell in CURVES if args.discriminant in (None,D)],indent=2))
