"""Reconstruct the spectra, verify the universal maps and check the specialization family."""
import json
from math import prod
from pathlib import Path

import sympy as sp

from descent import HERE, ODD_ISOGENY_PRIMES, XYZ, analyze, curve_trace, expr, parse_quartic


def require(test,message):
    if not test:
        raise ArithmeticError(message)


def verify_maps():
    a,b,c,m,n,o,t,s=sp.symbols("a b c m n o t s")
    F=a*s**4+b*t**4+c+m*s*s*t*t+n*s*s+o*t*t
    A,B,C=m*m-4*a*b,2*m*n-4*a*o,n*n-4*a*c
    w=2*a*s*s+m*t*t+n
    require(sp.expand(w*w-(A*t**4+B*t*t+C)-4*a*F)==0,"Sign-quotient identity")
    x,y=A*t*t,A*t*w
    require(sp.expand(y*y-(x**3+B*x*x+A*C*x)-4*a*A*A*t*t*F)==0,"Pointed map identity")
    # Generic 2-isogeny on the pointed target and its dual curve.
    x,y,r,u=sp.symbols("x y r u")
    X=x+r+u/x
    Y=y*(1-u/x**2)
    new=sp.together(Y*Y-(X**3-2*r*X**2+(r*r-4*u)*X)).as_numer_denom()[0]
    require(sp.rem(new,y*y-(x**3+r*x*x+u*x),y)==0,"Dual 2-isogeny identity")


def independent_count(A,B,C,p):
    # w^2=A t^4+B t^2+C as a weighted projective double cover of P^1.
    sums=[A*t**4+B*t*t+C for t in range(p)]+[A]
    return sum(sum(w*w%p==k%p for w in range(p)) for k in sums)


def family_form(t):
    return {"terms":[[4,0,0,1],[0,4,0,3],[0,0,4,63],[2,2,0,27*t]]}


def main():
    verify_maps()
    sample={}
    for case in ("diagonal","deformation","mixed"):
        data=json.loads((HERE/(case+".json")).read_text())
        result=analyze(data,recover=case=="mixed")
        require(result.get("all_rational_elliptic_degrees")==["4*n^2","8*n^2"],case)
        require(result.get("all_rational_genus_one_degrees")==["2*n^2","4*n^2"],case)
        require(result["source_index"]==4,case)
        for p,values in result["isogeny_sieve"]["frobenius_traces"].items():
            for model,t in zip(result["sign_quotients"],values):
                if t is not None:
                    require(p+1-independent_count(*model["binary_quartic"],p)==t,"Torsor point count")
                    require(curve_trace(model["degree_four_target"],p)==t,"Isogeny point count")
        sample[case]=result
    require(sample["diagonal"]["fermat_cubic_direction_orbit"]["galois_orbit_size"]==16,
            "Geometric cubic directions")
    require(sample["mixed"]["input_index_witness"] is None,"Mixed-frame control")
    require(sample["mixed"]["frame_index_witness"]["prime"]==3,"Recovered index obstruction")
    original=parse_quartic(json.loads((HERE/"mixed.json").read_text()))
    change=sp.Matrix(sample["mixed"]["old_coordinates_in_new_frame"])
    expected=expr(parse_quartic({"terms":sample["mixed"]["even_frame_terms"]}))
    ratio=sp.cancel(expr(original).subs(dict(zip(XYZ,change*sp.Matrix(XYZ))),simultaneous=True)/expected)
    require(ratio.is_Rational and ratio!=0,"Frame substitution")
    # Verify full scope of the arithmetic progression using integral polynomial coefficients.
    primes=sorted(sample["deformation"]["isogeny_sieve"]["frobenius_traces"])
    modulus=prod(set(primes)|{7})
    require(modulus==5*7*11*17*19*23*29,"Family comparison primes")
    control=analyze(family_form(1+modulus))
    require(control["all_rational_elliptic_degrees"]==["4*n^2","8*n^2"],"Second family member")
    for p,traces in sample["deformation"]["isogeny_sieve"]["frobenius_traces"].items():
        repeated=[curve_trace(m["jacobian"],p) for m in control["sign_quotients"]]
        require(repeated==traces,"Uniform local witnesses")
    missing=analyze(family_form(1),bound=5)
    require("all_rational_elliptic_degrees" not in missing,"Incomplete local evidence")
    weak=analyze({"terms":[[4,0,0,1],[0,4,0,1],[0,0,4,-1]]})
    require("source_index" not in weak,"Source with a rational point")
    general=analyze({"terms":[*family_form(1)["terms"],[3,0,1,27]]})
    require(general["source_index"]==4 and "sign_quotients" not in general,"Nondiagonal congruence class")
    # The only possibly tied terms for e=2 mod 4 require a nonsquare.
    for e in range(1,17):
        for f in range(1,9):
            if e*f%4 and e%2:
                require(len({0,e%4,2*e%4})==3,"Odd ramification")
            elif e*f%4:
                require(f%2==1 and e%4==2,"Even ramification")
    print(json.dumps({"verified_cases":list(sample),"complete_genus_one_degrees":["2*n^2","4*n^2"],
                      "complete_pointed_degrees":["4*n^2","8*n^2"],"progression_modulus":modulus,
                      "good_prime_traces":sample["deformation"]["isogeny_sieve"]["frobenius_traces"],
                      "local_non_diagonal_obstruction":general["input_index_witness"],
                      "diagonal_cubic_orbit":sample["diagonal"]["fermat_cubic_direction_orbit"]},indent=2))


if __name__=="__main__":
    main()
