"""Certify the split-prime obstruction at the five odd discriminants with ell=3."""
import json
from math import gcd
from pathlib import Path
import sys

from coprime import residue_field, scalar_orbits

sys.path.insert(0, str(Path(__file__).resolve().parents[1] / "hecke-observation"))
from neighbours import gram_matrix, pari
from verify import form_gram


def witness(D, dilation):
    ell = 3
    q = -D
    target = 2*q*ell*ell
    bound = 6*(q+1)
    field = residue_field(D,ell)
    centers = [(1,0),(1,1)]
    orbits = scalar_orbits(field)
    if not all(a in orbits for a in centers) or len({field.legendre(1+field.norm(a)) for a in centers}) != 1:
        raise ArithmeticError("Expected distinct scalar orbits in the same genus.")
    first,second = [gram_matrix(form_gram(field,a)).qfrep(bound**2) for a in centers]
    square_difference = [0]+[2*int(first[m*m-1]-second[m*m-1]) for m in range(1,bound+1)]
    lift = [0]*(bound+1)
    for d in range(1,bound+1,2):
        if gcd(d,q*ell)==1:
            character = int(pari.kronecker(D,d))
            for m in range(1,bound//d+1):
                lift[d*m] += character*square_difference[m]
    coeffs = lift.copy()
    for m in range(2,bound+1,2):
        coeffs[m] -= dilation*coeffs[m//2]
    space = pari.mfinit([target//2,2],1)
    cusp = pari.mflinear(space,pari.mftobasis(space,pari(coeffs)))
    if pari.mfcoefs(cusp,bound)!=pari(coeffs):
        raise ArithmeticError("The reconstructed cusp form does not give the square lift.")
    if not any(coeffs) or any(a and int(pari.kronecker(-3,m))!=-1 for m,a in enumerate(coeffs)):
        raise ArithmeticError("The cusp form must be nonzero and anti-invariant under the twist.")
    return {"D":D,"ell":ell,"centers":centers,"square_lift_level":target,"cusp_level":target//2,
            "sturm_bound":bound,"dilation_coefficient":dilation,"cusp_coefficients_through_16":coeffs[:17]}


if __name__ == "__main__":
    pari.allocatemem(32_000_000,512_000_000,silent=True)
    print(json.dumps([witness(D, -1 if D==-7 else 1) for D in (-7,-19,-43,-67,-163)],indent=2))
