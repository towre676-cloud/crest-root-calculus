"""Evaluate PARI's generic sextic invariant polynomials over a number field."""


def absolute_invariants(pari, coefficients):
    generic = pari("genus2igusa(a6*x^6+a5*x^5+a4*x^4+a3*x^3+a2*x^2+a1*x+a0)")
    j2, j4, j6, _, j10 = [
        pari.substvec(form, [f"a{i}" for i in range(7)], coefficients) for form in generic
    ]
    if j2 == 0 or j10 == 0:
        raise ArithmeticError("The sextic lies outside the smooth J2 != 0 chart.")
    return j2**5/j10, j2**3*j4/j10, j2**2*j6/j10
