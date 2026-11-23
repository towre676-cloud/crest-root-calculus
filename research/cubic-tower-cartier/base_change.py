"""Propagate a Weil-polynomial congruence through finite-field extensions."""
from math import factorial

import sympy as sp


def powered_congruence(coefficients, modulus, extension):
    """Return det(1−T F^extension) and the guaranteed new coefficient modulus.

    Starting coefficients are in ascending order for det(1−T F). The modulus
    gains p^v_p(extension) at each prime p already dividing the supplied modulus.
    """
    if modulus < 2 or extension < 1 or coefficients[0] % modulus != 1:
        raise ValueError("Expected a positive extension, modulus >=2 and constant one")
    target = modulus
    for p in sp.factorint(modulus):
        n = extension
        while n % p == 0:
            n //= p
            target *= p
    dim = len(coefficients)-1
    precision = target*factorial(dim)
    # Compute ghost components to enough precision for exact divisions by 1,...,dim.
    h = [int(c) % modulus for c in coefficients]
    traces = [0]*(extension*dim+1)
    for j in range(1,len(traces)):
        value = j*h[j] if j<=dim else 0
        value += sum(h[i]*traces[j-i] for i in range(1,min(j-1,dim)+1))
        traces[j] = -value % precision
    result = [1]
    for k in range(1,dim+1):
        total = -sum(result[k-i]*traces[extension*i] for i in range(1,k+1)) % precision
        if total % k:
            raise ArithmeticError("Lifted polynomial failed Newton integrality")
        result.append(total//k)
        precision //= k
    return {"extension_degree":extension,"coefficient_modulus":target,
            "numerator":[v%target for v in result]}
