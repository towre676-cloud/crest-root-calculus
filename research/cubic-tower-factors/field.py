"""Vector operations in a finite field using logs and a polynomial-basis table."""
from itertools import product

import numpy as np
import sympy as sp


def coefficients(value, p, degree):
    return [(value // p**j) % p for j in range(degree)]


def multiply(a, b, p, poly):
    """Polynomial-basis product; a can be a vector of encoded elements."""
    degree = len(poly)
    aa, bb = coefficients(a, p, degree), coefficients(b, p, degree)
    terms = [0] * (2*degree - 1)
    for j, v in enumerate(aa):
        for k, w in enumerate(bb):
            terms[j+k] = terms[j+k] + v*w
    for j in range(len(terms)-1, degree-1, -1):
        leading = terms[j] % p
        for k, v in enumerate(poly):
            terms[j-degree+k] = terms[j-degree+k] - leading*v
    return sum((terms[j] % p) * p**j for j in range(degree))


def power(a, exponent, p, poly):
    result = 1
    while exponent:
        if exponent & 1:
            result = multiply(result, a, p, poly)
        a = multiply(a, a, p, poly)
        exponent >>= 1
    return result


def primitive_modulus(p, degree):
    """A degree-m polynomial with X of exact order p^m-1 gives a field."""
    if degree == 1:
        return (0,), int(sp.primitive_root(p))
    order = p**degree - 1
    factors = list(sp.factorint(order))
    # Try sparse moduli first. The second pass includes all primitive polynomials.
    sparse = ((c, b) + (0,)*(degree-2) for b in range(p) for c in range(1, p))
    for candidates in (sparse, product(range(p), repeat=degree)):
        for poly in candidates:
            if poly[0] and power(p, order, p, poly) == 1 and all(
                    power(p, order//r, p, poly) != 1 for r in factors):
                return poly, p
    raise ArithmeticError("No primitive modulus found")


class Field:
    """Logs 0,...,q-2 represent nonzero elements. The field zero has log -1."""
    def __init__(self, p, degree, max_field_order=3_000_000):
        if p < 5 or not sp.isprime(p) or degree < 1:
            raise ValueError("Use a prime at least 5 and a positive extension degree")
        if p**degree > max_field_order:
            raise ValueError(f"Field order {p**degree} exceeds the table budget {max_field_order}")
        if 4*degree**2*p**2 >= 2**63 or p**degree >= 2**62:
            raise ValueError("Field products exceed the supported 64-bit arithmetic")
        self.p, self.degree = p, degree
        self.q = p**degree
        self.order = self.q-1
        modulus, generator = primitive_modulus(p, degree)
        self.modulus = [int(c) for c in modulus]
        self.exp = np.empty(self.order, dtype=np.int64)
        self.exp[0] = 1
        size, jump = 1, generator
        while size < self.order:
            count = min(size, self.order-size)
            # Small batches bound scratch memory in large extensions.
            for j in range(0, count, 32768):
                end = min(j+32768, count)
                self.exp[size+j:size+end] = multiply(self.exp[j:end], jump, p, modulus)
            jump = multiply(jump, jump, p, modulus)
            size *= 2
        self.log = np.full(self.q, -1, dtype=np.int64)
        self.log[self.exp] = np.arange(self.order)
        if np.any(self.log[1:] == -1) or self.log[0] != -1:
            raise ArithmeticError("Multiplicative group failed to enumerate the field")
        self.zech = self.log[self.translate(self.exp, 1)]

    def translate(self, encoded, constant):
        return encoded - encoded % self.p + (encoded + constant) % self.p

    def add(self, a, b):
        lift = self.zech[(b-a) % self.order]
        c = np.where(lift == -1, -1, (a+lift) % self.order)
        return np.where(a == -1, b, np.where(b == -1, a, c))

    def shift(self, logs, exponent):
        """Multiply by the field element g^exponent."""
        return np.where(logs == -1, -1, (logs + exponent) % self.order)
