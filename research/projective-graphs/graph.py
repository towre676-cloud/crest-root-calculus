"""Compute a projective graph and its bigraded Hilbert polynomial over Q.

The source ideal is assumed to define an integral projective scheme. Supply
homogeneous equations and a nonzero tuple of equal-degree homogeneous forms.
The CLI reads sparse homogeneous polynomials. Each term is a list of exponents
in source-variable order, followed by an integer or a rational string.
"""
import argparse
from collections import defaultdict
from fractions import Fraction
from itertools import chain
import json
from math import factorial
from pathlib import Path

import sympy as sp


def graph_basis(source, target, equations, coordinates, parameters=()):
    """Return the substitution kernel and whether the full basis specializes.

    Parameters are polynomial indeterminates over Q. If every leading
    coefficient in the source/target/eliminating variables is the constant 1,
    specialization preserves the computed basis at every characteristic-zero
    parameter value. Integrality of the specialized source must still be proved
    before interpreting its kernel as a prime projective graph.
    """
    source, target, parameters = tuple(source), tuple(target), tuple(parameters)
    variables = source + target + parameters
    if len(set(variables)) != len(variables) or min(len(source), len(target)) < 2:
        raise ValueError("Use distinct variables and projective spaces of dimension at least one.")
    if len(coordinates) != len(target):
        raise ValueError("The number of forms must equal the number of target coordinates.")
    polynomials = [sp.Poly(f, *source, *parameters, domain=sp.QQ)
                   for f in chain(equations, coordinates)]
    degree_sets = [{sum(term[:len(source)]) for term in p.monoms()} if not p.is_zero else set()
                   for p in polynomials]
    if any(len(degrees) > 1 for degrees in degree_sets):
        raise ValueError("Every equation and coordinate must be homogeneous in the source variables.")
    map_degrees = set().union(*degree_sets[len(equations):])
    if len(map_degrees) != 1:
        raise ValueError("Nonzero map coordinates must have the same homogeneous degree.")
    source_basis = sp.groebner(equations, *source, *parameters, order="lex", domain=sp.QQ)
    if not any(source_basis.reduce(sp.sympify(h))[1] for h in coordinates):
        raise ValueError("All map coordinates vanish in the supplied source ring.")

    scale = sp.Dummy("scale")
    relations = list(equations) + [y - scale*h for y, h in zip(target, coordinates)]
    full = sp.groebner(relations, scale, *variables, order="lex", domain=sp.QQ)
    graph = [f for f in full if scale not in f.free_symbols]
    # Parameter variables are ordered last. A full leading term free of
    # parameters implies the leading coefficient in Q[parameters] is 1.
    specialize = not parameters or all(
        not any(p.monoms()[0][-len(parameters):]) for p in full.polys)
    return graph, specialize


def hilbert_polynomial(graph, source, target):
    """Compute the eventual Hilbert polynomial of the bihomogeneous quotient.

    Each initial monomial defines an orthant of multiples. Inclusion/exclusion
    on their least common multiples gives the exact Hilbert numerator. The
    denominator is (1-S)^len(source) (1-T)^len(target).
    """
    source, target = tuple(source), tuple(target)
    variables = source + target
    basis = sp.groebner(graph, *variables, order="lex", domain=sp.QQ)
    width = len(source)
    zero = (0,)*len(variables)
    numerator = {zero: 1}
    for polynomial in basis.polys:
        exponents = polynomial.monoms()
        degrees = {(sum(mon[:width]), sum(mon[width:])) for mon in exponents}
        if len(degrees) != 1:
            raise ValueError("Expected bihomogeneous graph equations.")
        lead = exponents[0]
        updated = defaultdict(int, numerator)
        for monomial, coefficient in numerator.items():
            multiple = tuple(max(a, b) for a, b in zip(lead, monomial))
            updated[multiple] -= coefficient
        numerator = {key: value for key, value in updated.items() if value}

    shifts = defaultdict(int)
    for monomial, coefficient in numerator.items():
        shifts[(sum(monomial[:width]), sum(monomial[width:]))] += coefficient
    m, n = sp.symbols("m n")
    def contribution(variable, degree, count):
        return sp.prod(variable-degree+j for j in range(1, count))/factorial(count-1)
    return sp.expand(sum(value*contribution(m, i, width)*contribution(n, j, len(target))
                         for (i, j), value in shifts.items() if value))


def read_sparse(polynomial, variables):
    if not isinstance(polynomial, list):
        raise ValueError("Each polynomial must be a list of terms.")
    total = 0
    for term in polynomial:
        if not isinstance(term, list) or len(term) != 1 + len(variables):
            raise ValueError("Each term needs one exponent per source variable and a coefficient.")
        if any(type(degree) is not int or degree < 0 for degree in term[:-1]):
            raise ValueError("Exponents must be nonnegative integers.")
        if type(term[-1]) not in (int, str):
            raise ValueError("Each coefficient must be an integer or a rational string.")
        coefficient = Fraction(term[-1])
        total += sp.Rational(coefficient.numerator, coefficient.denominator) * sp.prod(
            variable**degree for variable, degree in zip(variables, term[:-1]))
    return sp.sympify(total)


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("input", type=Path)
    args = parser.parse_args()
    try:
        data = json.loads(args.input.read_text(encoding="utf-8"))
        width = data["source_coordinates"]
        if type(width) is not int or width < 2:
            raise ValueError("source_coordinates must be an integer of at least two.")
        source = sp.symbols(f"x:{width}")
        forms = [read_sparse(form, source) for form in data["map"]]
        target = sp.symbols(f"y:{len(forms)}")
        equations = [read_sparse(form, source) for form in data["equations"]]
        graph, _ = graph_basis(source, target, equations, forms)
        hilbert = hilbert_polynomial(graph, source, target)
        print(json.dumps({"source": list(map(str, source)), "target": list(map(str, target)),
                          "graph_basis": list(map(str, graph)),
                          "bigraded_hilbert_polynomial": str(hilbert)}, indent=2))
    except (KeyError, ValueError, TypeError, OSError, ZeroDivisionError) as error:
        parser.exit(2, f"Input error: {error}\n")


if __name__ == "__main__":
    main()
