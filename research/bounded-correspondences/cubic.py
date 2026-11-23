"""Exact reduction for all genus-two maps totally ramified over elliptic infinity.

Inputs are odd quintic and cubic branch polynomials in increasing powers.
The returned single-parameter equations retain every geometric map and permit
the generic finite-incidence engine to recover all rational maps by local lifting.
"""
import sympy as sp

from finite import polynomial


X = sp.Symbol("x")
SCALE = sp.Symbol("scale")


def cubic_incidence(source, target):
    """Return equations, variables, U and V for u=U(x), v=y*V(x).

    A Poly input may use any variable. Coefficient-list inputs use ascending
    powers. The source and target must be squarefree, of degrees five and three.
    """
    source = source if isinstance(source, sp.Poly) else polynomial(source, X)
    target = target if isinstance(target, sp.Poly) else polynomial(target, X)
    if source.degree() != 5 or target.degree() != 3:
        raise ValueError("The cubic reduction requires a quintic and a cubic.")
    if sp.gcd(source, source.diff()).degree() or sp.gcd(target, target.diff()).degree():
        raise ValueError("Both branch polynomials must be squarefree.")
    x = source.gens[0]
    source_lead = source.LC()
    source_shift = source.nth(4) / (5 * source_lead)
    centered = sp.Poly(source.as_expr().subs(x, x - source_shift) / source_lead, x)
    A, B, C, D = (centered.nth(i) for i in (3, 2, 1, 0))

    target_lead = target.nth(3)
    target_quadratic = target.nth(2)
    P = target_lead * target.nth(1) - target_quadratic**2 / 3
    Q = (2 * target_quadratic**3 / 27
         - target_lead * target_quadratic * target.nth(1) / 3
         + target_lead**2 * target.nth(0))

    centered_x = x + source_shift
    base_u = centered_x**3 + 3 * A * centered_x / 7 + B / 3
    base_v = centered_x**2 + A / 7
    base_p = 4 * A**3 / 343 - B**2 / 3
    base_q = 8 * A**3 * B / 1029 + 2 * B**3 / 27
    U = sp.expand((source_lead * SCALE**2 * base_u - target_quadratic / 3)
                  / target_lead)
    V = sp.expand(source_lead * SCALE**3 * base_v / target_lead)

    if C != 12 * A**2 / 49 or D != 4 * A * B / 7:
        return [sp.Integer(1)], [SCALE], U, V
    equations = [P - source_lead**2 * base_p * SCALE**4,
                 Q - source_lead**3 * base_q * SCALE**6]
    # Smoothness of the target implies (P,Q) != (0,0). Thus scale=0
    # cannot be a solution and an extra inverse variable is unnecessary.
    return equations, [SCALE], U, V
