"""Independent symbolic rank argument and regression for the permuted-fibre inverse."""
from __future__ import annotations

import argparse
import json
from pathlib import Path

import sympy as sp

from etale_fibres import coefficients, normalize, pari, recover, recover_kernel, t, x
from reconstruct_four import certify_identity, odd_chart

HERE = Path(__file__).resolve().parent
PERMUTED = x**6-6*x**5-11*x**4+14*x**3+12*x*x-12*x+1
SPLIT = 11*(x*x-10*x+13)*(2*x*x+2*x-7)*(7*x*x+18*x+3)


def check_rank_identities():
    a, b, c, d, A = sp.symbols("a b c d A")
    # Coefficient of x in q_b*q_c modulo q_a, for each finite marked
    # critical point of a quartic polynomial with zero cubic term.
    qa, qb, qc = (x*x+2*v*x+(A+3*v*v) for v in (a, b, -a-b))
    r = sp.Poly(sp.rem(qb*qc, qa, x), x).nth(1)
    assert sp.expand(r+2*a*(a-b)*(2*a+b)) == 0

    # With poles at 0 and infinity, check the coefficient of 1/x in
    # q_b*q_c/x² modulo q_a. The derivative polynomial lacks x².
    def parameters(v):
        return 2*(a+b+c+d-v), -a*b*c*d/v**2

    sa, pa = parameters(a)
    sb, pb = parameters(b)
    sc, pc = parameters(c)
    r = (sa-sb)*(pc-pa)+(sa-sc)*(pb-pa)+(pb-pa)*(pc-pa)*sa/pa
    sigma = a*b+a*c+a*d+b*c+b*d+c*d
    assert sp.cancel(r+2*d*(a-b)*(a-c)*(a+d)/a
                     +2*d*(a-b)*(a-c)*sigma/(b*c)) == 0


def packing(cover):
    return {key: coefficients(getattr(cover, key))
            for key in ("numerator", "half_poles", "lift", "cubic")}


def check_graph(first, second, two, four):
    assert second.cubic.compose(two).rem(first.cubic).is_zero
    for factor, A, B in four:
        assert (first.cubic*B**2-second.cubic.compose(A)).rem(factor).is_zero

        def double(g, coordinate):
            numerator = (coordinate**2*3+g.nth(1))**2
            return (numerator*(g.compose(coordinate)*4).invert(factor)-coordinate*2).rem(factor)

        z = normalize(t, t)
        assert (two.compose(double(first.cubic, z))-double(second.cubic, A)).rem(factor).is_zero


def run():
    check_rank_identities()
    assert recover(normalize(x**6-x-1, x))[2] == []  # S6 has no stable pairing.
    source, center, solutions = recover(normalize(PERMUTED, x))
    assert center is None and source.is_irreducible and len(solutions) == 2
    (first, first_pairing), (second, second_pairing) = solutions
    assert first_pairing == second_pairing and first_pairing.modulus.is_irreducible
    for cover in (first, second):
        certify_identity(source, cover)
        assert cover.cubic.is_irreducible
    two, four = recover_kernel(source, first, second, first_pairing)
    check_graph(first, second, two, four)
    assert len(four) == 1 and four[0][0].degree() == 6
    # A quadratic change of chart alters the equations but preserves the target invariants.
    changed = odd_chart(source, sp.QQ.convert(0), 6, sp.QQ)
    _, _, other_solutions = recover(changed)
    def elliptic(cover):
        return pari.ellinit([pari(str(cover.cubic.nth(k))) for k in (1, 0)])

    assert sorted(str(elliptic(c).j()) for c, _ in solutions) == sorted(
        str(elliptic(c).j()) for c, _ in other_solutions
    )
    curves = [elliptic(first), elliptic(second)]
    assert source.discriminant() % 5 != 0 and all(e.disc() % 5 != 0 for e in curves)
    traces = [int(e.ellap(5)) for e in curves]
    assert traces == [-3, 1]  # Ordinary with distinct geometric endomorphism fields.
    # The old paired chart is an independent source of the pole equations.
    _, _, split_solutions = recover(normalize(SPLIT, x))
    assert sorted(str(c.half_poles.as_expr()) for c, _ in split_solutions) == ["x**2 - 3", "x**2 - 5"]
    return {
        "source": coefficients(source),
        "separation_prime": 5,
        "elliptic_traces": traces,
        "label_cubic": coefficients(first_pairing.modulus),
        "maps": [packing(c) for c, _ in solutions],
        "order_two_x_map": coefficients(two),
        "order_four": [{k: coefficients(v) for k, v in zip(("factor", "x_map", "y_multiplier"), row)}
                       for row in four],
    }


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write", action="store_true")
    args = parser.parse_args()
    receipt = HERE / "galois-receipt.json"
    result = run()
    if args.write:
        receipt.write_text(json.dumps(result, indent=2)+"\n")
    else:
        assert result == json.loads(receipt.read_text())
    print("Both rational maps, the signed 16-point graph, the reciprocal chart and both rank identities passed.")
