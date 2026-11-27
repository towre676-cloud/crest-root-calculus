"""Equation-only recovery and independent arithmetic checks on Howe's examples."""
from __future__ import annotations

import argparse
import json
from importlib.util import module_from_spec, spec_from_file_location
import sys
from math import lcm, prod
from pathlib import Path

import sympy as sp

from reconstruct import (
    coefficients, cover_identity, discover, pari, pari_polynomial, recover_graph, serialise, x,
)

EXAMPLES = (
    ("C70", 2*x**3-3*x*x-41*x+110, x**3-51*x*x+425*x+179, (7, 10)),
    ("C28,1", x*x+x, x**6+3*x**5+5*x**4-4*x*x-10*x+4, (4, 7)),
    ("C28,2", x*x+x, x**6+3*x**5+3*x**4+13*x**3-6*x*x+18*x, (4, 7)),
    ("C28,3", x*x+x, 4*x**6-2*x**5+18*x**4+3*x**3+13*x*x+23*x-11, (4, 7)),
    ("C28,4", x*x+x, 28320768*x**6+167100960*x**5+213557586*x**4-35302844*x**3
     +154134546*x*x-155174208*x+40064896, (4, 7)),
    ("C28,5", x*x+x, 25*x**6-455*x**5+1675*x**4+2494*x**3+570*x*x-1210*x, (4, 7)),
    ("C27 split", 2*x**3+3*x*x-3*x+2, 6*x**3+6, (3, 9)),
    ("C39", 6*x**3+6*x*x-7*x-9, -2*x-2, ()),
    ("C27 no cubic factor", 6*x**3+9*x*x+6*x-1, -3*x*x, ()),
)


def source_model(h, f):
    return sp.Poly(h*h+4*f, x, domain=sp.QQ)


def elliptic(cover):
    return pari.ellinit([pari(str(cover.cubic.nth(k))) for k in (1, 0)])


def genus_two_counts(source, primes):
    model = pari_polynomial(source)
    disc = int(source.discriminant())
    counts = {}
    for p in primes:
        assert p > 2 and disc % p and source.LC() % p
        counts[p] = int(pari.hyperellcharpoly(model * pari.Mod(1, p)).subst("x", 1))
    return counts


def torsion_bound(counts):
    """Intersect prime-to-residue-characteristic bounds one prime at a time."""
    orders = {p: sp.factorint(n) for p, n in counts.items()}
    primes = set().union(*orders.values())
    return prod(ell**min(f.get(ell, 0) for p, f in orders.items() if p != ell) for ell in primes)


def image(cover, point):
    coordinate, ordinate = point
    polys = (cover.numerator, cover.poles, cover.lift)
    n, d, h = [p.nth(3) if coordinate is None else p.eval(coordinate) for p in polys]
    return [pari(str(n/d)), pari(str(h*ordinate/d**2))] if d else [0]


def elliptic_order(target, point):
    assert target.ellisoncurve(point)
    return int(target.ellorder(point))


def certify_howe_divisor(source, h, covers, targets):
    points = [(a, 2*b+h.subs(x, a)) for a, b in ((2, 17), (4, 23))]
    assert all(source.eval(a) == b*b for a, b in points)
    counts = genus_two_counts(source, (17, 19, 23))
    assert counts == {17: 280, 19: 280, 23: 630}
    assert torsion_bound(counts) == 70
    orders = []
    for cover, target in zip(covers, targets):
        leading = sp.sqrt(source.LC())
        plus = image(cover, (None, leading))
        minus = image(cover, (None, -leading))
        assert target.elladd(plus, minus) == pari([0])
        orders.append(elliptic_order(target, target.elladd(*(image(cover, p) for p in points))))
    # The pushforward kernel is killed by 3, so [210]D = 0 before using reductions.
    assert sorted(orders) == [7, 10]
    return orders, counts


def certify_order_27(source, covers, targets):
    point, opposite = (-1, 6), (-1, -6)
    assert source.eval(-1) == 36
    assert [elliptic_order(target, target.ellsub(image(cover, point), image(cover, opposite)))
            for cover, target in zip(covers, targets)] == [3, 9]
    # Thus [27]([point]-[opposite]) = 0 over Q, by the dual kernel.
    counts = genus_two_counts(source, (5, 7))
    assert torsion_bound(counts) == 27
    sys.path.insert(0, str(Path(__file__).resolve().parents[2] / "torsion-denominators"))
    from jac import Jacobian, odd_model

    p = 7
    coeffs, root = odd_model([int(c) for c in reversed(source.all_coeffs())], p)
    jac = Jacobian(coeffs, p, 1)
    u = jac.F(point[0]-root)**-1
    divisor = (jac.x-u, jac.R(u**3*point[1]))
    difference = jac.mul(2, divisor)
    assert not jac.is_zero(jac.mul(9, difference))
    assert jac.is_zero(jac.mul(27, difference))
    return counts


def certify_twists(source, covers, targets):
    p, q = 17, 23
    assert sp.legendre_symbol(5, p) == sp.legendre_symbol(5, q) == -1
    counts = genus_two_counts(source * 5, (p, q))
    assert counts == {17: 352, 23: 486}
    assert torsion_bound(counts) == 2
    # One target has a rational 2-torsion point. Twisting preserves it, and the
    # degree-9 quotient of the product preserves its order.
    assert any(int(e.elltors()[0]) % 2 == 0 for e in targets)
    traces = [int(e.ellap(p)) for e in targets]
    assert sorted(traces) == [-2, 4]
    assert [a*a-4*p for a in traces] in ([-64, -52], [-52, -64])
    # The distinct ordinary reduction fields exclude even geometric isogeny.
    for cover in covers:
        twisted = type(cover)(cover.numerator*5, cover.poles, cover.lift*5,
                              cover.cubic.compose(sp.Poly(cover.cubic.gen/5, cover.cubic.gen))*125)
        assert cover_identity(source*5, twisted)
    return counts


def local_degrees(source):
    path = Path(__file__).resolve().parents[2] / "torsion-denominators"
    sys.path.insert(0, str(path))
    from certify import certify_curve

    spec = spec_from_file_location("verify_torsion_denominators", path / "verify.py")
    verifier = module_from_spec(spec)
    spec.loader.exec_module(verifier)
    coefficients = [int(c) for c in reversed(source.all_coeffs())]
    for p in (17, 23):
        certificate, _ = certify_curve(coefficients, p)
        certificate["curve"] = coefficients
        assert verifier.check_certificate(certificate) == [3, 3]
        print(f"Independent local degree certificate at {p}: 3, 3")


def run(check_local=False):
    records = {}
    for name, h, f, expected in EXAMPLES:
        source = source_model(h, f)
        covers = discover(source)
        assert len(covers) == len(expected)
        record = records[name] = {"sextic": coefficients(source), "maps": [], "kernel": []}
        if not covers:
            print(f"{name}: no rational degree-three map")
            continue
        # Search and kernel reconstruction receive only the source and derived maps.
        components = recover_graph(source, *covers)
        targets = [elliptic(cover) for cover in covers]
        assert targets[0].j() != targets[1].j()
        orders = [int(e.elltors()[0]) for e in targets]
        assert sorted(orders) == list(expected)
        for target, order in zip(targets, orders):
            generators = target.elltors()[2]
            assert lcm(*(elliptic_order(target, point) for point in generators)) == order
        record.update(serialise(covers, components))
        record["elliptic_torsion_orders"] = orders
        if name == "C70":
            assert int(targets[0].ellglobalred()[0]) == 66 and int(targets[1].ellglobalred()[0]) == 858
            record["divisor_image_orders"], record["jacobian_counts"] = certify_howe_divisor(source, h, covers, targets)
            record["twist_by_5_jacobian_counts"] = certify_twists(source, covers, targets)
            record["jacobian_torsion_order"] = 70
            if check_local:
                local_degrees(source)
        elif name == "C27 split":
            record["jacobian_counts"] = certify_order_27(source, covers, targets)
            record["jacobian_torsion_order"] = 27
        else:
            # Find upper bounds; even when a bound is loose the coprime-to-3
            # factor generators give a rational point of order 28.
            model = pari_polynomial(source)
            good = [p for p in (5, 7, 11, 13, 17, 19, 23, 29, 31, 37)
                    if int(model.poldisc()) % p and source.LC() % p]
            counts = genus_two_counts(source, good)
            bound = torsion_bound(counts)
            # The degree-nine isogeny identifies prime-to-three rational torsion.
            assert bound % 3 != 0 and bound % 28 == 0, (name, counts, bound)
            record["jacobian_counts"] = counts
            record["jacobian_torsion_order"] = 28
        print(f"{name}: two exact covers, full 3-torsion graph, elliptic orders {orders}, "
              f"Jacobian torsion C{record['jacobian_torsion_order']}")
    return records


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--receipt", type=Path, help="Write the verified maps and graph over Q.")
    parser.add_argument("--local", action="store_true", help="Also run the existing local degree certifier.")
    args = parser.parse_args()
    records = run(args.local)
    if args.receipt:
        args.receipt.write_text(json.dumps(records, indent=2) + "\n")
