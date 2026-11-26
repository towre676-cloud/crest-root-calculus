"""Replay algebraic Richelot equations and their exact Humbert forms."""

import argparse
from itertools import combinations
import json

from flint import acb, acb_mat, acb_poly, arb, ctx

from verify_d163 import (
    D, HERE, PARI, compute, elliptic_data, embeddings, even_sextic, read, require,
    targets, torsion_data,
)
from mestre_link import PA, certify_cm_order, humbert_from_order
from certify_link import gram, isometry
from verify_d43 import pair_spectrum

PRECISION = 2200
SELECTOR_SCALE = 2**100
PATHS = HERE / "richelot_models_-163.json"
RESULT = HERE / "richelot_certified_-163.json"


def evaluate(expression, value):
    out = acb(0)
    for a in PARI.Vec(PARI(expression)):
        out = out*value + acb(arb(int(a.numerator()))/int(a.denominator()))
    return out


def select(roots, rectangle):
    x, y = rectangle
    region = acb(arb(x, 2)/SELECTOR_SCALE, arb(y, 2)/SELECTOR_SCALE)
    hits = [v for v in roots if region.contains(v)]
    require(len(hits) == 1 and all(not region.overlaps(v) for v in roots if v is not hits[0]),
            "The selector does not identify a unique algebraic root.")
    return hits[0]


def distinct(roots):
    return all(r.is_finite() for r in roots) and all(not r.overlaps(s) for r,s in combinations(roots, 2))


def chart_roots(roots, chart):
    """Send three roots to 0, 1, infinity, then apply t -> 1/(t-2-3i)."""
    require(len(chart) == len(set(chart)) == 3 and all(0 <= j < 6 for j in chart),
            "Invalid Möbius chart.")
    a, b, c = (roots[i] for i in chart)
    remaining = [(z-a)*(b-c)/((z-c)*(b-a)) for i,z in enumerate(roots) if i not in chart]
    z0 = acb(2, 3)
    result = [1/(z-z0) for z in [acb(0), acb(1)] + remaining] + [acb(0)]
    require(distinct(result), "The chart does not isolate six distinct roots.")
    return result


def radical(value):
    return (acb(0, 1)*(-value).sqrt(analytic=True)
            if value.real < 0 else value.sqrt(analytic=True))


def richelot_roots(roots, pairing):
    """Return the roots of the product of the three Richelot brackets."""
    require(len(pairing) == 3 and all(len(pair) == 2 for pair in pairing) and
            sorted(i for pair in pairing for i in pair) == list(range(6)), "Invalid Richelot matching.")
    factors = [acb_poly([-roots[i], 1])*acb_poly([-roots[j], 1]) for i,j in pairing]
    determinant = acb_mat([[f[i] for i in range(3)] for f in factors]).det()
    require(not determinant.contains(0), "The Richelot determinant is not certified nonzero.")
    result = []
    for a,b in zip(factors, factors[1:]+factors[:1], strict=True):
        bracket = a.derivative()*b-a*b.derivative()
        c,b,a = [bracket[i] for i in range(3)]
        require(not a.contains(0), "The Richelot bracket is not certified quadratic.")
        delta = b*b-4*a*c
        require(not delta.contains(0), "The Richelot bracket may have a repeated root.")
        square_root = radical(delta)
        require(square_root.is_finite(), "The quadratic radical could not be isolated.")
        result.extend([(-b-square_root)/(2*a), (-b+square_root)/(2*a)])
    require(distinct(result), "The Richelot sextic is not certified squarefree.")
    return result


def period_data(roots, bits=128):
    """Use the existing integration proof directly on enclosed branch points."""
    with ctx.workprec(bits):
        options = []
        for zr,zi in PA.Z0_CHOICES:
            z = acb(arb(zr.numerator)/zr.denominator, arb(zi.numerator)/zi.denominator)
            if any((a-z).contains(0) for a in roots):
                continue
            moved = [1/(a-z) for a in roots]
            midpoints = list(map(PA.mid, moved))
            if len(set(midpoints)) != 6 or not distinct(moved):
                continue
            span = max(abs(a-b) for a,b in combinations(midpoints,2))
            score = max(PA.clearance(midpoints, i) for i in range(6))/span
            options.append((score, moved))
        for _,moved in sorted(options, key=lambda option: option[0], reverse=True):
            try:
                data = PA.star_periods_from_roots(moved, acb(1))
                P = list(map(list, zip(*data["periods"][:4])))
                E = [row[:4] for row in data["intersection"][:4]]
                period, form = acb_mat(P), acb_mat(E)
                relation = period*form.inv()*period.transpose()
                require(all(relation[i,j].contains(0) for i in range(2) for j in range(2)),
                        "Riemann relation failed.")
                return {"Pi": P, "E": E}
            except (ArithmeticError, ValueError):
                continue
    raise ArithmeticError("No certified period paths for the stated chart.")


def verify_paths(paths):
    initial = read(HERE / "models_certified_-163.json")
    require(initial == compute(), "The exact anchor models differ from their certificate.")
    archive = targets()
    computed = {name: archive[name] for name in initial["models"]}
    all_roots = {}
    depth = {}
    results = {}
    with ctx.workprec(PRECISION):
        embedding = select(embeddings(PARI(initial["field"])), paths["field_root"])
        for name, record in initial["models"].items():
            rhs = acb_poly([evaluate(expr, embedding) for expr in record["model_constant_first"]])
            roots = rhs.roots(tol=arb(2)**-440, maxprec=4*PRECISION)
            require(len(roots) == 6 and distinct(roots), "The anchor roots do not isolate.")
            all_roots[name] = [select(roots, box) for box in paths["anchor_roots"][name]]
            require(len(all_roots[name]) == 6 and distinct(all_roots[name]), "An anchor selector was repeated.")
            depth[name] = 0
        field, torsion, kernels = torsion_data(read(HERE / "d163_torsion_field.json"))
        require(field == PARI(initial["field"]), "The auxiliary gluing uses a different field.")
        _, elliptic = elliptic_data(torsion, kernels)
        for name, seed in paths["seeds"].items():
            i, j = seed["indices"]
            perm = seed["permutation"]
            require(name not in depth and 0 <= i < j < len(elliptic) and sorted(perm) == [0, 1, 2],
                    "The auxiliary gluing record is invalid.")
            factors = even_sextic(elliptic[i][0], [elliptic[j][0][k] for k in perm])
            branch = []
            for factor in factors:
                value = evaluate(str(PARI.lift(-factor.polcoef(0, "x"))), embedding)
                root = radical(value)
                branch.extend([-root, root])
            require(distinct(branch), "The auxiliary gluing is not a smooth sextic.")
            all_roots[name] = branch
            depth[name] = 0

        for edge in paths["edges"]:
            name, parent = edge["name"], edge["parent"]
            require(parent in depth and name not in depth, "Invalid path order.")
            depth[name] = depth[parent] + 1
            model_roots = richelot_roots(chart_roots(all_roots[parent], edge["chart"]), edge["pairs"])
            all_roots[name] = model_roots
            data = period_data(chart_roots(model_roots, edge["period_chart"]))
            multiplier = 6 * 2**depth[name]
            scalar, order = certify_cm_order(data, multiplier, D)
            q = gram(humbert_from_order(data, order))
            matches = [label for label, target in archive.items()
                       if q.matdet() == target.matdet() and PARI.qfisom(q, target)]
            expected = [name] if name in archive else []
            require(matches == expected, f"Unexpected isometry classes for {name}: {matches}.")
            if name in archive:
                computed[name] = q
                results[name] = {
                    "humbert_form": [[int(q[i,j]) for j in range(3)] for i in range(3)],
                    "isometry": isometry(q, archive[name]),
                    "scalar_cm_matrix": [[int(scalar[i,j]) for j in range(4)] for i in range(4)],
                    "richelot_steps": depth[name],
                }
                print(f"{name}: verified exact model after {depth[name]} steps.", flush=True)
    require(set(computed) == set(archive), "The path certificate does not cover all 38 twin lattices.")
    spectra = {}
    for pair in range(19):
        names = tuple(f"p{pair}_q{j}" for j in (1,2))
        if all(name in computed for name in names):
            spectra[f"p{pair}"] = pair_spectrum({name: computed[name] for name in names}, names, discriminant=D)
    return {"models": results, "spectra": spectra}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args()
    results = verify_paths(read(PATHS))
    if args.write_certificate:
        RESULT.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    else:
        require(results == read(RESULT), "The exact path certificate differs from the replay.")
    print(f"Verified {len(results['models']) + 2} equations and {len(results['spectra'])} complete twin pairs.")
