"""Certify the seven discriminant -67 twin pairs from their quartic equations."""

import argparse
from collections import Counter
import json

from flint import acb, acb_poly, arb, ctx

from certify_link import HERE, gram, isometry, read, require
from mestre_link import certify_cm_order, humbert_from_order, normalised_quintic, period_matrix_acb
from verify_bielliptic import build, certify
from verify_d43 import pair_spectrum
from verify_d163_paths import chart_roots, distinct, period_data
from verify_mestre_descent_d43 import replay as verify_models
from verify_richelot_d43 import PARI, verify_paths


def humbert(model, field, embedding, steps):
    f = normalised_quintic(model["model_constant_first"], field, embedding, 2100)
    try:
        data = period_matrix_acb(f, 400)
    except (ArithmeticError, ValueError, ZeroDivisionError):
        with ctx.workprec(1800):
            # The quintic has five finite roots. Move the root at infinity
            # into the affine chart before using the branch-point integrator.
            finite = acb_poly(f).roots(tol=arb(2)**-500, maxprec=6000)
            require(len(finite) == 5 and distinct(finite), "The quintic roots do not isolate.")
            move = 2+3j
            require(all(not (v-move).contains(0) for v in finite), "The moved poles meet a branch point.")
            branch = [1/(v-move) for v in finite] + [acb(0)]
            data = period_data(chart_roots(branch, (0, 1, 2)), 300)
    scalar, order = certify_cm_order(data, 6 * 2**steps, -67)
    return humbert_from_order(data, order), [[int(scalar[i, j]) for j in range(4)] for i in range(4)]


def replay():
    anchors = read(HERE / "certified_bielliptic_-67.json")
    field, classes, models = build(-67)
    require(certify(-67, field, classes, models) == anchors, "The bielliptic receipt differs from the replay.")
    verify_models(-67)
    distances = verify_paths(-67)
    archive = read(HERE / "twin_lattices_-67.json")
    sextics = read(HERE / "mestre_quartic_models_-67.json")
    require(len(archive) == 14 and set(archive) == set(sextics) | set(models), "The archive is incomplete.")
    rows = {}
    forms = {name: gram(archive[name]["q"]) for name in models}
    for name, record in sorted(sextics.items()):
        row = {"richelot_steps": distances[name], "embeddings": []}
        for embedding in range(4):
            q, scalar = humbert(record, field.replace("x", "w"), embedding, distances[name])
            form = gram(q)
            matches = [label for label, item in sorted(archive.items()) if PARI.qfisom(gram(item["q"]), form)]
            require(matches == [name], f"{name} at embedding {embedding} matches {matches}.")
            row["embeddings"].append({"humbert_form": q, "scalar_cm_matrix": scalar,
                                      "isometry": isometry(form, gram(archive[name]["q"]))})
            forms[name] = form
        rows[name] = row
        print(f"{name}: exact Humbert form certified at all four embeddings.", flush=True)
    spectra = {}
    for i in range(7):
        pair = f"p{i}"
        names = (f"{pair}_q1", f"{pair}_q2")
        spectra[pair] = pair_spectrum({name: forms[name] for name in names}, names, discriminant=-67)
    return {"models": rows, "spectra": spectra}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--write-certificate", action="store_true")
    args = parser.parse_args()
    receipt = replay()
    path = HERE / "mestre_link_quartic_-67.json"
    if args.write_certificate:
        path.write_text(json.dumps(receipt, indent=2) + "\n", encoding="utf-8")
    else:
        require(receipt == read(path), "The lattice receipt differs from the replay.")
    counts = Counter(row["first_degree"] for row in receipt["spectra"].values())
    print(f"Seven certified pairs. First-difference degrees: {dict(sorted(counts.items()))}.")
