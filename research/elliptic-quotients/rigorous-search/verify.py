"""Reconstruct the exact lattices and compare with a certified period shortlist."""
import argparse
import json
from pathlib import Path
from collections import Counter
from time import monotonic


from curves import CASES, qm_tangents, split_case
from exact_lattices import (I, check_end_shortlist, check_qm, check_stored_ns, generated_order, homology_action,
                            lattices, norm_endomorphism, plane_set, rosati_count)
from isotropic import spectrum_from_lattice
from periods_arb import period_matrix_arb, set_prec
from search import elliptic_shortlist, frame


def replay(name, precision, shortlist=True):
    degree, trace = CASES[name]
    set_prec(precision)
    if name == "QM":
        f, tangents = qm_tangents()
        rank = 4
    else:
        f, maps, rank, stored = split_case(name)
    data = period_matrix_arb(f, bits=precision)
    J, G = frame(data)
    if name == "QM":
        qm = [[homology_action(data, M) for M in embedding] for embedding in tangents]
        anchors = [I, *qm[0]]
    else:
        anchors = [I, *(norm_endomorphism(data, G, a, b, d) for a, b, d in maps)]
    span = generated_order(anchors, rank)
    end, ns, gram, theta, trace_gram = lattices(data["E"], span)
    spec, classes = spectrum_from_lattice(gram, theta, degree)
    exact = plane_set(data["E"], ns, classes)
    receipt = {"equation": f, "end_rank": rank, "ns_rank": len(ns), "ns_discriminant": int(gram.det()),
               "degree_bound": degree, "spectrum": spec, "trace_bound": trace,
               "endomorphisms_up_to_sign": rosati_count(trace_gram, trace)}
    if name == "QM":
        moduli = {check_qm(end, embedding, gram) for embedding in qm}
        if exact or gram.det() != 12:
            raise ArithmeticError("Unexpected QM lattice invariants.")
        receipt.update(reduced_discriminant=6, anisotropy_modulus=moduli.pop())
    else:
        expected = {int(k): v for k, v in stored["spectrum_through_X"].items() if int(k) <= degree}
        if spec != expected or gram.det() != stored["NS_discriminant"]:
            raise ArithmeticError("The exact homology lattice disagrees with the torsion certificate.")
        check_stored_ns(anchors, maps, ns, stored)
    if shortlist:
        possible, count = elliptic_shortlist(data, J, G, degree)
        if possible != exact or Counter(exact.values()) != spec:
            raise ArithmeticError(f"Period shortlist differs from exact planes: {set(possible)^set(exact)}")
        receipt["witness_vectors_up_to_sign"] = count
        receipt["matched_planes"] = len(exact)
        accepted, rejected = check_end_shortlist(G, J, end, trace_gram, trace)
        if accepted != receipt["endomorphisms_up_to_sign"]:
            raise ArithmeticError("Ambient enumeration differs from the exact endomorphism count.")
        receipt["ambient_commutator_rejections"] = rejected
        if name == "C":
            # Lower the commutator penalty to exercise the exclusion path.
            accepted, rejected = check_end_shortlist(G, J, end, trace_gram, trace, penalty=1)
            if accepted != receipt["endomorphisms_up_to_sign"] or rejected == 0:
                raise ArithmeticError("The low-penalty exclusion control failed.")
            receipt["low_penalty_commutator_rejections"] = rejected
    return receipt


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("cases", nargs="*", metavar="CASE")
    parser.add_argument("--precision", type=int, default=300)
    parser.add_argument("--lattices-only", action="store_true")
    parser.add_argument("--output", type=Path)
    args = parser.parse_args()
    if any(name not in CASES for name in args.cases):
        parser.error("Unknown case. Choose from " + ", ".join(CASES))
    results = {}
    for name in args.cases or CASES:
        start = monotonic()
        results[name] = replay(name, args.precision, not args.lattices_only)
        r = results[name]
        print(f"{name}: End rank {r['end_rank']}, NS discriminant {r['ns_discriminant']}, "
              f"{sum(r['spectrum'].values())} planes through degree {r['degree_bound']} "
              f"({monotonic()-start:.1f}s)", flush=True)
    if args.output:
        args.output.write_text(json.dumps(results, indent=2) + "\n", encoding="utf-8")
    elif not args.cases and not args.lattices_only and args.precision == 300:
        expected = json.loads(Path(__file__).with_name("results.json").read_text(encoding="utf-8"))
        if json.loads(json.dumps(results)) != expected:
            raise ArithmeticError("The complete replay disagrees with the stored results.")


if __name__ == "__main__":
    main()
