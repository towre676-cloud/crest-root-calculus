# Cartoner stroke trial

[Open the deployed application](https://towre676-cloud.github.io/crest-root-calculus/cartoner/) · [Download the offline trial](../../downloads/linkwright-cartoner-trial.zip) · [Trial guide](https://towre676-cloud.github.io/crest-root-calculus/cartoner/guide.html)

Run the extracted `index.html` in a current browser. The program reads local CSV and JSON, checks full-turn or limited-angle four-bar travel, checks analog process requirements across continuous encoder-phase intervals, and checks event separation. Every cycle retains its own evidence. The exported JSON contains the exact inputs, results, version and optional paired inspection notes. Reopening the record recomputes the mathematical results. Record replay checks consistency. It supplies no digital signature or independent validation of the original measurements.

## Research used in the trial

The source review used the integrated monograph through §148, the research register, the current application studies and the research branches present at the review. The release extends three results to cartoner evaluation.

| Existing work | Engineering use |
| --- | --- |
| [Four-bar configuration curves](../../research/configuration-and-certificates/mechanisms/README.md) and [Linkwright engineering verification](../verification/README.md) | The extrema of the triangle inequalities decide strict full-turn closure. The trial accepts independently bounded lengths, returns an attaining build, and extends the rectangle calculation to a specified operating cosine range. A geometry violation is a counterexample to the guarantee across permitted builds. |
| [Crossing v3 continuous motion](../../research/deployable-crossing/v3/README.md) and [compatible evidence](../../research/deployable-crossing/v3/evidence-and-operations.md) | Carry the entire declared interval through the calculation, bind evidence to the revision, and retain ambiguity and contradictory evidence. The new phase-trace solver intersects all reading-error bands with a slope constraint to enclose unsampled motion. |
| [Measurement applications](../../applications/README.md) and [solver comparisons](../../research/solver-comparison/README.md) | Retain all compatible outcomes, compare against the existing inspection on the same inputs, and distinguish reference-case replay from held-out field performance. |

The [coupled six-bar study](../../research/configuration-and-certificates/sixbar/README.md) distinguishes loop arrangements with the same input angle. Each configured loop in this trial must have a specified drive range. Shared links, tool bodies and contact geometry need a validated task-specific model. The later quotient, finite-kernel, gluing, torsion and spectral results concern exact algebraic invariants. Using them to infer carton timing, strength, or contact from an encoder trace would require an additional physical model and a proved connection.

## Acceptance checks

The synthetic opener has a sampled first-cycle margin of +2.5° and a compatible minimum in [−4°, +3°]. Its second cycle gives [+0.5°, +7.5°]. The synthetic folder yields a return-clearance minimum in [−1.35, −0.25] mm on cycle 1. The full-turn regression gives an attained −0.00006 mm reach margin at 180.25°. These outcomes follow from the supplied fixture inputs.

The offline trial guide specifies units, CSV formats, uncertainty, bounds, report interpretation, the factory comparator and physical follow-up. Machine-specific proof of clearance requires calibrated clearance data or a separately verified swept-body model. Further empirical work includes obtaining independent OEM drawings, traces, fault labels and a buyer-defined acceptance decision.

With Node.js 18 or newer, run:

```sh
node --test test.cjs
node cli.cjs examples/opener.json examples/opener.csv > record.json
node cli.cjs --replay record.json
```

The independent numerical control evaluates every observation cone at all candidate phases on 240 seeded cases. Boundary, inconsistent-data, equality, incomplete coverage, handoff, restricted-travel and replay tests run alongside it. Use `python release.py` from the repository checkout to recreate the deterministic ZIP and embedded example data.
