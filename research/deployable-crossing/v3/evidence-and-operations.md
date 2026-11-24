# Evidence compatibility and model operation states

This module connects the reference assembly and its continuous-motion receipts to an explicit evidence model and a finite transition graph. It makes two different statements separately: a declared mathematical condition can be established within a supplied model, and a transition can be admitted by the declared simulation. **Neither statement authorizes physical deployment or public opening.** All packaged sensor observations are synthetic fixtures. There are no measured bridge sensor data in this release.

## What runs

From the repository root, with Python 3.11 or later:

```sh
python -S research/deployable-crossing/v3/evidence.py --write
python -S research/deployable-crossing/v3/operation.py --write
python -S research/deployable-crossing/v3/verify_evidence.py
```

Generate the assembly and motion receipts first. The operation program reads `design.json`, `assembly-report.json`, `motion-report.json` and `evidence-report.json` from its own directory. All paths are independent of the current working directory. Python's standard library is sufficient. `verify_evidence.py --unit-only` checks the algorithms and explicit test fixtures without requiring the generated release files.

The APIs are:

```python
solve_evidence(case)
build_evidence_examples(design_sha256)
run_operation(
    evidence_report,
    assembly_receipt=None,
    motion_receipt=None,
    design_sha256=None,                 # keyword argument
    evidence_case_id="independent_channels_synthetic",
    unsafe_graph=False,
    include_controls=True,
)
```

The final five arguments to `run_operation` are keyword-only. Missing assembly or motion receipts produce unresolved guards and block progress beyond the staged/hold alternatives. A missing or malformed design digest raises an input error. The program emits no hardware command.

## Exact compatibility, including uncertainty

A case supplies a finite, nonempty library of candidate states. Every state has a prediction interval for every registered measurement and an explicit Boolean assessment of **every declared hazard obligation**. A candidate may represent a joint combination of engagement, support state and other modeled quantities. It is the same candidate across all readings; different readings cannot silently choose different states.

For state \(s\) and channel \(j\), let its prediction be the closed rational interval \(P_{sj}\). A recorded interval \([a_j,b_j]\) with a bounded nonnegative absolute error \(e_j\) permits

\[
O_j=[a_j-e_j,b_j+e_j].
\]

The compatible states are exactly

\[
C=\{s:P_{sj}\cap O_j\ne\varnothing\text{ for every observed }j\}.
\]

Fractions and finite decimal strings are evaluated as rational numbers; binary floating-point inputs are rejected. Touching endpoints count as compatible. Prediction intervals already describe the allowed model variation. The separately supplied observation error is applied once. Nuisance values inside different channel intervals are independent in this model; correlated nuisance constraints require a richer joint model or a suitably refined state library.

| Outcome | Meaning |
|---|---|
| `condition_established` | Required measurements are present, the compatible set is nonempty, and every compatible candidate discharges every declared hazard obligation. |
| `condition_fails` | Required measurements are present, the compatible set is nonempty, and every compatible candidate fails at least one declared obligation. Different candidates may fail different obligations. |
| `ambiguous` | Compatible candidates disagree, or a required measurement is absent. A missing measurement conservatively prevents either complete protocol decision. |
| `model_mismatch` | No candidate explains all readings, or the input, units, provenance, freshness or expected bindings are invalid. |

In particular, an empty compatible set never proves the condition by vacuous universal quantification. A missing sensor never silently becomes a passing reading. `hazard_witnesses` retains the discharging and failing candidate IDs separately for each obligation. `excluded_states` retains disjoint prediction/observation intervals as exact exclusion witnesses.

When compatible candidates disagree, the program can rank unused measurements. It expands each predicted interval by that channel's declared planning error, counts opposite-decision candidate pairs whose expanded intervals are disjoint, and divides that count by the measurement cost. The result is a finite-library separation criterion. It is not a probability, a guarantee of finding the true physical state, or a proof that adaptive measurement is faster than a fixed protocol.

## Bound sources and explicit missing information

A case contains `model`, `source`, `bindings` and `expected_bindings`. The bindings include the canonical-JSON content digests of the design, model and source record. Canonical JSON uses sorted keys, compact separators and unescaped UTF-8, matching the assembly and motion modules. The source also identifies the model and design to which its observations apply. `solve_evidence` recomputes the model and source content digests, compares them to the expected bindings and checks all cross-references before using any observation. The operation command-line program checks the design digest against the actual design file.

The caller supplies the expected bindings as trust anchors. Matching hashes establish byte/content consistency; they do not authenticate the author, prove calibration or stop a caller from replacing every trust anchor. The private `_seal` helper creates test fixtures. It is not an evidence authentication service.

Every channel declares its unit, maximum allowed error, planning error, maximum age in logical steps and permitted observation provenance. Required observations cannot be satisfied by command-record channels. Unknown measurement IDs, omitted candidate predictions, unknown required channels, negative/excessive errors, nonfinite scalars, reversed bounds, stale or future readings and unit mismatches are rejected. A command record cannot be relabeled as an independent engagement measurement on that channel.

Freshness uses explicitly supplied integer steps. This exercises stale-data behavior without claiming to integrate a real clock, synchronized sensor network or replay-resistant control protocol.

`expected_hazards` must exactly match the model's declared obligations, and each state must explicitly assess all of them. This prevents omissions **within the stated contract**. It does not show that the supplied hazard list covers every physical hazard, or that the finite library contains the actual bridge state. Those remain modeling and empirical obligations.

Provenance distinguishes `synthetic_fixture`, `independent_measurement` and command records. An actual measurement would still require an independently justified source and error budget. Every observation packaged here is labeled synthetic or a synthetic command record; the report records **zero physical measurements**.

## The lock ambiguity examples

The synthetic library is the Cartesian product of three engagement states—full, partial and missing—with two support states—established and missing. Every one of the six states predicts the same actuator command. The command therefore supplies no information about which engagement state occurred.

A separately defined engagement-probe channel distinguishes those three modeled states. A separate support-reaction channel distinguishes the two support states. With the narrow illustrative intervals, both readings together leave one candidate that satisfies the two declared obligations. An engagement reading alone leaves support uncertainty. A broad but still permitted error retains all three engagement states. An impossible reading produces a model mismatch.

The fixture's millimetre and newton values are authored signal examples. They are **not** a calibrated sensor response, a specified connector stroke, a measured support reaction or a hardware dimension derived from the reference assembly. The design digest prevents mixing this study's reports with another assembly; it does not turn the fixture numbers into physical properties. The independent channel is independent of the command in the observation model; no statistical independence or independent failure mode has been demonstrated.

Fourteen recorded cases cover command-only ambiguity, a complete condition, partial/missing engagement, missing support, a missing sensor, broad error, impossible readings, stale readings, wrong units, excessive error, the wrong design, changed model bytes and an omitted hazard. These are regression controls for conditional inference.

## Composition with assembly and motion

The operation checker accepts the shared receipt contract:

- Assembly: a matching `design_sha256`, `physical_approval: false`, and explicit checks for mass conservation, expected model solutions, invalid-mode mechanisms and controlled release.
- Motion: the same `design_sha256`, a matching `geometry_and_mass_source_sha256`, `physical_approval: false`, and explicit checks for continuous support containment, half separation, nominal capture, capture using the loaded structural shape, stage endpoint agreement, release/service support, and service reference clearance.
- Evidence: a matching design digest, complete case count and exact replay of every saved case result. The selected lock case must establish the condition and include both full connector engagement and support-state obligations.

The checker hashes the actual supplied JSON content. It recomputes the hash of the assembly report's whole `assembly` object and requires both reports to reference that same geometry/mass model. It also requires the loaded-capture source hash to match the assembly's actual `capture_reference` object and the motion receipt's `structural_report_sha256` to match the whole assembly report. Its command-line report additionally binds the exact bytes of the design and all three receipt files. Missing Boolean checks remain `unresolved`; the string `"true"` is rejected as a Boolean check. A mismatched design blocks progression. Receipt validation here establishes the linkage and use of recorded conclusions. The assembly and motion generators and their independent verifiers establish their mathematics; this program does not recompute a structural solve behind each guard.

## Directed transition graph

| Transition | Essential model guard |
|---|---|
| Staged → stabilized | Bound reports, mass conservation, solved declared assembly models and continuous support containment. |
| Stabilized → deploying | Established model support state, continuous support containment and matching stage endpoints. |
| Deploying → captured | Continuous support containment, half separation, and both nominal and structurally loaded capture within the declared envelope. |
| Captured → locked | The selected evidence condition establishes both declared lock/support obligations; support remains established. |
| Locked → transferred | Checked controlled release, support containment during release, lock evidence, established support and the discrete secured-lock state. This transition releases the modeled cantilever restraint. |
| Transferred → service | The model supports and locks remain established, service support and reference clearance are satisfied, transfer has occurred, the declared assembly cases solve and the endpoints match. |
| Pre-service mode or recovery → hold | An accepted model hold is available; no hardware actuation is implied. Model service is an endpoint of this bounded deployment graph. |
| Hold → recovery → staged | A recovery-path certificate is required and is currently unresolved. |

The graph retains separate flags for established support, secured locks and released cantilever restraint. It enumerates every state reachable under the fixed guard valuation and checks that release never precedes both locks and support. Model service additionally requires a completed transfer. Concrete transition-ID witness paths are retained for reachable modes and invariant violations.

The report uses one synthetic evidence case as a guard fixture; it is not a time series of observations acquired during deployment. A physical time-varying implementation would need to reacquire and revalidate evidence at the relevant transitions and model changes. That implementation is absent.

The continuous side consists of the supplied whole-path receipt conclusions. The discrete side checks the declared guarded transitions. This is a bounded hybrid-model composition, **not** a proof that every real motion, contact, component failure, observation delay or emergency action is represented by the graph.

## Safety, hold and liveness are different questions

The normal graph can satisfy its release invariant while failing to reach model service. An accepted hold is an explicit stopping state. It is not counted as an unexpected deadlock and is not counted as service. An unexpected terminal state in another mode is reported with a witness path.

The controls deliberately add an unguarded staged-to-transferred release. Reachability produces a concrete violating path with no locks or supports established. Removing a required support certificate blocks progress; removing the accepted hold escape as well produces an unexpected deadlock. A failed controlled-release certificate blocks transfer even if the preceding model stages can be reached.

`model_service_state_reached` is existential reachability within this graph. The checker does not claim that every execution reaches service. Hold is permitted, and recovery is unresolved.

`public_opening.allowed` is **always false** in this research module, including a case where model service is reachable. Site permission, structural capacity/foundations, an accessible route, calibrated field evidence, hardware/actuation details, work-zone operation and inspection/recovery procedures remain unresolved gates. This list is a study interface, not a claim of exhaustive regulatory or hazard coverage.

## Verification and relationship to earlier applications

The standard-library verifier includes an independently enumerated integer-grid compatibility oracle, exact touching-endpoint and broad-error controls, malformed-input cases, source/model/design tampering, evidence-result replay, missing/incorrect receipt checks, an unsafe-release witness and a deadlock-versus-hold distinction. The final checks rebuild the saved evidence and operation reports and check source-file hashes.

The evidence design follows the existing [wafer application's](https://github.com/towre676-cloud/crest-root-calculus/blob/main/applications/wafer/README.md) separation of candidate compatibility from a decision and the [water application's](https://github.com/towre676-cloud/crest-root-calculus/blob/main/applications/water/README.md) warning that an overly narrow model can exclude the true state. The wafer replay's measured dies did not fit its illustrative tolerance; the calibrated water replay did not demonstrate an adaptive advantage over its fixed comparison. Those are reasons to retain mismatch and calibration as first-class outcomes here. Neither dataset validates a bridge sensor.

The reusable result is an executable, uncertainty-aware contract between a declared model, its evidence and its transition checks. Physical hardware design, calibration, validation and authorization are still separate work.

## Separate outcome accounting from model feasibility

`outcomes.py` connects the study to the recovered Safe Roads `local_crossing_harm` model. The inspected source is `gap_closure/model_repair/models.py` in *Safe_Roads_Recovery_2026-11-23*, revision 2, with SHA-256 `d2c6a478d8276f376367ae80dfff1e439a9a214844e22d93ce3413979283412f`. The local adapter preserves its formula and scope, uses exact rational arithmetic, and adds explicit unknown inputs and design/site bindings. It does not import an empirical benefit estimate.

For **one site, one declared period and one outcome**, the model computes the range enclosure

\[
H_0=N_0r_0,\qquad
H_1=N_1\bigl[au\,r_b+(1-au)r_f\bigr]+D+S,\qquad
\Delta H=H_0-H_1.
\]

Here \(N_0\) and \(N_1\) are baseline and intervention person-crossing exposures; they are not assumed equal. Availability \(a\) is the fraction of intervention crossing opportunities for which the facility is available. It is not automatically equal to clock-time uptime. Uptake \(u\) is conditional on that availability. The baseline, facility and fallback event risks are \(r_0,r_b,r_f\). Deployment harm \(D\) and displaced-site harm \(S\) must have explicit bounds and a supplied explanation establishing their nonoverlap with each other and main-site events. The existence or truth of that explanation still needs evidence.

Deaths and serious nonfatal injuries are evaluated separately. All nine inputs carry units, site, design, outcome, period and provenance. The interval calculation bounds an independently variable input box; it does not estimate causality or attach statistical confidence. It cannot establish that a bridge helps by itself. Extra harm can reverse the sign of an otherwise positive synthetic calculation.

The actual report reads the six site IDs and labels from the current site report and binds its complete canonical content hash. It records twelve site/outcome evaluations. **Every one remains unknown:** all site-specific exposure, availability, uptake, risk, deployment and displacement quantities are null, as are the target period and extra-harm boundary evidence. Unknown values are never replaced by zero. No geometry result, nominal operation path or selected scenario supplies these missing inputs.

Three clearly separated synthetic controls exercise the accounting: two exact point calculations, one for each outcome, and an interval example whose net range includes zero. Their illustrative values are unrelated to Tillamook measurements and are not predictions of deaths or injuries avoided.

```sh
python -S research/deployable-crossing/v3/outcomes.py --write
python -S research/deployable-crossing/v3/verify_outcomes.py
```

The verifier enumerates all 512 corners of each of the three nine-input boxes, including duplicate vertices for point intervals, and checks malformed units, mixed sites/outcomes/periods, missing inputs, sign reversal and exact report rebuilding. An optional `--source-model /path/to/recovered/models.py` checks the original file's pinned hash and compares all three controls directly with its `local_crossing_harm` implementation. That optional comparison passed when the adapter was created; normal offline replay does not require the external recovery archive.
