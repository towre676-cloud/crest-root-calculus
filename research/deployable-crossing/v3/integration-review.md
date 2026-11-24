# One reference bridge, connected evidence

This version joins the structural assembly, prescribed deployment, Tillamook site screen, evidence inference and operating-state checks around one design identity. The [offline workbench](workbench.html) makes those connections inspectable. The [integrated report](integrated-report.json) retains the component results and their bindings.

**Research scope:** this is a shared reference assembly and an executable model of its interfaces. Dimensions, materials, mass allowances, loads and motion are declared research inputs. Mapped Tillamook features are source data with the limits recorded in v2. A software pass is not a capacity rating, site selection, actuator design or public-use authorization.

## What the integrated run establishes

The common design identity is `c82b8359a7e5fcb07c52bfb33628b5fd386c07b23b4d37f9b43b25fbc221358f`. Every component report refers to it, while separate content hashes bind the actual assembly, structural capture state, continuous motion, site results and evidence record.

| Interface result | Computed finding | Boundary |
| --- | --- | --- |
| Structural assembly | 40 physical nodes, 114 bars and 10,206.794 kg of declared moving bridge mass; six stable linear cases and two deliberate mechanism cases | Response and rank, not selected capacities; fixed bases are additional mass |
| Continuous prescribed deployment | 24 support paths and 68 conservative swept boxes for the three motion stages | Declared force box and rigid reference motion; omitted hardware and actual engagement remain unresolved |
| Complete structural base wrenches | 48 combinations across six stable structural cases, two bases and four fixture-force corners have a resultant inside the pad rectangle | Existence of nonnegative coplanar vertical pad reactions; no bearing, sliding, torsional or hardware capacity |
| Controlled support transfer | Eight continuous release-support paths, using the full FE force and moment resultants | Commanded quasi-static ideal restraint motion; not instantaneous release or a physical controller |
| Structurally loaded capture | Maximum absolute coordinate mismatch is approximately 1.029 mm within the declared ±20 mm capture box | Independent-half elastic dead load only; foundation error, fabrication, contact and lock geometry are absent |
| Reference-plane clearance | Minimum declared loaded member surface is approximately 5.46567 m, leaving 82.47 mm above the reference plus overlay | Flat model road plane and small-displacement axial members, not surveyed or dynamic clearance |
| Site study | 216 access scenarios at 18 placements; two have neither tested direct building/vehicle-centerline intersections nor possible swept-bound building intersections | Both remain research hypotheses; zero approved sites |
| Evidence and operation | 14 synthetic evidence cases; guarded model service is reachable and physical opening remains disabled | No field observations and no physical-controller completeness claim |
| Local outcomes | Twelve site/outcome records remain unknown across six sites and two disjoint outcomes | No computed real-world benefit or zero imputation |

The corresponding receipts are [assembly-report.json](assembly-report.json), [motion-report.json](motion-report.json), [site-report.json](site-report.json), [evidence-report.json](evidence-report.json), [operation-report.json](operation-report.json) and [outcome-report.json](outcome-report.json).

## The architectural decision has a measurable access consequence

The reference span between the two root positions is 23 m. Each half is a fully braced box/deck truss. Its bracing sits **below** the walking surface. This resolves the geometric conflict that would arise if transverse diagonals occupied a through-truss pedestrian passage; it does not establish the strength of the floor, rails or their connections.

The shared truss depth is `138/55 m`, approximately 2.5091 m. The lowest truss-node plane is 5.5332 m above the model road plane. A declared 0.05 m member-envelope radius puts the nominal lowest steel envelope at 5.4832 m: the earlier 5.2832 m clearance reference plus the declared 0.1 m overlay and 0.1 m reserve. Adding the real truss depth and a 0.1 m floor to the node plane places the walking surface at exactly `223913/27500 m`, approximately **8.14229 m**. The envelope radius is a research bound, not a selected member profile.

That value replaces v2's 6.0832 m illustrative walking level, which assumed only 0.6 m of structure below the deck. The increase is approximately 2.0591 m. Reusing the old approaches would make the integrated design internally inconsistent. Even on a flat reference plane, the existing 7.5% research ramp slope would require approximately **108.56 m of sloped horizontal projection per side**, before landings. Terrain and an actual approach arrangement must be applied separately.

The service elevation remains a reference design input, not a measured Tillamook clearance. The Oregon source applicability and unresolved design-basis questions remain in [oregon-design-basis.md](../v2/oregon-design-basis.md).

## What each interface must carry

| Interface | Shared content | What agreement establishes | What remains outside it |
| --- | --- | --- | --- |
| Design → assembly | Node coordinates, member connections, material/section hypotheses, floor and rail mass items | Every structural case is derived from one declared assembly | Selected member, connection, deck and rail capacities |
| Assembly → motion | Physical node IDs, member geometry, mass ledger and centres of mass | The moving geometry and support moments refer to the same mass and assembly | Real lifting equipment, actuator forces, restraint hardware and tolerances |
| Motion → site | Model coordinates, whole-half transforms, continuous swept bounds and final footprint | The fixed map placement is checked against the declared moving envelope | Surveyed clearance, complete utilities/obstacles, legal access and work-zone approval |
| Reports → evidence | Design/model/source identity and explicit candidate-state predictions | A result uses a single compatible candidate state across the supplied observations | Completeness of physical candidate states, calibrated measurements and sensor independence |
| Evidence + receipts → operation | Bound support/separation/transfer results and lock-evidence outcome | The finite state graph uses the intended guard values | Physical controller completeness, hardware interlocks and field authorization |

A hash binds bytes or canonical content to an identity. It does not establish that a source describes the as-built object, that a sensor is accurate or that every hazard has been modeled. Missing or inconsistent bindings must remain an explicit failure or unresolved condition.

## The deployment architecture is explicit

The current movement is an **ideal lift, slew and short-slide reference**, using two rigid half-trusses. It is not a telescoping structural design. Both halves initially point back along the road. The prescribed sequence raises them, turns them toward the center while their roots are retreated outward, then translates them inward for nominal engagement.

This choice makes a bounded kinematic question testable using the actual assembly geometry. It does not claim that the packaged trusses, masts, chassis or actuators fit a transport envelope. Those components require a subsequent design using the derived forces, space and mass requirements.

The workbench draws the motion report's exported poses. Its playback rate has no physical timing meaning. Motion samples show the construction; they do not prove safety between samples. Continuous support and separation conclusions must come from the analytic checks in the motion receipt. Half-to-half separation of the inflated member envelopes is scoped outside a reserved 0.15 m joint-contact zone per half. Actual engagement and contact inside that zone remain unresolved; no collision-free joint geometry is inferred. Likewise, a two-dimensional side view collapses the along-road coordinate; apparent overlap in that projection is not a three-dimensional collision test.

The support and connection conditions change across deployment and service. A final pose with matching connector locations does not establish engagement, locking or a safe support transfer. The structural release certificate and the operation guard are separate interfaces for that transition; they do not supply a physical release actuator or a recovery procedure.

## Site screening retains inconvenient evidence

The viewer embeds the real v2 Atlas subset: 395 building outlines, 401 parcel parts and 244 street/path records, with source provenance. It uses a fixed selected research placement. It does not quietly move a base or shorten an approach to avoid an intersection.

The displayed placement is Main Avenue at 7th Street, with a −30 m station offset and the elevator/stair plan allowances extending in opposite along-road directions. It has no direct building/vehicle-centerline hit and no possible swept-bound building hit in the tested horizontal screen. **Two buildings enter the additional 3 m sensitivity**, and five parcel parts intersect the base/access footprints. The 3 m distance is a sensitivity test, not a measured mapping-error bound. Ten of the 18 tested placements have possible swept-bound building intersections. These records do not identify a selected or approved crossing.

The site's road-envelope width is a hypothesis. The source street widths are authored display proxies, not surveyed pavement, right-of-way or usable roadside widths. Parcel contours do not grant access, and missing mapped features do not establish clear ground. Source terrain elevations and their quantization do not supply a survey accuracy or resolved vertical datum.

The motion envelope, base footprints and access allowances must be checked as separate, explicitly bounded objects. A continuously enclosing swept box can conservatively intersect a mapped building even when not every point in the box is occupied. Such a hit is unresolved spatial evidence requiring refinement or an explicit design change; it is not a physical collision reconstructed from building heights the source does not contain.

The preceding [v2 design review](../v2/design-review.md) remains the record for the broader 648-scenario access study. This version does not relabel those earlier calculations as a new search at the higher walking surface.

## Evidence is more than a commanded state

The evidence module distinguishes four outcomes: `condition_established`, `condition_fails`, `ambiguous` and `model_mismatch`. A command to close a lock is not an independent observation that it closed. Every surviving candidate must satisfy the declared condition before the finite model can establish it. Incompatible measurements or an invalid binding must not be turned into success through an empty candidate set.

The supplied cases are synthetic fixtures. They test the inference rules, intervals, missing measurements and data bindings. They are not sensor trials, a calibrated error model or evidence that an actual connector is locked. Candidate-library completeness and physical channel independence remain external obligations.

The operation module uses the selected synthetic case as a guard fixture while preserving that provenance. A reachable model `service` state and public opening are different conclusions. Public-use gates remain unresolved, and the interface must continue to display them even when the supplied finite graph has no invariant violation.

## The recovered Safe Roads work supplies the outcome boundary

The earlier road-safety program sought usable travel for people of all ages and abilities, with traceable evidence. The bridge's relevant outcome is a complete accessible journey. Ascent, descent, queueing, detour, opening hours, failures and the fallback route belong in the evaluation alongside the span.

The corrected local harm accounting is:

\[
H_0-H_1=N_0r_0-N_1\left[au\,r_b+(1-au)r_f\right]-D-S.
\]

Here `N0` and `N1` are baseline and intervention person-crossings over one site and period. Availability `a` is weighted by demand and requires accessible use; `u` is uptake conditional on that availability. The baseline, facility-user and fallback risks are separate. `D` and `S` account for additional, nonoverlapping deployment and displacement harm. Deaths and serious nonfatal injuries are evaluated separately.

This bridge model supplies none of the missing local exposure or causal-effect inputs automatically. Unknown quantities remain unknown. The recovered national mortality totals, correlation table and policy targets cannot stand in for the selected crossing's risk or benefit. Scenario intervals are not confidence intervals.

The compact [outcomes.py](outcomes.py) adapter implements this accounting and adds explicit unknown-input and design/site binding rules. Its [report provenance](outcome-report.json) identifies the inspected recovered `local_crossing_harm` source in `gap_closure/model_repair/models.py`, SHA-256 `d2c6a478d8276f376367ae80dfff1e439a9a214844e22d93ce3413979283412f`. The original recovered Safe Roads archive is outside this packet. Its historical evidence consisted of attributed retrieval notes and saved corrections, not complete original conversations. No private archive link or assumed empirical effect is a dependency of this adapter.

## The next physical iteration must answer the exposed interfaces

| Decision or measurement | Concrete input required | Downstream use |
| --- | --- | --- |
| Lift/slew/slide equipment | Actuator geometry, rated forces, pressure/drive limits, mass, backlash, failure and stopping behavior | Replace ideal prescribed motion and update the common mass/support model |
| Central connection | Engagement geometry, bearing/contact paths, lock tolerances, rated forces and independent state measurements | Replace nominal capture and synthetic lock evidence |
| Support transfer | Physical restraint/release mechanism, controllable reaction path, tolerances and interrupted-release response | Validate the structural release path and operation transition |
| Member and foundation design | Selected material/sections, imperfections, connection behavior, 3D stability, pad/bearing/sliding and soil inputs | Establish capacities for the forces and reactions already derived |
| Site and access | Surveyed road/ground/utilities, available space and permissions; a complete accessible approach | Resolve recorded intersections and the higher walking surface |
| Service and recovery | Load envelope, crowd/vibration response, wind limits, inspection, power loss, evacuation and reverse/recovery sequence | Replace unresolved operating gates with reviewed procedures and evidence |
| Public benefit | Observed demand, desire lines, access availability, uptake, fallback behavior and a feasible comparison | Evaluate the complete crossing service without inheriting unsupported national forecasts |

The useful next prototype is a controlled off-road apparatus that replaces one ideal interface at a time while preserving the shared design identity and its evidence. A new actuator, changed truss depth, heavier lock or different approach must propagate into the same assembly, motion, site and operation records.

## Viewer verification boundary

The page is self-contained and makes no network requests. It reads embedded report/map data, exposes keyboard controls and includes a responsive mobile layout. Browser rendering and touch behavior require a real browser check; source syntax and a DOM/Canvas adapter do not establish visual or mobile correctness. Eleven data/binding/static checks passed, and the embedded JavaScript executed all 27 poses, three stage buttons, five layer toggles, seven view controls, playback and keyboard/pointer handlers in two explicitly simulated canvas sizes without runtime errors. These are adapter checks, not actual browser or mobile rendering tests. The final [viewer receipt](verification-viewer.json) records the tested scope and exact hashes.

An independent GDAL/GEOS replay in [verification-site.json](verification-site.json) checks all 216 stored static scenarios and 18 swept placements against the pinned source geometries. All 1,170 geometric hit sets and 1,404 hazard-feature sets agree, as do the separate 1,404 coverage and 1,404 status comparisons. No invalid geometry or mismatch was found. This independent check does not rerun layout construction, establish terrain extrema, convert conservative boxes into physical collisions or certify source accuracy.
