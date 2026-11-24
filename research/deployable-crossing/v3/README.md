# Tillamook crossing v3 — one design across the studies

**Version 3 connects the bridge assembly, motion, support reactions, site placement, measurement evidence and operating sequence to one shared design.** Changes can no longer quietly use one mass in the truss solver, another in the support model and a lower deck in the access study.

[Download the complete offline packet](https://github.com/towre676-cloud/crest-root-calculus/raw/refs/heads/main/downloads/crest-root-tillamook-crossing-v3-2026-11-23.zip) · [Open the workbench file](workbench.html) · [Integration review](integration-review.md) · [All results](integrated-report.json)

The new connection exposed a practical problem: keeping the fully braced reference truss below the walking surface raises that surface to **8.1423 m**. Only **two of 216** tested access/base layouts avoid the mapped building and vehicle-centerline intersections; both are elevator-and-stair space allowances near the Main/7th reference, and both have buildings within the 3 m sensitivity screen. They are unresolved hypotheses. **No site, bridge capacity, deployment time or public opening is approved.**

## Open it

Extract the ZIP and double-click **Open-Workbench.cmd** on Windows, or open **v3/workbench.html** in a browser. The page is self-contained and needs no server or network connection. Its 3D schematic, side elevation and actual Tillamook map use the same exported assembly and poses. The play control has no real deployment-time meaning.

The offline packet includes the exact reports, source code, independent checks, captured hazard layers and relevant v2 dependencies. An ordinary browser opens the workbench; Python is needed only to rerun the calculations. The recorded JavaScript checks use a DOM/Canvas adapter. Actual browser rendering, mobile touch behavior and the Windows launcher have not been tested in this execution environment.

## The shared reference

| Item | Declared or derived value |
| --- | ---: |
| Span between ideal roots | 23 m |
| Physical nodes / bars | 40 / 114 |
| Length of each rigid half | 11.5 m |
| Floor width | 2 m |
| Truss width / depth | 2.63455 / 2.50909 m |
| Bridge mass from the actual ledger | 10,206.794 kg |
| Fixed chassis, fixture and ballast allowance per base | 20,500 kg |
| Declared base pad rectangle | 5 × 12 m |
| Nominal lowest member envelope above reference plane | 5.4832 m |
| Walking surface above reference plane | 8.14229 m |

The mathematical reference uses an Euler brick with panel:width:depth ratio 275:252:240. All axial member lengths are rational, permitting exact stiffness, rank, force and energy calculations. **That ratio is a computational reference, not an optimized structural proportion or selected transport envelope.** A shallower or otherwise better physical structure would require a new shared design and analysis; its benefit must not be assumed from this reference.

Member mass derives from density, area and exact length. Floor, rail and node-fitting masses are explicit. True elevated mass centroids feed the moving-support calculations even when gravity is transferred to structural nodes. Each half weighs 5,103.397 kg; the declared system mass including both fixed base allowances is 51,206.794 kg before any omitted hardware or access equipment.

The prescribed movement is a **lift–slew–short-slide variant**. Both rigid halves begin along the road, lift, rotate inward with a 0.2 m root retreat and then approach the seam. It supplies a concrete kinematic reference for comparison with the original extending-half idea. It does not supply a telescoping mechanism, selected actuators, a mobile chassis design or deployment timing.

## What the structural solver found

Every case preserves the same physical member and mass inventory. Restraints and ideal seam constraints change by mode.

| Mode | Supported rank / free degrees of freedom | Result |
| --- | ---: | --- |
| Independent clamped halves, dead load | 96 / 96 | Unique linear solution |
| Joined while roots remain clamped | 84 / 84 | Unique linear solution |
| Joined service bearings; dead, symmetric and asymmetric crowd | 101 / 101 | Unique linear solution in each case |
| Only lower seam connections; temporary restraints released | 106 / 107 | One first-order mechanism |
| No seam connection; temporary restraints released | 108 / 113 | Five first-order mechanisms |
| Actual midpoint of commanded restraint transfer | 84 / 84 | Matches the exact affine release identity |

Both invalid modes have nonzero gravity virtual work and no compatible equilibrium in the declared linear system. The result makes a center-interface requirement concrete: nominal deck continuity or a lower-only join cannot justify the modeled release.

Each independent half requires a gravity root shear of approximately 50.0472 kN and pitch moment of 287.7716 kN m. The joined-clamped case also has horizontal root reactions. The full root wrench is retained, and its vertical-load and roll/pitch effects are transported into the same pad-plane screen; moments are not counted twice. Horizontal and yaw restraint remain outside that screen.

The controlled transfer prescribes temporary support translations proportional to the final service displacement. Within one fixed linear stiffness model, every response is affine in the release parameter. An actual midpoint solve verifies this identity, and all temporary reactions vanish at the endpoint. That proves an ideal quasi-static load-transfer path. Real jacks, constraints, lock engagement and interrupted release remain to be designed.

The largest computed vertical displacement in the symmetric crowd case is approximately 17.534 mm, with maximum absolute axial stress approximately 65.013 MPa. These are response values. No strength, buckling, fatigue, floor, connection or vibration limit is supplied or passed.

## Continuous paths and complete base wrenches

The motion module provides 24 complete support paths: two halves, three stages and four corners of a declared horizontal fixture-force box. Rational proper rotations and exact quadratic extrema cover every path parameter. The horizontal fixture forces are base-interface test loads, **not a bridge wind model or a wind-speed limit**.

It also checks 48 complete structural base-wrench cases and eight continuous commanded-release support paths. The supplied reference remains strictly inside the declared support rectangle in these calculations. The statement means that nonnegative vertical pad reactions exist for coplanar pads. It supplies no unique load at each pad, ground-bearing capacity, sliding resistance, anchorage, mast capacity or torsional restraint.

Three analytic separation certificates cover the nominal half-truss prisms. The declared 0.05 m member envelopes stay apart outside a reserved 0.15 m joint-contact zone per half. Actual contact and engagement inside that zone are unresolved. Rail entries are mass lines, not full guard collision bodies; actuator and lock bodies are absent. Sixty-eight conservative interval boxes bound the specified geometry continuously and feed the site screen. The 27 animation poses do not serve as a proof between frames.

Four nominal seam locations agree exactly. Using the independent-half dead-load displacement field gives a largest coordinate mismatch of approximately 1.0294 mm, within the hypothetical ±20 mm capture box. Foundation error, backlash, connector compliance and physical self-alignment are not part of that conclusion.

The lowest member envelope in the symmetric service-load model is approximately 5.46567 m above the flat reference road plane. It remains above the captured 5.2832 m clearance reference plus a hypothetical 0.1 m overlay, leaving approximately 82.47 mm of the nominal reserve. This is a linear-model result; a surveyed road envelope, settlement, dynamic movement and applicable owner criteria remain unresolved. [Captured Oregon design basis](../v2/oregon-design-basis.md)

## Access and Tillamook

The unchanged source Atlas contributes 395 building outlines, 401 parcel parts, 244 street/path records and six candidate references. The v3 study uses three station offsets, three approach types and four direction pairs: 216 layouts on 18 common placements. The assumed roadway envelope is 12 m wide, derived from the 23 m span, base width and declared setback. It is not a surveyed road width.

At the flat reference rise and the inherited 7.5% planning slope:

| Access quantity, per side | Result |
| --- | ---: |
| Horizontal projection of sloped ramp portions | 108.564 m |
| Runs / landings | 11 / 12 |
| Straight layout including landings | 126.852 m × 1.8 m |
| Unstacked switchback envelope | 12.917 m × 19.8 m |

Terrain sampling can enlarge those hypotheses. Two expansion passes are not a certified terrain envelope. An elevator-and-stair option is only a **20 × 6 m plan allowance**; there is no elevator, stair, queue or complete accessible-route design behind that rectangle.

Of the 216 layouts, 214 intersect a mapped building or vehicle centerline. The remaining two also avoid possible swept-bound building intersections, but each lies within 3 m of two buildings and overlaps several parcel parts. Ten of 18 deployment placements have possible building conflicts with the conservative sweep bounds. Such bounds can produce false positives, and missing height/utility data prevent a physical clearance conclusion. Source road widths are authored proxies, and parcel geometry grants no permission.

Captured FEMA and DOGAMI layers are applied to the actual v3 footprints and sweep bounds. An empty layer intersection is only a result for the queried source coverage. The independent GDAL/GEOS replay agrees on all 1,170 geometry hit sets, 1,404 hazard feature sets and 1,404 coverage/status checks. It does not improve the source survey accuracy. [Source notices](../v2/data/NOTICE.md)

## Evidence and operating sequence

The evidence module retains a common candidate state across every supplied observation. It returns one of four outcomes: `condition_established`, `condition_fails`, `ambiguous`, or `model_mismatch`. A command alone admits full, partial and missing locks. Synthetic independent channels can distinguish them; missing data, broad errors, stale bindings and incompatible observations cannot silently become success.

All 14 examples are **synthetic** and contain zero physical measurements. The finite candidate library is not claimed to contain every real failure. The operating graph consumes the actual structural and continuous-motion reports plus that explicit synthetic evidence fixture. It checks safe ordering, reachable states, accepted hold and unexpected deadlock separately. A deliberately unsafe graph produces a premature-release witness.

Model service is reachable under the supplied assumptions. Public opening remains blocked by unresolved site, capacity/foundation, access, field-measurement, hardware, work-zone and inspection/emergency obligations. The reverse/recovery path is unresolved. No actuator commands are issued.

The Safe Roads outcome adapter carries the recovered local harm-accounting formula into the same design/site record. Deaths and serious injuries remain separate. All local exposure, uptake, availability, risk and additional-harm inputs are unknown, so no casualty reduction is estimated. Exact synthetic controls demonstrate the accounting without pretending to measure a benefit. [Outcome report](outcome-report.json)

## Reproduce and inspect

From the extracted packet root, with Python 3.11 or newer:

```sh
python -B v3/verify_release.py
python -B v3/integrate.py --refresh
python -B v3/build_viewer.py
```

From the repository root, substitute `research/deployable-crossing/v3/` for `v3/`. The first command validates the saved artifacts and runs the independent exact-mechanics, motion, evidence and outcome checks. Refreshing recomputes the component reports in dependency order. Optional independent GIS replay requires GDAL/GEOS, and JavaScript execution checks require Node:

```sh
python -B v3/verify_site.py
python -B v3/independent_review.py
```

The packet records the actual verification environment and any untested surfaces. SHA-256 bindings and the archive manifest establish file identity, not physical truth or external approval.

## The next decision

The useful next physical target is an instrumented, non-occupied, off-road interface prototype. Select and measure the center lock and controlled support-transfer mechanism, then replace those ideal interfaces in this common model. In parallel, a real site/access envelope should decide whether a lower-depth structure, prepared docking locations or a different crossing arrangement can meet the intended rapid-deployment goal. A motion duration or lives-saved claim should follow evidence from the complete usable crossing.

This work demonstrates a strength of Crest–Root: explicit constraints, exact boundary and rank checks, continuous enclosures, checkable failure witnesses and traceable composition. It is not a new benchmark against PARI/GP or a replacement for structural design software and professional engineering review.
