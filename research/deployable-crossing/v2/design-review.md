# Tillamook v2 — engineering decisions and remaining work

## Decision from the mapped test bed

The two mobile half-spans remain a research concept. This release turns the saved Tillamook Atlas into an auditable engineering test bed and adds actual structural analysis. It does not select a site or specify a bridge for construction.

The first site sweep covers six real mapped junctions: Main and Pacific Avenues at 7th, 9th and 11th Streets. The choice is supported by the crossing-needs corridor in ODOT's March 2025 planning memorandum, which summarizes street-level crossing improvements. Those recommendations are not endorsements of an overhead bridge. See [the Oregon basis](oregon-design-basis.md).

The sweep contains 648 deliberately declared placement hypotheses: six junctions, three road-envelope widths, three station offsets, three access arrangements and four combinations of approach directions. Of these, 636 intersect at least one mapped building footprint or vehicle centerline. The other 12 all use the elevator-and-stair allowance. All 12 have a mapped building within the additional 3 m separation sensitivity. That sensitivity is a hypothetical allowance, not a measured position-error bound.

These counts describe this finite sweep. They do not prove that no bridge can fit Tillamook: other streets, setbacks, orientations, access structures, properties and machine dimensions were not exhausted. Equally, the 12 cases do not identify available sites. They still encounter parcels, lack verified accessible connectivity, and require utility, ground, hazard and permission checks. The independent GDAL/GEOS replay reproduced all 5,184 conflict sets. See [the site receipt](verification-site-independent.json).

**Design consequence:** make the complete roadside system the subject of the next hardware design. Optimizing a telescoping span in isolation would miss the binding footprint and access questions visible in the map.

## What the Atlas establishes

The saved V5 source is pinned by its archive SHA-256 and a checked 707-file ledger. The extraction retains local metre coordinates in NAD83 / UTM zone 10N, original feature IDs, rings, source tags and a crop of the unblended terrain grid. It excludes authored façades and room geometry from obstacle evidence. The downtown subset contains 395 building footprints and 401 parcel parts, including 78 interior parcel rings.

Every one of the subset's 244 street/path widths is authored. In the full Atlas only two of 2,044 path widths come from attributes, and those are runways. Therefore a rendered road surface cannot supply our required span. The 8, 12 and 16 m road envelopes are explicit hypotheses; with the declared setbacks and support envelopes they produce 19, 23 and 27 m base-center spans. They are not measured dimensions.

Half reach in the placement report means base-center to midpoint distance. The selected coupling and bearing design would add its own engagement length, offsets and tolerances; those dimensions have not been supplied by the map.

The source stores road coordinates to 0.01 m and terrain samples at 0.025 m increments on a 2 m grid. Encoding precision does not establish survey accuracy. In particular, the terrain's vertical datum lineage and error bounds remain unresolved. Source terrain interpolation is useful for footprint sensitivity; it is insufficient for pavement high points, curbs, individual pad levels or clearances. Source sidewalk tags exist where the separate footway geometry is incomplete; missing lines do not mean sidewalks are absent. [Atlas audit](atlas-audit.md)

## Hazard overlays use the proposed footprints

All 2,592 separate base/access footprints in the sweep lie inside the captured query area. GDAL/GEOS performed 33,696 polygon-pair tests against the six official layers. The resulting [hazard report](hazard-report.json) preserves each source feature's complete attributes once in a catalogue and references the particular footprint intersections and their areas.

| Source classification | Layouts with an intersecting footprint |
|---|---:|
| FEMA AE | 213 |
| FEMA X, 0.2% annual-chance flood hazard | 166 |
| FEMA X, minimal flood hazard | 547 |
| DOGAMI evacuation layer, Outside Hazard Area | 648 |
| DOGAMI evacuation layer, Local Tsunami Evacuation Zone | 0 |
| DOGAMI local XXL1 inundation polygon | 0 |
| Detailed landslide study coverage, CalhNC2020a | 648 |
| Deep landslide susceptibility, Low | 648 |

FEMA counts overlap because a layout can touch several categories. No mapped landslide deposit was returned for the query. The zero intersection counts and the source labels are not ground or hazard-free certifications. In particular, the broader query rectangle contains local-tsunami and XXL1 polygons even though the tested footprints do not intersect them; this is why the analysis intersects the actual layouts rather than assigning a single label to the town.

Among the 12 layouts with no tested building or vehicle-centerline overlap, ten touch only the FEMA X/minimal category. Two at Main and 11th also touch other flood classifications; one includes AE. Every one of the 12 still encounters the separate 3 m building sensitivity. The map does not establish an available or authorized parcel. [Source capture, meanings and limitations](data/hazards/README.md)

## Architecture being compared

| Arrangement | What the current model evaluates | What determines whether it is practical |
|---|---|---|
| Two units with straight accessible ramps | Complete plan envelopes for the units and long approaches | Roadside corridor length, entrances, driveways, property and continuity |
| Two units with unstacked switchback ramps | Side-by-side ramp runs with terminal/intermediate landing allowances | Much wider roadside footprint, supports and changes in grade |
| Two units with elevator and stairs | A declared 20 × 6 m combined access allowance per side, plus bases | Actual compliant elevator/stair design, power, throughput, queues and rescue |
| A planned street-level crossing improvement | Documented comparator in the local planning evidence | Current pedestrian demand, traffic conditions, accessible journey and agency design |

The elevator/stair rectangle is an allowance, not a designed elevator. The switchback template is unstacked; it does not assume that folding ramps can occupy the same space without headroom and structural consequences. None of the templates includes every barrier, queue, turning path, service vehicle or utility setback.

### Height has a substantial access cost

The starting Oregon pedestrian-overpass reference is 17 ft 4 in, or exactly 5.2832 m, to the underside, excluding future-overlay allowance. For a flat-ground illustration we explicitly add 0.10 m overlay, 0.60 m of structure below the walking surface and 0.10 m clearance reserve. These are study assumptions; the last two do not come from a solved bridge deformation or foundation model.

The resulting illustrative walking-surface rise is 6.0832 m. At the 7.5% ramp slope used in ODOT's temporary accessible-route reference, the sloped horizontal run is 81.1093 m per side. Eight runs are needed with a 0.762 m maximum rise per run. Including nine 1.524 m end/intermediate landing allowances gives a 94.8253 m straight footprint per side. The unstacked eight-lane switchback envelope is 13.208 × 14.4 m, using an assumed 1.8 m outside lane width.

These are bounding geometry calculations. They do not establish compliant cross slope, turning space, rail details, headroom, load capacity or an approved choice of accessibility standard. Actual approach elevation changes the rise; the site sweep samples the terrain and expands the approach allowance in two declared passes, without claiming an extrema or convergence proof. [ODOT references and applicability](oregon-design-basis.md)

## Structural and mechanical concept review

### Treat the base as part of the bridge

The wheels provide transport. A proposed occupied configuration needs a designed load path through its supports to verified ground. A wheel brake by itself does not establish resistance to overturning, ground settlement or load redistribution. The new support-polygon calculation checks whether the resultant admits nonnegative vertical reactions at coplanar pad locations. It does not determine the actual reactions of four flexible pads, their contact pressures, sliding, torsion or anchorage capacity.

The discrete support challenge enumerates 198 combinations of reach, pad spread and ballast. It checks four bounding horizontal-force combinations across every point of a prescribed extension path. The favorable results are conditional on the declared weights, centre-of-gravity motion and forces. They are useful for finding tradeoffs and eliminating configurations, not for specifying ballast on a real machine. Even a positive contact margin needs individual-pad stiffness, ground bearing, friction and dynamic analysis.

### The unlocked half and the connected bridge are different structural states

During extension, each half must carry its own weight, attachments and environmental/dynamic loads before any help arrives from the other side. The illustrative 12 m plane is therefore modeled as a cantilever held at both root chord nodes. In its prescribed 1 kN/m per-plane load case, the root must supply a 72 kN·m resisting couple. That calculation illustrates the demand; the mast, vehicle and soil that must supply it are not designed.

The separate 24 m service example uses a pin at one end and a roller at the other, with ideal continuity at both chord levels at midspan. It is a second idealization, not the automatic service state of the modeled cantilever. Its 8.3025 mm calculated midspan deflection applies only to that exact small-displacement axial model and its inputs. The 12 m cantilever's 0.97625 mm tip result likewise excludes root flexibility, changing overlap, real joint slip and second-order response. Neither is a rated allowable deflection.

**An actual support-state transition remains to be designed and analyzed.** Holding both bases fully fixed while also locking the seam can create restraint forces under temperature changes, foundation movement or manufacturing error. Releasing a temporary restraint too early can create a mechanism. The next three-dimensional model must specify every temporary restraint, service bearing and transition, including asymmetric extension and retraction. A shared node in a solver is not evidence that a physical lock can enforce that compatibility.

### Center connection requirements

The v2 capture check applies an exact finite rotation and translation to each declared connector point. Four connectors at separated positions can miss their capture envelopes even when the midpoint looks aligned. The example deliberately demonstrates that a small pitch error can move the upper/lower connectors beyond a declared axial capture allowance.

A practical connection concept needs distinct answers for:

1. **Guiding and capture:** allowable approach position and orientation, funnels/rollers, collision exclusion, and a recoverable misalignment state.
2. **Structural transfer:** vertical and lateral shear, axial force, chord-force couple, torsion and reversals, with local bearing, tear-out, weld, pin and fatigue analysis.
3. **Positive retention:** mechanical load holding and retention under loss of actuator pressure or power; engagement geometry that can be inspected.
4. **Verification:** evidence of each required lock state and disagreement handling. A sensor saying “closed” is not proof of adequate pin engagement or connection capacity.
5. **Walking surface:** continuous support, trip/slip transitions, drainage, guarding and accessible clearances throughout the operating tolerance envelope.

The mathematical capture box is not a specified pin tolerance or a load-rated interface. No controller or automatic public-access authorization is supplied.

### Stability of the actual members and whole assembly

The new truss solver provides axial force, stress, displacement, reaction, exact rank and mechanism witnesses for rational-length 2D/3D models. Its 3D test confirms the implementation's spatial coordinate handling. It does not establish the three-dimensional stability of a bridge represented by two independent plane models.

The physical design still needs member bending, in-plane and out-of-plane buckling, local section stability, connection eccentricity, overlap/lock stiffness, geometric nonlinearity, deck action and torsional/lateral bracing. A catalogue section must provide more than cross-sectional area: its moments of inertia, radii of gyration, material properties, effective lengths, connections and corrosion/fatigue allowances affect capacity. Crowd-induced vibration, lateral synchronization, wind response and vehicle-strike scenarios require their own models and accepted criteria. The present exact solver deliberately does not turn axial stress into a capacity verdict. [Mechanics scope and verification](mechanics.md)

## Closing the operational gaps

The release does not claim a deployment time. “Rapid” must cover the full arrival-to-access process: staging, traffic control, ground checks, support setup, elevation/extension, docking, lock inspection, accessible approaches and opening checks. Retraction and abnormal recovery also need timing. Record both typical completion and adverse-case delays; a span-extension stopwatch alone would conceal much of the system work.

An engineering test specification should cover at least these independent failure conditions:

| Condition | Evidence to obtain before public operation |
|---|---|
| One unit stops or retracts unexpectedly | A stable bounded state, restraint/load-path analysis, detection and recovery demonstration |
| Electrical or hydraulic power loss | Demonstrated load holding and safe recovery without relying on unverified pressure retention |
| One pad settles or loses contact | Defined inspection limits, modelled redistribution, response procedure and ground acceptance evidence |
| One connector is short of full engagement | Independent detection/inspection and inability to authorize use from the incomplete state |
| Position sensors disagree | An explicit fault state and verified response, including avoidance of further motion into the other half |
| Wind increases during setup or service | A measured and validated operating envelope, forecast/monitoring method and evacuation/retraction plan |
| Elevator unavailable or queues exceed space | Accessible operating/recovery plan, capacity evidence and alternate route |
| Vehicle strike or excessive height vehicle | Authority-reviewed protection, load consequences, closure response and inspection criteria |
| Emergency evacuation, fire or medical event | Egress and rescue routes that remain usable in the relevant fault states |
| Repeated transport and deployment cycles | Fatigue basis, wear limits, corrosion protection, inspection intervals and traceable maintenance |

These are proposed project verification items. They are not tests already performed. The program cannot close physical evidence gaps by changing a status label.

## Recommended next engineering milestone

Use the mapped findings to define a **controlled off-road prototype requirement**, then choose a representative trial area with permission and suitable ground. The Atlas can locate candidates, but this release has not identified an authorized trial parcel. The trial should be designed and reviewed before anyone occupies the device or it moves over public traffic.

For that milestone, freeze the complete assembly: two units, real members and connections, masts, bearings, actuators, supports, deck and the intended access method. Build one spatial model with consistent mass, stiffness, load combinations and support states. Replace each major assumed dimension/load with a drawing, certified component property, survey or measurement, and independently check it. This is the shortest defensible path from the software to a machine, because it targets the unresolved items the test bed has exposed.

Public-benefit and economic work should proceed alongside mechanics: actual crossing demand and destinations; accessible journey distance and waiting time; deployment disruption; transport/crew/maintenance costs; and comparison with planned street-level improvements. Existing patent disclosures cited in the previous release mean novelty or freedom to operate cannot be assumed from this concept. No lives-saved, cost-saving, patentability or reliability claim is established here.

## Strengths and the next computational challenge

Crest–Root now has a credible role as an exact model checker and discrete design explorer. It can preserve a source-backed geometry hypothesis, reject a mechanism with an exact witness, solve supported rational-length trusses, expose intermediate path failures, and enumerate support tradeoffs without numerical rank tolerances. That complements engineering judgement and measured inputs.

No PARI/GP performance comparison was run for this v2 work. A fair next benchmark would give both programs identical rational truss matrices and support states, require the same exact ranks/solutions/witnesses, and separate graph construction, solving and certificate verification. The practical structural comparison should also include a verified conventional finite-element model with equivalent restraints and loading. Correctness and scope come before a speed claim.
