# Two mobile halves of a pedestrian bridge

**A feasibility study for a deployable crossing.** The proposed machine places a wheeled unit on each roadside, establishes stable supports, raises and extends two half-spans, and connects them above the road. Pedestrians then use a complete elevated crossing. Crest–Root can help analyze the geometry, deployment states, and mathematical margins of such a system. The code in this directory implements a first set of exact research screens; it is not a structural design, an operational controller, or a construction approval.

[Download the independent study and verifier](https://github.com/towre676-cloud/crest-root-calculus/raw/refs/heads/main/downloads/crest-root-deployable-crossing-2026-11-23.zip) · [Equations and proof obligations](proof.md) · [Primary sources](sources.md)

## A concrete architecture to investigate

The starting architecture has two transportable support towers, an extending truss or girder assembly at each tower, a guided center connection, and an approach system on each roadside. Wheels provide transport. An engineered support system transfers deployment and service loads to the ground. That support system might use stabilizers and pads, ballast, anchors, or prepared docking foundations; none is selected or sized here.

During deployment, each half must have a defined load path before the other half can help it. After joining, the final support conditions and the center connection must agree with the structural model. The approach system is a major part of the product: ramps, suitable elevators, or existing elevated ground have very different space, transport, reliability, and operating implications.

A sensible first use case to investigate is a temporary crossing that remains deployed for an event or works period. Deployment and retraction would occur within an approved traffic exclusion or closure. Duration of closure, roadway clearance, walking width, accessibility provisions, ground conditions, wind limits, and crowd demand remain project inputs. No deployment time or capacity has been established.

## What already exists elsewhere

The broad two-unit concept has a close published precedent. **CN106012799A**, published in 2016, describes opposed vehicles with support legs, raised platforms, extending bridge bodies, counterweights, and connectors at the meeting ends. Its Figure 6 shows the two halves forming an overpass. This is a patent disclosure, not evidence that its load capacity, deployment time, or commercial performance was demonstrated. [Original document and translation](https://patents.google.com/patent/CN106012799A/en)

**US12480266B1**, published in 2025, describes a related arrangement in which one truck raises and extends a scissor span to a receiving truck. Its document includes subsequently installed access and deck components. Its asserted deployment time is not used as a measured benchmark here. [Patent document](https://patents.google.com/patent/US12480266B1/en)

For built mechanisms, the designers of Copenhagen's **Inderhavnsbroen** describe two sliding cantilevers that meet and interlock. Their account explains changing support reactions during movement. Its substantial permanent foundations differ from portable roadside bases. APS also documents an actual modular temporary pedestrian overpass at the Australian Grand Prix, assembled using transport and cranage. These examples motivate comparison with existing engineering rather than a claim that the broad arrangement is new. [Designers' account](https://ingenia.org.uk/articles/copenhagens-kissing-bridge/) · [Temporary overpass project](https://apsystems.net.au/project-examples/grand-prix-overpass-1-temporary-modular-overpass/)

## The constraints that decide feasibility first

### Reaching the deck

Let the walking surface rise by $H$ above its approach ground level. At a declared maximum ramp slope of $1:12$, the sloping portions require at least $12H$ of horizontal run on **each** side. A hypothetical **18 ft** walking-surface rise therefore requires **216 ft** of ramp run per approach. With a 30-inch maximum rise per run, at least eight runs and seven intermediate landings are needed. Five-foot intermediate landings add 35 ft of level travel, giving 251 ft before the top and bottom landings. These dimensions illustrate the Access Board's ramp guidance; the applicable public-roadway requirements must be established for the project. [Access Board ramp guidance](https://www.access-board.gov/ada/guides/chapter-4-ramps-and-curb-ramps/)

That 18 ft is an assumed **approach-to-walking-surface rise**, not a recommended vehicle clearance. The road owner determines the applicable clearance envelope, and the bridge's lowest component and deformed geometry must respect it. Road profile, shoulders, structural depth below the deck, and approach-ground elevations affect the calculation. Switchbacks reduce the parcel length but do not remove the climb or the route length. A suitable elevator is a different system to engineer, not an automatic substitution of any available platform lift. [Public-right-of-way guidelines, R203.7 and R407](https://www.access-board.gov/prowag/complete.html)

### Stability before joining

For an ideal horizontal cantilever of reach $\ell$ carrying uniform downward line load $w$, its gravity shear and root bending moment have magnitudes

$$V=w\ell,\qquad M_{\rm root}=\frac{w\ell^2}{2}.$$

For the same line load, reducing the unsupported reach from a full span $L$ to a half-span $L/2$ reduces this particular root moment by a factor of four. It does not prove that the roadside units can resist their moments, that another deployment state is less severe, or that the final bridge is adequate.

For a planar support screen, take vertical pads at $x_L<x_R$, total downward force $W>0$, and total clockwise external moment $M$ about $x=0$. Their vertical reactions are

$$R_L=\frac{x_RW-M}{x_R-x_L},\qquad
R_R=\frac{M-x_LW}{x_R-x_L}.$$

Both reactions are strictly positive exactly when $x_LW<M<x_RW$. At equality a support has unloaded; this study labels the result **marginal**, not acceptable by default. Negative reactions require a hold-down load path that this unilateral-pad model does not provide. This screen does not check sliding, bearing pressure, soil failure, three-dimensional overturning, support strength, or a safety factor. Ground slope and load directions must also be represented appropriately in a later model.

The root moment above is a local member result. It must not be added twice when assembling the total external moment of the entire unit: the distributed bridge weight already contributes its force times its center-of-gravity coordinate. A horizontal force at height $z$ contributes a signed moment $Fz$ under the stated clockwise convention.

FHWA's guidance on moving and erecting bridges treats temporary supports, ground bearing, bracing, wind, and moving-state loads as explicit engineering problems. It supplies relevant principles, not a ready-made standard for this proposed machine. [FHWA design guidance](https://www.fhwa.dot.gov/bridge/pubs/07022/chap04.cfm)

### A structural connection at the center

The center interface has to accommodate initial alignment error and then establish the intended structural restraint. Guiding surfaces, capture features, positive retention, and verification of engagement are functions to develop; this study supplies no connector drawing or rated lock.

The required load transfer depends on the whole system. A final bridge modeled as simply supported at its roadside ends needs a connection capable of transferring the internal forces of a continuous span, including bending. Two independently restrained cantilevers could instead meet through a connection with different rotational behavior, but their bases must then supply the cantilever restraint. A plate that merely covers the deck gap establishes neither model.

The exact rigid-body rank examples demonstrate this distinction. Collinear half-spans supported by outer pins and joined by a center hinge have a **first-order mechanism**. Fixing the bases changes the answer. A rigid center joint also changes the constraints. These examples do not calculate elastic stiffness, member buckling, or finite-amplitude motion; a zero first-order stiffness at a straight configuration does not necessarily mean a finite folding motion exists with the endpoints fixed.

### Time to usable crossing

The meaningful time is from arrival to an inspected, accessible crossing and reopened road, with retraction measured separately. Parking, traffic control, stabilizing, approach installation, extension, alignment, locking, and inspection all contribute. The deployment stroke alone is not a useful product benchmark. A small prototype may establish repeatability and alignment performance but cannot establish full-scale pedestrian capacity or bridge dynamics.

## What this program now checks

The new functions use Python's exact rational arithmetic and closed intervals. There are no default bridge capacities, roadway clearances, allowable misalignments, or wind limits. Decimal inputs must be strings so their intended values are retained exactly; binary floats are rejected.

| Screen | Computed result | Explicit scope |
|---|---|---|
| Uniform cantilever | Gravity root shear and moment bounds | Declared line load and horizontal reach; no beam deflection or capacity |
| Planar pad reactions | Exact reaction bounds and a strict no-lift classification | Two vertical unilateral supports; declared force/moment intervals |
| Center kinematics | Exact constraint rank and first-order degrees of freedom | Two ideal rigid planar halves and named support/joint constraints |
| Joining geometry | Longitudinal overlap and vertical-mismatch bounds | No lateral offset, yaw, pitch, roll, elastic deformation, or connector mechanics |
| Ramp slope | Exact comparison with a declared maximum slope | No complete accessibility determination |
| Abstract sequence record | Order of declared state labels in one fixed sequence | No sensor verification, machinery control, or operational authorization |

An **inside** result means only that every tested margin is strictly positive throughout the declared input box. **Outside** means at least one tested margin is strictly negative throughout the box. **Marginal** identifies an identically zero margin when no other check is outside or indeterminate. **Indeterminate** covers the remaining unresolved enclosures, including an interval that merely touches zero without being identically zero. All classifications remain subject to the model and its stated omissions. Independent input boxes may overestimate uncertainty when their quantities are physically correlated.

## Reproduce the study

The focused ZIP is self-contained and uses Python 3.11 or later, with no third-party packages. Extract it, open a terminal in the extracted directory, and run:

```sh
py verify.py
py study.py example.json
```

On platforms where Python is named `python`, use that command instead. From the repository root, the equivalent commands are:

```sh
python research/deployable-crossing/verify.py
python research/deployable-crossing/study.py research/deployable-crossing/example.json
```

The examples are algebra and geometry controls, not a proposed bridge specification. The reported reactions and tolerance windows must not be interpreted as selected hardware properties or a load rating. The study's recorded output is in `report.json`; the verifier checks independent expected values and adverse cases. The ZIP includes a SHA-256 manifest for its contents.

The default verifier also runs with `python -S verify.py`, which disables installed site packages. An additional independent symbolic cross-check is available with SymPy:

```sh
py -m pip install "sympy>=1.12,<2"
py verify.py --symbolic
```

The archived `verification.json` and `verification-symbolic.json` identify their respective modes. The symbolic pass checks ranks, equilibrium equations at justified interval corners, and the distributed-load integral using a separate implementation.

## The next engineering challenge

The first project-specific model should declare a road and access envelope, a transport envelope, support locations and ground model, moving masses, an actual extension geometry, center-joint behavior, and service load cases. It should analyze deployment, first contact, partial engagement, the completed connection, occupied service, interrupted movement, and retraction. Every claimed pass needs a stated load path and model.

For mathematical design, let $p$ describe the design, $s$ a point on its deployment path, and $\xi$ an allowed uncertainty. For each named constraint define a margin $g_j(p,s,\xi)$ and a positive reference scale $G_j$ with the same units. Then

$$C(p)=\inf_{s,\xi}\min_j \frac{g_j(p,s,\xi)}{G_j}.$$

Proving $C(p)>0$ would establish all the **modeled** strict inequalities throughout the specified domain. Separate discrete modes require separate constraints and justified transitions. This is a formulation for research, not a computation completed by the small kernel above. Its connection to Crest–Root is exact boundary analysis, structured reduction, interval enclosure, and checkable witnesses. No new theorem or speed advantage over PARI or structural solvers is claimed.

The physical development sequence is a site-envelope study, a coupled structural/mechanism model, an instrumented non-occupied prototype, and then a professionally designed and independently reviewed full-scale test program. Public use would require the road owner's acceptance and the applicable structural, accessibility, machinery, traffic, inspection, and operating provisions. The particular standards and approvals depend on location and use.

The first inputs that would change this study most are the intended road width and clearance, available roadside space, means of reaching the deck, transport limits, and whether prepared docking sites are acceptable. Those choices determine whether rapid deployment is practically valuable at the intended sites.
