# One physical assembly across structural states

This packet is an exact research model of a **fully braced deck truss below a walking surface**. It is not a selected bridge design, a load rating, a capacity check, a validated deployment controller or permission to fabricate or operate equipment. The dimensional, material, mass and load inputs in `design.json` are explicit reference hypotheses. They are not measured properties of a vehicle or bridge and are not inferred from the Tillamook map.

The purpose is to close a specific modeling gap: independent cantilever, joined and partially joined analyses must use the same physical members and mass. Connecting the halves changes kinematic constraints; it does not make bars, floor or joint mass disappear.

## Geometry and inventory

The installed reference centerlines have global coordinates (0\leq x\leq23\) m, with (y) across the bridge and (z) upward. Each half has four panels. A scaled Euler brick gives

\[
(p,w,h)=\frac{23}{2200}(275,252,240)
=\left(\frac{23}{8},\frac{1449}{550},\frac{138}{55}\right)\mathrm{m}.
\]

Its face diagonals are rational because

\[
275^2+252^2=373^2,\qquad
275^2+240^2=365^2,\qquad
252^2+240^2=348^2.
\]

There are no body-diagonal bars. All longitudinal-face rectangles and station cross-sections are triangulated. The transverse diagonals remain below the walking plane; this is a deck-truss reference, not a through-truss with diagonals across a pedestrian opening. Deck plate and guard behavior are not represented by axial bars.

`L_i_j_k` and `R_i_j_k` identify physical nodes, with station (i=0,\ldots,4), side corner (j=0,1), and height corner (k=0,1). The local coordinates of either half are

\[
(ip,(j-1/2)w,kh).
\]

The left placement is the identity. The right placement is the proper rigid rotation and translation

\[
(x,y,z)\longmapsto(23-x,-y,z).
\]

Consequently the four ideal seam pairs are `L_4_j_k` with `R_4_(1-j)_k`. This orientation is a rotation, not a reflection. Each half has 20 nodes and 57 physical bars: 16 longitudinal, 10 transverse, 10 vertical, eight diagonals in each longitudinal face family and five transverse diagonals. The complete assembly has **40 physical nodes and 114 physical bars in every mode**.

Perfect ties alias solver coordinates only. In particular, two physical end-frame bars that become coincident after joining remain distinct bars with their own area, stiffness and mass. A physical joint fitting is retained on each side of every paired connector node. The solver is not a detailed center-lock model.

### Elevations and access penalty

The 2 m wide declared floor is 0.1 m thick and has its walking plane at local (z=h+0.1=287/110) m. The input named `truss_bottom_elevation_m` means the **bottom node centerline elevation**, not the underside of a real section. The separately declared member-envelope radius is 0.05 m. It is a geometric hypothesis for subsequent swept-envelope screens, not a section property derived from axial area or a strength selection.

The resulting reference elevations are

\[
z_{\rm lower\ centerline}=5.5332\mathrm{m},\qquad
z_{\rm nominal\ lower\ envelope}=5.4832\mathrm{m},
\]

\[
z_{\rm walk}=5.5332+\frac{138}{55}+0.1
=\frac{223913}{27500}\mathrm{m}\approx8.14229\mathrm{m}.
\]

The deeper rational truss therefore increases the reference access rise. No earlier shallow-understructure assumption is retained. A 0.15 m local seam zone is reserved per half for separately unresolved joint contact and engagement. The actual section envelope, joint geometry, fabrication tolerances and any contact inside that zone remain open engineering inputs.

## One true mass ledger and one gravity prescription

The elastic inputs are (E=200{,}000{,}000\ \mathrm{kN/m^2}), declared positive areas by bar role, density (7850\ \mathrm{kg/m^3}), and gravity (9.80665\ \mathrm{m/s^2}). These are hypotheses for the reference computation; the areas have not been checked for strength, slenderness or buckling.

For every physical bar,

\[
m_e=\rho A_eL_e.
\]

Its true mass point is its reference midpoint and its vertical weight is divided equally between its endpoint nodes. Each node has a declared 5 kg fitting. Floor mass is 60 kg per square metre of the 2 m floor; rail mass is 12 kg per metre on each side. Floor mass is at the plate mid-thickness. The declared rail mass line is 0.55 m above the walking surface. Neither mass prescription establishes a floor or railing design.

| Physical category | Exact full-bridge mass, kg |
|---|---:|
| Axial members | (3347397/500=6694.794) |
| Floor | (2760) |
| Rails | (552) |
| Physical node fittings | (200) |
| **Total bridge** | **(5103397/500=10206.794)** |

The two halves have identical mass, (5103.397\) kg each, and local true centers of mass

\[
\left(\frac{23}{4},0,\frac{480060393}{280686835}\right)\mathrm{m}.
\]

Floor and rail forces transfer to top-chord nodes through exact tributary weights. For each mass item these weights conserve mass and both horizontal first moments. Thus they conserve the total vertical gravity force and its reference moment. They do **not** put elevated floor or rail mass physically at a truss node. The true three-dimensional mass point remains in `local_mass_items`; motion and support calculations must use that point, not a centroid reconstructed from force-node heights.

Each base also has a separate declared ledger: 8000 kg chassis, 2500 kg stationary fixture and 10000 kg ballast, at explicit coordinates in its base frame. These point masses are not included in bridge stiffness or bridge nodal loads. The complete declared two-base-plus-bridge mass is (51206.794\) kg. The ledger does not assert that unmodeled actuators, locks or other actual hardware weigh nothing; their properties remain unknown.

The applied load cases are dead load, dead plus 5 kPa over the complete 23 m by 2 m floor, and dead plus 5 kPa over its global-(y>0) half. The additional crowd resultants are respectively 230 kN and 115 kN downward. The asymmetric crowd centroid is at (y=0.5) m. Exact transverse tributary weights transfer that eccentric load to the wider truss; neither a deck plate capacity nor the physical transfer detail is assumed proved.

## Linear mechanics and its scope

For a bar joining reference nodes (a,b), set (d=x_b-x_a), (L=\sqrt{d\cdot d}). The unnormalized rigidity row (r_e) has entries (-d) at node (a) and (d) at node (b). The unchanged v2 solver assembles

\[
K=R^T\operatorname{diag}(E_eA_e/L_e^3)R.
\]

With prescribed displacements eliminated, positive areas and moduli give

\[
u_f^TK_{ff}u_f
=\sum_e\frac{E_eA_e}{L_e^3}(r_{e,f}u_f)^2.
\]

Therefore full column rank of the supported rigidity matrix certifies a unique solution of these linear axial equations. Its nullspace certifies first-order mechanisms. Neither result proves nonlinear stability, compression-member buckling resistance, a second-order equilibrium path, fatigue life, connection capacity or dynamic performance.

From a displacement solution,

\[
\delta_e=\frac{d\cdot(u_b-u_a)}{L},\qquad
N_e=\frac{E_eA_e}{L}\delta_e,\qquad
\sigma_e=N_e/A_e.
\]

Tension is positive. Reactions are (Ku-f) on prescribed coordinates. Their sign is the force **on the bridge**; the reported `bridge_on_base` is its negative. Moments use the right-handed convention (r\times F), in kN m. Base moments are about each outer root's reference center ((0,0,0)) or ((23,0,0)); separate base/world translations must use the usual moment-transport formula. They are not already moments about the vehicle center or a ground pad.

## Restraint and join states

Permanent service restraints are four bottom-root vertical translations, two lateral translations along the same global edge at opposite ends, and one longitudinal translation. The ideal mast-clamp state additionally holds all translations of the four root nodes of each half.

| State | Seam ties | Free DOFs | Exact supported rank | Linear conclusion |
|---|---:|---:|---:|---|
| Independent clamped halves | 0 | 96 | 96 | Unique solution |
| Joined, roots still clamped | 4 node pairs | 84 | 84 | Unique solution |
| Joined, permanent service bearings | 4 node pairs | 101 | 101 | Unique solution |
| Lower-only join, temporary restraints released | 2 lower node pairs | 107 | 106 | One first-order mechanism |
| Unjoined, temporary restraints released | 0 | 113 | 108 | Five first-order mechanisms |
| Halfway commanded support release | 4 node pairs | 84 | 84 | Unique solution with prescribed translations |

The fully joined service state is solved for all three load cases. The invalid modes have nonzero gravity virtual work in their mechanism space, so the declared linear static loads are incompatible with equilibrium; no displacement solution is fabricated for them.

There is also a direct witness for the lower-only join. For each physical point define the velocity

\[
v_L(x,y,z)=(z,0,-x),\qquad
v_R(x,y,z)=(-z,0,x-23).
\]

These are rigid-half pitch velocities. They satisfy every permanent restraint and agree at both lower seam nodes. They produce no first-order bar extension but disagree at the upper seam nodes. Their gravity virtual work is

\[
f\cdot v=W_{\rm bridge}\,\frac{23}{4}>0.
\]

This independently proves the first-order failure of that declared support/join state. It does not assert a finite folding path while fixing the exact reference endpoint distances; infinitesimal mechanisms and finite motion are different questions.

For the independent dead-load halves, a separate statics check gives each root vertical reaction (W_{\rm bridge}/2), with root pitch moments

\[
M_{y,L}=-\frac{W_{\rm bridge}}2\frac{23}4,
\qquad M_{y,R}=+\frac{W_{\rm bridge}}2\frac{23}4.
\]

These are checks on the complete same-design member/floor/rail/fitting ledger, not a substituted uniform dead-pressure assumption.

## A controlled ideal restraint-release identity

Joining the halves while roots are clamped and suddenly removing every temporary restraint are different operations. The report provides a narrower, checkable mathematical construction.

Let (u_c) be the fully joined dead-load solution with every root translation held at zero, and (u_s) the fully joined dead-load solution with only permanent bearings. For each of the 17 temporary root degrees of freedom prescribe

\[
u_t(\lambda)=\lambda(u_s)_t,\qquad0\leq\lambda\leq1,
\]

while permanent bearing displacements remain zero. Keep the same reference stiffness, ties, masses and dead loads. Linearity and uniqueness then give

\[
u(\lambda)=(1-\lambda)u_c+\lambda u_s.
\]

Member forces and reactions are affine in the same parameter. Temporary reactions satisfy (r_t(\lambda)=(1-\lambda)(r_c)_t), since they vanish in the service endpoint. The implementation performs an actual solve at (\lambda=1/2) and checks all displacements, bar forces and reactions against this identity. The verifier independently checks the endpoint's free-root equilibrium, including the zero temporary reaction result.

This is a quasi-static commanded-support-displacement path in one linear model. It supplies neither an actuator arrangement nor a force/stroke/control certification. It does not prove that instantaneous release is acceptable or that locks can carry the extracted demands. Dynamic effects, yielding, buckling, contact and geometric nonlinearity are outside the proof.

## Connector demand and capture output

For each physical seam node the module recovers the sum of incident member force gradients (g_n) before aliasing the two halves. The connector force on that physical half is

\[
c_n=g_n-f_n
\]

at the seam, where no support reaction is prescribed. A tied pair must have (c_L+c_R=0); an untied free seam node must have (c_n=0). The report checks these identities and gives each force vector and the relative displacement. This is a force demand under ideal perfect tying, not a lock or contact capacity.

The report also supplies the independently clamped dead-load half shapes in their local frames, under `capture_reference.deformed_local_nodes_by_half`. The separate motion module can transform those calculated tips instead of assuming undeformed perfect alignment. These are small-displacement axial-model shapes, not surveyed hardware or a contact simulation. Horizontal yaw preserves gravity direction; other transformations or changed support loads would require the applicable structural state to be recalculated.

## Data contract

All numerical result scalars are exact rational strings. Counts, indices and booleans remain ordinary JSON primitives. No floating-point value is accepted as a design scalar.

`load_design(path=...)` validates and returns the input dictionary. `build_assembly(design)` returns:

- `design_sha256`: SHA-256 of sorted compact UTF-8 JSON for the declared design.
- `dimensions`: centerline span, panel width/depth, floor dimensions, declared member envelope and derived walking/lowest-member elevations.
- `halves.L` and `halves.R`: `local_nodes`, installed `global_nodes`, physical `members`, `local_mass_items`, mass, true local/installed CG, and proper installed rotation/translation.
- `local_mass_items`: each has `id`, `half`, `category`, `mass_kg`, `position_m` and `local_center_m` (identical local true CG), `installed_center_m`, and exact `gravity_mass_allocation_kg` to structural nodes.
- `physical_nodes`, `physical_members`, `mass_ledger`, `nodal_dead_forces`, `connector_pairs` and explicit `mode_definitions`.

The stationary base mass ledger in `design.base_interfaces.fixed_mass_items_per_base` is in base coordinates and remains separate from each movable bridge half. The bottom-centerline elevation, local pivot offset, nominal pad rectangle and fixture-force hypotheses are explicit inputs for other modules; the truss solver itself does not prove their contact feasibility or capacity.

`run_assembly(design)` returns `assembly`, eight `cases`, `release_certificate`, `capture_reference`, source hashes and scope flags. `geometry_and_mass_source_sha256` is the same canonical hash applied to the complete `assembly` object. Every case retains the common design digest and physical inventory. Stable cases contain physical node displacement/reaction dictionaries, physical member results, per-base resultants, connector demands and response extrema. Member results are matched by physical member ID; solver endpoint names may be aliased at the seam. Mechanism cases contain exact physical-node witnesses and load virtual work instead of a fictitious equilibrium displacement.

## Independent verification

Run from this directory:

```sh
python -B assembly.py --output assembly-report.json
python -B verify_assembly.py
python -B verify_assembly.py --replay
```

The default verifier does not call the truss solver. It checks the physical inventory, rational metrics, independent category masses, true mass CG and gravity moments; reconstructs the nodal gravity/crowd load prescription from closed-form tributary formulas; recomputes every bar's extension, axial force, stress and energy; verifies nodal equilibrium, total force/moment balance, work-energy equality, reaction signs and connector balance; and checks the release identity and local capture mapping.

For rank, the verifier builds rigidity rows directly and obtains a lower bound by elimination modulo the prime 1000003. A nonzero minor modulo that prime is a nonzero rational minor. In regular states this lower bound equals the number of free DOFs. In singular states an independently checked exact nullspace basis supplies the matching upper bound. Thus these are rank certificates, not floating-point singular-value tolerances. The analytic rigid-half pitch field and cantilever moments are additional controls with a different derivation.

Twelve negative controls reject malformed exact inputs, an irrational bar metric, tampered member force, altered mechanism restraint, changed gravity load, false release displacement and inconsistent source-design bindings. `--replay` additionally regenerates the entire deterministic report and compares it field for field. A successful arithmetic receipt establishes the stated equations for this declared assembly; it does not convert the reference hypotheses into verified physical hardware.
