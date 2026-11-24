# Exact axial mechanics in v2

This module adds a small, usable **2D/3D linear elastic pin-jointed truss solver**. It calculates displacements, axial forces, axial stresses and support reactions from a completely declared model. Its answers are exact rational numbers. It does not determine whether a real bridge, vehicle, mast, connection or deployment is safe.

The implementation and default verification use Python's standard library. Optional SymPy checks are independent cross-checks; they do not change the solver's supported geometry.

## Run it

```bash
python research/deployable-crossing/v2/verify_truss.py
python research/deployable-crossing/v2/verify_truss.py --symbolic
python research/deployable-crossing/v2/bridge_examples.py
python research/deployable-crossing/v2/truss.py research/deployable-crossing/v2/examples/service-24m-plane.json
```

`bridge_examples.py --write` regenerates the two illustrative input models, their exact outputs and the example report under `examples/`. It does not modify the preceding release. Running the verifier with `python -S` confirms that its default mode does not require site-packages.

## Input contract

The public entry point is `solve_truss(case)` in `truss.py`. The case contains exactly:

| Field | Meaning |
|---|---|
| `id` | Nonempty model identifier |
| `dimension` | Integer `2` or `3` |
| `units` | Explicit `length` and `force` labels |
| `nodes` | Node identifiers mapped to coordinate lists of the selected dimension |
| `members` | Records with `id`, `start`, `end`, `young_modulus`, and `area` |
| `restraints` | Records with `node`, zero-based `axis`, and explicit `displacement` |
| `loads` | Records with `node` and a force `components` list |

Coordinates and prescribed displacements use the declared length unit. Area uses its square; Young's modulus uses force/length². Unit labels are recorded without conversion: the caller must provide consistent quantities. Duplicate loads at a node are added. A support can have a nonzero prescribed displacement, allowing an exact imposed-settlement calculation. Supports constrain selected translational coordinates; general directional restraints and multipoint constraints are not implemented.

All scalars must be integers, `Fraction` values, or exact rational/decimal strings. A JSON value such as `0.003` is a floating-point number and is rejected; use `"0.003"` or `"3/1000"`. There are no material, cross-section, load or allowable-stress defaults.

For each member, the squared length is calculated from its coordinates. Its numerator and denominator must both be perfect squares. An unsupported irrational length raises `UnsupportedGeometryError` with the exact squared length; it is never silently rounded. General rational-coordinate trusses can have irrational lengths, so this is a substantive limitation. The 3–4–5 panels used here deliberately remain inside the rational kernel.

The solver does not infer self-weight, crowds, wind, equipment loads, load factors or combinations. An absent nodal force is zero in the declared model, not evidence that no physical load exists there.

## Equations and certificates

For member \(e=(i,j)\), let

\[
d_e=x_j-x_i,\qquad L_e=\lVert d_e\rVert,
\]

and let the unnormalized rigidity row \(r_e\) contain \(-d_e\) at node \(i\), \(+d_e\) at node \(j\), and zeros elsewhere. The first-order extension is

\[
\Delta_e=\frac{r_eu}{L_e}.
\]

Consequently,

\[
K=\sum_e\frac{E_eA_e}{L_e^3}r_e^{\mathsf T}r_e,
\qquad
N_e=\frac{E_eA_e}{L_e}\Delta_e,
\qquad
\sigma_e=\frac{N_e}{A_e}.
\]

Positive \(N_e\) denotes tension. Negative \(N_e\) denotes compression; it is not a buckling check.

Partitioning the coordinates into free \(f\) and prescribed \(c\) gives

\[
K_{ff}u_f=f_f-K_{fc}u_c.
\]

Every modulus and area must be strictly positive. Therefore, \(K_{ff}\) is positive definite exactly when the supported rigidity matrix \(R_f\) has full column rank. The implementation checks that rank before solving. A unique result certifies a unique displacement solution to these linear axial equations, with these ideal restraints.

If the rank is deficient, the result has `status="mechanism"` and includes an exact basis of supported first-order motions. Each witness has zero displacement at restrained coordinates and zero first-order extension in every bar. Its virtual work against the free load distinguishes two cases:

- Nonzero work for any witness: the declared load is incompatible with linear static equilibrium.
- Zero work for every witness: equilibria exist, but the displacement is not unique.

No arbitrary particular displacement is presented as a unique solution. These are infinitesimal-mechanism statements; neither rank result establishes finite-motion behavior or nonlinear stability under prestress.

For a unique solution, the output includes all node displacements, constrained-coordinate reactions \(Ku-f\), member extensions/strains/forces/stresses, the total elastic energy, and free-coordinate residuals. It also checks total external force and moment equilibrium using the **reference geometry**, as appropriate to the linear model. Member energy is checked against the assembled energy:

\[
U=\sum_e\frac{E_eA_e}{2L_e}\Delta_e^2
=\frac12u^{\mathsf T}Ku.
\]

Prescribed support movement can perform work. The reported applied-force work and support-reaction work therefore satisfy \(2U=u^{\mathsf T}(f+R)\); the support-work term must not be dropped in a settlement case.

`supported_self_stress_dimension` is the dimension of member-force combinations balanced at the free coordinates, with support reactions allowed. It is a static redundancy count, not a reserve-capacity measure.

## Independent verification

`verify_truss.py` builds its comparison stiffness from normalized unit directions and \(EA/L\), separately from the kernel's unnormalized rigidity construction. For small controls it solves by determinant expansion and Cramer's rule, rather than reproducing the kernel's row reduction. It also checks member compatibility, constitutive equations, joint equilibrium, total force/moment balance and energy.

The hand-derived controls include:

1. The supported 3–4–5 triangle, with known displacements, forces and reactions.
2. A 3-by-4 rectangle with both 5-unit diagonals: one redundant member-force degree, whose response requires compatibility.
3. A spatial control with three non-axis-aligned 5-unit bars and coupled 3D displacement components.
4. Two collinear bars with imposed support movement, including a force applied at a restrained node.

Additional controls cover an incompatible loaded mechanism, an unloaded mechanism, a free body's rigid motions, member orientation reversal, equivalent nodal-load decomposition, malformed inputs, unsupported irrational lengths, and altered solution/mechanism certificates. Optional SymPy checks independently solve and determine ranks of the same small systems.

## Illustrative bridge planes

These examples exercise the solver on a larger structural topology. **The 24 m span and 12 m reach are not dimensions extracted from Tillamook mapping.** They are declared, synthetic geometry controls. The mapped site model must establish its own required span, clearance, approaches and placement constraints.

Both examples use 3 m panels, a 4 m chord separation, 5 m diagonals, \(E=200{,}000{,}000\ \mathrm{kN/m^2}\), and area \(A=0.003\ \mathrm{m^2}\). The modulus and area are illustrative inputs, not a verified member selection.

### Connected 24 m service plane

The eight-panel plane has lower and upper chords, verticals at all panel stations, and diagonals mirrored about midspan. The lower left bearing is a pin and the lower right bearing a vertical roller.

The declared pressure inputs are 5 kPa crowd and 1 kPa dead load over a 2 m deck, shared equally between two independent planes:

\[
w_{\text{plane}}=(5+1)\frac{2}{2}=6\ \mathrm{kN/m}.
\]

Tributary nodal loads are 9 kN at each end and 18 kN at the seven internal lower nodes. They total 144 kN per plane. The bearing reactions are 72 kN each.

Top and lower midspan nodes are each shared in the input graph. This enforces ideal translational compatibility at **both chord levels** and allows an axial-force couple across the depth. It is not a single center hinge, nor a verification of any real latch, lock or joint.

The computed lower-midspan displacement is

\[
(u_x,u_y)=\left(\frac{459}{400000},-\frac{3321}{400000}\right)\mathrm m.
\]

The downward component is **8.3025 mm**. No permissible-deflection or capacity comparison is implied.

### Separate 12 m prejoin half

The four-panel half uses the same panel geometry. Both translations are held at each of its lower and upper root nodes, representing an ideal rigid mast interface. Only the prescribed 1 kPa dead-load input is used, giving 1 kN/m per plane. There is no public crowd load in this prejoin case. The dead input is not a mass calculation from member area or equipment; density and an actual dead-load inventory remain separate requirements.

End nodal loads are 1.5 kN and internal loads are 3 kN, totaling 12 kN per plane. Reactions are

\[
R_{B0}=(18,1.5)\ \mathrm{kN},\qquad
R_{T0}=(-18,10.5)\ \mathrm{kN}.
\]

Their horizontal forces form a 72 kN·m counterclockwise couple over the 4 m depth. This balances the clockwise load moment

\[
\frac{w\ell^2}{2}=\frac{1\times12^2}{2}=72\ \mathrm{kN\,m}.
\]

The lower-tip displacement is

\[
(u_x,u_y)=\left(-\frac{27}{160000},-\frac{781}{800000}\right)\mathrm m,
\]

whose downward component is **0.97625 mm**. The ideal mast has been assumed capable of imposing its restraints; its capacity and the supporting vehicle's ground reactions have not been calculated by this truss example.

### A separate check of the bridge results

`bridge_examples.py` verifies all 50 reported member forces by the method of sections. For the service plane, panel-station moments are

\[
M_j=\frac{wp^2}{2}j(n-j),
\]

and panel shears are \(wp((n-1)/2-i)\). Corresponding cantilever moments are \(wp^2(n-j)^2/2\). Chord forces follow from moment divided by the truss depth; diagonal and vertical forces follow from vertical equilibrium.

A separate unit downward load is then placed at the service midspan or prejoin tip. The independently derived member forces \(n_e\) give the downward displacement through virtual work:

\[
\delta=\sum_e\frac{N_en_eL_e}{E_eA_e}.
\]

Both deflections agree exactly with the stiffness solver. The example report records these checks, load totals, reaction moments and exact displacements; complete inputs and outputs are retained alongside it.

## Remaining mechanics scope

Two independent 2D planes do not constitute a spatial bridge model. Deck load distribution, transverse bracing, lateral stability, torsion and the actual joining mechanism remain unmodeled in these examples. The generic 3D truss control verifies coordinate handling; it does not supply that missing bridge topology.

Neither the solver nor these two endpoint configurations cover member buckling, beam bending, material strength, second-order effects, vibration, fatigue, joint/mast/chassis capacity, foundations, partial extension, lifting or joining transients. Those are specific remaining analyses. Exact arithmetic removes numerical rounding from the declared equations; it does not validate the equations, measured inputs, omitted loads or hardware.
