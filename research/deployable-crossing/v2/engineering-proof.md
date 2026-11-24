# Exact support paths and finite connector geometry

These statements concern the declared mathematical models in `engineering.py`. They neither infer missing physical loads nor establish capacity. The independent verifier checks alternative analytic derivations and finite Hamilton-product rotations.

## Vertical contact existence on a plane

Use horizontal coordinates \((x,y)\), with \(z\) upward. Let \(W>0\) be total downward force and let \((Q_x,Q_y)\) be its required support first moments. A downward force \(w\) at \((x,y)\) contributes \((wx,wy)\). A horizontal force \((F_x,F_y)\) at height \(h\) contributes \((hF_x,hF_y)\). With this convention the external Cartesian moments are \(M_x=-Q_y\) and \(M_y=Q_x\).

Vertical pad reactions \(R_i\ge0\) at \(p_i=(x_i,y_i)\) must satisfy

\[
\sum_iR_i=W,\qquad \sum_i x_iR_i=Q_x,\qquad \sum_i y_iR_i=Q_y.
\]

Dividing by \(W\) proves that a solution exists exactly when

\[
c=(Q_x/W,Q_y/W)\in\operatorname{conv}\{p_i\}.
\]

For every counterclockwise hull edge \(a\to b\), this is equivalent to

\[
g_{ab}=(b_x-a_x)(Q_y-Wa_y)-(b_y-a_y)(Q_x-Wa_x)\ge0.
\]

All inputs and hull predicates use rational arithmetic. At least one negative margin means no nonnegative vertical-reaction solution exists in the model. All positive margins place the resultant strictly inside the hull. Zero margins mark a boundary condition. A triangulation/barycentric calculation supplies an explicit nonnegative reaction witness supported on at most three pads.

This witness is an existence proof. Four real pads have a stiffness/contact problem that determines their actual forces; the witness does not solve it. The edge margins have dimensions force × length² and must not be reported as a dimensionless safety factor. Horizontal reaction distribution, friction, torsion, anchorage, finite pad area, ground bearing and settlement are absent.

## Entire-path checks

The supported path model keeps every weight constant and specifies its horizontal coordinates by polynomials of degree at most two in \(s\in[0,1]\). Horizontal forces are constant within a case and their application heights are also quadratic. Hence \(W\) is constant, the first moments are quadratic, and each support-edge margin has the form

\[
g(s)=a_0+a_1s+a_2s^2.
\]

Its minimum on the closed interval occurs at an endpoint, or at \(-a_1/(2a_2)\) when \(a_2>0\) and that point is internal. Every candidate location and value is rational. The implementation checks precisely this finite set, so it covers the continuous parameter interval for the given polynomial model; it is not a dense-grid approximation.

The recorded counterexample uses a unit weight with

\[
x(s)=\tfrac34+2s-2s^2,\qquad y(s)=0,
\]

over the square support hull \([-1,1]^2\). Both endpoints have \(x=3/4\), while \(x(1/2)=5/4\) lies outside. Endpoint-only checking misses this failure.

In the support challenge, each edge margin is affine in the two independently bounded horizontal-force components. At a fixed \(s\), the minimum over their rectangular domain therefore lies at a corner. Checking all four corners over the entire path proves the stated envelope result. It does not convert the assumed force bounds into a validated wind speed or account for inertial loads.

The discrete search retains all 198 enumerated combinations and their four path receipts. Its six Pareto points are nondominated with respect to ballast weight and support-rectangle area **within that finite grid**. There is no claim of a globally optimal physical mechanism or a continuous design-space optimum.

## Finite connector transforms

For an exactly unit quaternion \(q=(w,x,y,z)\), the standard quaternion rotation matrix has rational entries when \(q\) is rational. The code requires exact unit norm and does not silently normalize a supplied approximation.

For connector position \(p_i\) and relative translation \(t\), the finite relative displacement is

\[
\delta_i=t+(R(q)-I)p_i.
\]

With declared axis-aligned capture half-widths \(b_j\), each geometric margin is \(b_j-|\delta_{ij}|\). Checking every connector and axis handles translation plus finite rotation at the specified points. The independent test uses Hamilton products \(q(0,p_i)q^*\), not the implementation's matrix formula.

This is a geometric envelope test. It does not solve engagement forces, guide contact, elastic deformation, actuator force, interference away from the modeled points, load transfer or physical lock verification. The capture dimensions in the recorded examples are synthetic inputs, not manufacturing tolerances selected for this bridge.

## Reproduce

```bash
python -S verify_engineering.py
python -S verify_engineering_independent.py
python verify_engineering_independent.py --symbolic
```

The independent receipt checks all 198 design entries, 792 wind-corner paths, 3,168 edge minima, support sign controls, exact quadratic extrema, finite rotations, ramp boundaries and modified-report rejection cases. GIS source uncertainty remains outside these exact algebraic proofs.
