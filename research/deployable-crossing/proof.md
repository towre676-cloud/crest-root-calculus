# Exact screens and their proof obligations

This note proves the algebra used by the first deployable-crossing screens. It does not validate a physical structure. All distances must use one consistent length unit and all forces one consistent force unit. Moments use their product. Positive vertical reactions act upward; positive external moment is clockwise in the stated planar coordinate system.

## 1. Uniform horizontal cantilever

Let $w\ge0$ be uniform downward force per length and $\ell\ge0$ horizontal reach. Integrating the load gives $V=\int_0^\ell w\,dx=w\ell$ and $M=\int_0^\ell wx\,dx=w\ell^2/2$. Both expressions are monotone in nonnegative $w,\ell$, so independent closed intervals give exact extremal bounds by their endpoints. Extending toward negative $x$ reverses the signed external gravity moment while preserving the shear magnitude.

The formula is about a load distributed along the current reach. The program does not infer the moving mass distribution of a telescoping, folding, or scissor mechanism. Stored material and its changing center of gravity need their own load description. An elastic cantilever's end displacement also depends on bending rigidity and boundary conditions, which this screen does not supply.

For a horizontal force $F$ acting toward positive $x$ at height $z$ above the moment origin, the clockwise moment is $Fz$. The function encloses this signed product. A negative force reverses its sign. Wind force is an input, not a value obtained from a wind-speed model.

## 2. Two unilateral vertical supports

Let $R_L,R_R$ be upward reactions at fixed horizontal locations $x_L<x_R$. Force and moment equilibrium give

$$R_L+R_R=W,\qquad x_LR_L+x_RR_R=M.$$

The unique solution is

$$R_L=\frac{x_RW-M}{x_R-x_L},\qquad
R_R=\frac{M-x_LW}{x_R-x_L}.$$

Since the denominator is positive, both reactions are positive exactly when $x_LW<M<x_RW$. The reaction equations describe a unilateral-pad equilibrium only while nonnegative reactions are possible. A negative algebraic reaction is an indication that this support model cannot realize the required equilibrium. It is not a verified load in an unspecified anchor.

For independent intervals in $W,M$ and fixed support positions, each numerator is affine and its extrema occur at endpoints. Interval multiplication by a known signed coordinate therefore encloses the full box without sampling. At an exact zero reaction the screen reports a boundary. If an interval includes both favorable and unfavorable states, its classification must not assert a universal pass. Intervals enclosing a correlated physical set may be conservative; they are not a completeness theorem for that smaller set.

The input $M$ is the **total external** clockwise moment about the declared origin, including every downward load at its horizontal coordinate and any applied couples or horizontal forces at height. The local root moment of a bridge half is not an extra external moment to add after its distributed weight has already been included for the complete unit.

No-lift is only one necessary screen. A three-dimensional support polygon, lateral force equilibrium, friction, soil pressure, settlement, actuator strength, structural capacity, and prescribed design margins are different calculations.

## 3. Rigid-half constraint rank

Represent each planar rigid half by infinitesimal base translation $(u_x,u_y)$ and rotation $\theta$. A point at relative coordinates $(r_x,r_y)$ has velocity

$$v=(u_x-\theta r_y,\;u_y+\theta r_x).$$

A base pin imposes $u_x=u_y=0$; a fixed base also imposes $\theta=0$. A vertical roller restrains vertical base translation while permitting horizontal translation and rotation. A center hinge imposes equality of the two end-point velocities. A rigid center additionally imposes equality of their rotations. The rank of these homogeneous rational constraints is computed by exact Gaussian elimination; its nullity is the number of first-order freedoms of this ideal model.

For two horizontal halves, with relative tips $(a,0)$ and $(-b,0)$ and $a,b>0$, base pins reduce the center-hinge condition to

$$a\theta_L+b\theta_R=0.$$

The horizontal velocity condition is then automatically satisfied, leaving one first-order freedom. A rigid center adds $\theta_L=\theta_R$, so $(a+b)\theta_L=0$ and both rotations vanish. With fixed bases, both rotations vanish before a center joint is introduced. The vertical-roller example with a left pin has the same distinction after its center horizontal-velocity condition restrains the remaining base translation.

This is a kinematic rank statement, not a member-stiffness or buckling calculation. It assumes the endpoints are geometrically in contact; compatibility of their actual coordinates is a separate check. A first-order freedom at a straight singular configuration need not integrate into a finite motion with all geometric distances held fixed. Conversely, full rank here does not prove adequate elastic stiffness or strength.

## 4. Scalar joining envelope

Let $S$ be the horizontal separation between the two declared reference origins, and $r_L,r_R$ their inward reaches. Define the scalar overlap

$$o=r_L+r_R-S.$$

Positive $o$ means nominal overlap under this convention; negative $o$ means a gap. Given explicitly declared allowed overlap $[o_{\min},o_{\max}]$ and maximum height mismatch $d_{\max}$, the tested inequalities are

$$o-o_{\min}\ge0,\qquad o_{\max}-o\ge0,\qquad
d_{\max}-|h_L-h_R|\ge0.$$

An interval is fully within the stated geometry window only if its complete overlap and mismatch enclosures satisfy these limits. Boundary cases are distinguished from strict interior clearance. These inequalities do not design a connector or determine whether insertion force, angular mismatch, lateral offset, pin engagement, or material stress is acceptable. The overlap inputs are an abstract coordinate contract, not instructions to collide two solid bridge noses.

## 5. Ramp arithmetic

For nonnegative rise $H$, positive horizontal run $r$, and positive declared maximum slope $m$, the slope screen is $H/r\le m$, equivalently $H\le mr$. Exact fractions preserve equality at a limit. This comparison alone does not establish accessible width, landings, edge protection, handrails, cross slope, surface quality, or code applicability.

The example uses $H=18$ ft, $m=1/12$, maximum rise per run $h_{\max}=5/2$ ft, and intermediate landing length $b=5$ ft. Therefore $r=216$ ft, $n=\lceil H/h_{\max}\rceil=8$, and the aggregate intermediate landing length is $(n-1)b=35$ ft. The sum $251$ ft is horizontal route travel before the terminal landings, not a minimum straight parcel length. The actual inclined travel on the ramp segments is greater than their horizontal run.

## 6. Logical records and physical evidence

The implemented record checker recognizes prefixes of one fixed sequence of declared state labels. It rejects known labels that are skipped, repeated, or out of order, and rejects unknown labels as malformed input. It does not verify evidence, handle emergency transitions, or infer that a support is engaged, a pin is carrying load, a roadway is clear, or an instrument is accurate. This small example must not be connected to actuators or used as a safety controller. A physical system needs an engineered control architecture, appropriate fault behavior, independent validation of observations, and a justified operating procedure.

## 7. What a future whole-path certificate would need

A continuous interval enclosure can cover all states of a specified motion segment. To certify a full deployable crossing, each enclosure would need to cover the correct structural configuration, moving masses, support contacts, loads, tolerances, and collision geometry. Discrete engagement and release events require separate model changes. Every relevant constraint must be included; proving a subset of inequalities cannot establish the omitted ones. The current kernel proves only its named elementary screens.
