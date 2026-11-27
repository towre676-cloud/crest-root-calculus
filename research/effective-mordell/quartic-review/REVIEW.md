# Degree-four quotients and Effective Mordell

The four solver repairs from this review are incorporated in the [number-field reconstruction](../../elliptic-quotients/quartic-gluing-reconstruction/number-fields.md).


The submission reviewed is commit **093779c145d062dbbc1b0698f858f5b9f2da4a59**. The supplied Git bundle names **7303ec15aeeaa858e7f976d46e78430939cedff6** as its prerequisite. The complete Git tree computed from the project ZIP is 26dd4eaa1f6e8d6de387e5ec6ee71b512059ef73. It matches the tree named in the bundle commit. The evidence below concerns explicit cases of ProofAtlas #87. [Exact problem](https://www.proofatlas.ai/open-problems/)

## Stoll's external benchmark

Stoll's Example 5.3(4) is \(C:y^2=6x^5+28x^3+54x\). He gives it as a minimal-degree-four example on which an earlier Magma Chabauty procedure got stuck. Lemma 5.5 in the same paper develops another route for the rank-one case. The control run on Magma V2.29-10 returned both rational points and certified completeness. [Stoll, Example 5.3 and Lemma 5.5](https://mathe2.uni-bayreuth.de/stoll/papers/rank1.pdf)

Running Claude's inverse solver on the quintic gave two rational maps. In the original coordinates, their abscissas and their ordinates divided by \(y\) are as follows.

| Elliptic target | \(X\) | \(Y/y\) | \(\mathbb Q\)-rank |
|---|---|---|---|
| \(Y^2=X^3-108X\) | \(-32x/(x^2+1)^2\) | \(8(3x^2-1)/(x^2+1)^3\) | 0 |
| \(Y^2=X^3-12X\) | \(-32x^3/(x^2+9)^2\) | \(-8x(x^2-27)/(x^2+9)^3\) | 1 |

The first target has rational torsion of order two. Its only affine rational point is \((0,0)\), so any finite rational source point has \(x=0\). The horizontal denominator is everywhere nonzero on the rational affine line. The quintic model has one rational point at infinity. Thus \(C(\mathbb Q)=\{(0,0),\infty\}\). The files include the solver's original output in its reciprocal chart and the exact identity in the original coordinates.

For the quadratic twist \(C_d:y^2=d(6x^5+28x^3+54x)\), multiply the listed \(X\) and \(Y/y\) by \(d\). The targets become \(Y^2=X^3-108d^2X\) and \(Y^2=X^3-12d^2X\). Each of the two identities holds symbolically in \(d\). For every positive squarefree \(d\le30\), the verifier computes both rank intervals and torsion orders. For 17 of these 19 parameters at least one target has rank zero and torsion order two. The same proof yields exactly \((0,0)\) and \(\infty\). For \(d=14,22\) both ranks are one. A bounded search on \(d=22\) also finds \((1,\pm44)\). The rational point sets of these two members remain uncertified here. These are related twists of one published example. They measure when the arithmetic method applies on this family.

## A rank-two comparison for Magma

The most informative comparative input from the twist panel is **\(d=26\)**. The two pullback differentials are \(4\,dx/y\) and \(-4x\,dx/y\), so the induced map from the Jacobian to the product of the two targets is a rational isogeny. The respective rank intervals are \([0,0]\) and \([2,2]\), giving Jacobian rank two. The rank-two target has the independent points returned by PARI at \((-48,528)\) and \((208,2704)\). The other quotient determines the complete point set. Stoll's minimal rational cover degree four persists under this simultaneous quadratic twist. In detail, translate a hypothetical rational lower-degree map to send the rational Weierstrass point at infinity to the origin. The hyperelliptic involution then acts by negation on the map, by its induced Jacobian homomorphism. Twisting both source and target back by 26 would give a rational map of the same lower degree on the published curve, a contradiction. The documented Magma routine includes degree two and three subcovers and rank-at-most-one Chabauty, as well as descent and sieving. On 27 November 2026 the public Magma V2.29-10 calculator returned the same two points, a false completeness flag, and a certified search height of 20,000 for this twist. It returned a true completeness flag on the original curve. The comparison transcript is included. Crest–Root also recovered the two maps directly from the twisted source equation in a fresh run. Its scaled target coefficients were −4563 and −507, respectively. The comparison establishes an advantage in certification on this specific input. [Magma handbook](https://magma.maths.usyd.edu.au/magma/handbook/text/1613)

## The two repeated-factor curves

Both \(T_1,T_2\) map over \(\mathbb Q\) to \(E_-:Y^2=X^3-105903X-11165202\). PARI's lower and upper bounds on the rank are both zero. Its torsion group has order four, and the three affine points are \((-246,0),(-123,0),(369,0)\). [PARI rank and torsion functions](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html)

| Curve | Rational finite source abscissas to check | Poles | Leading coefficient |
|---|---|---|---:|
| \(T_1\) | \(-2/5,\ 7/3,\ 16\) | \(x^2-10=0\) | 154 |
| \(T_2\) | \(-11,\ -14/5,\ 8/3\) | \(x^2+2=0\) | 28 |

For each displayed abscissa, \(f(x)\) is a nonsquare in \(\mathbb Q\). Neither pole polynomial has a rational root, and neither leading coefficient is a rational square. This proves \(T_1(\mathbb Q)=T_2(\mathbb Q)=\varnothing\).

Both curves have real points and points over **every** \(\mathbb Q_p\). The sextic discriminants have prime support \(2,3,5,41\). For \(p=2,3,5,7,11,13,41\), respectively, the following integers \(x\) give a square value of \(f(x)\) in \(\mathbb Q_p\).

| \(T_1\) | 0 | 11 | 2 | 5 | 0 | 3 | 3 |
|---|---:|---:|---:|---:|---:|---:|---:|
| \(T_2\) | 3 | 1 | 0 | 2 | 1 | 1 | 0 |

For every remaining prime \(p\ge17\), good reduction and the Hasse–Weil lower bound \(p+1-4\sqrt p>0\) give a smooth point that lifts. The positive leading coefficients give real points. These are explicit failures of the Hasse principle detected using rational degree-four quotients. We have not established priority or compared these cases with modern full descent. [Hasse–Weil bound for hyperelliptic curves](https://www.mathe2.uni-bayreuth.de/stoll/teaching/ArithHypKurven-SS2014/Skript-ArithHypCurves-pub-print.pdf)

The geometric computation contributes different information. The extra maps have conjugate differentials and require \(\mathbb Q(i)\). Their displayed elliptic targets are defined over \(\mathbb Q\), and their rational ranks are 0 and 1 for the two signs of the constant term. The targets are quadratic twists of each other by \(-1\), so each has rank \(0+1=1\) over \(\mathbb Q(i)\), by the two eigenspaces for conjugation after tensoring with \(\mathbb Q\). Thus the finite-fibre computation above genuinely uses the rational rank-zero quotient. A separate example is needed to demonstrate an arithmetic gain from a map available only over an extension.

## Findings for the solver author

**A positive example for the missing pairing case.** On \(y^2=x^6+1\) the original solver reaches model scaling and fails in PARI with a valuation at \(p=-1\). In the diagnostic run the sole change was to use the absolute value of the rational content of the lifted coefficient polynomial before factoring it in normalize_scale. The replay then returned six quadratic map orbits, or twelve distinct maps. Four of the orbits have pairing field \(K=\mathbb Q(w)\), \(w^2-w+1=0\). For one of them,
\[
S=x^2+(w-2)x+2-2w,\qquad E:Y^2=X^3+15wX+22,
\quad \sigma(x)=((w-1)x^5-2wx^3-2x)/3 .
\]
The verifier checks the full map identity, the differential, all twelve distinct half-poles, and the three branch images of multiplicity two for the four selected representatives. In \(K[x]/(x^6+1)\), their abscissas are fixed by \(\sigma\) and fail to be fixed by its conjugate \(w\mapsto1-w\). Thus the pairing itself needs the quadratic field, a stronger property than the map needing that field. Full numerators and ordinates are in the input records. The same source also has the elementary rational double cover \((x,y)\mapsto(x^2,y)\) to \(Y^2=X^3+1\). The value of this test is the geometric recovery of a nonstable pairing.

1. **One-variable chart.** Apply the rational change of coordinates \(x=1+1/u,\ v=u^3y\) to the archived PERMUTED sextic. Its ascending coefficients become \([1,0,-26,-70,-57,-14,-1]\). The original has two maps and isomorphism transports both. In the submitted code the new equation raises an error about positive dimension in the chart with one pole parameter. The polynomial \(x^2-1\) of one original half-pole becomes \(2u+1\). In one variable the nonconstant saturated common divisor represents isolated points. Factor it, check residuals and degeneracies at its roots, and branch to the plane-component logic only for two variables.

2. **Rational input.** For the monic \(f(x)=\prod_{j=0}^{5}(x-j/2)\), the first separating weight is 1. The \(z^{14}\) coefficient of the exact resolvent is \(-418377/64\). Rounding this to an integer cannot recover that resolvent. Use a primitive integral multiple for the root invariants while retaining the original source equation for the cover identity. Include separate cases of rational coefficients and nonsquare leading coefficient.

3. **Distinct pairing orbits.** Exact involutions, generation of each individual field, and the sum 15 establish completeness once different accepted orbits are proved disjoint. Numerical membership in the resolvent factors does not supply that last exact check. For every tested main record and the Stoll curve, the weighted sum of the six involution coefficients with weights \((1,2,3,5,7,11)\) has degree equal to the proposed orbit size and a minimal polynomial coprime to those of all other orbits. That proves all 15 pairings are distinct on these inputs. A general implementation can try separating linear combinations with an exact acceptance criterion and refine if a combination fails. Alternatively, recover each pairing's original resolvent invariant exactly in its stated factor.

## Regression coverage

On the 36 Richelot neighbours the equation-only reconstruction and independent map and disjointness checks passed for all rows. There were 32 rows with two geometric maps and four with four. The rational half-poles and point fields matched the archived results. The six main cases and the rational Galois verifier were also replayed. The three further controls in review-inputs.json are Stoll’s original, its rank-two twist, and the nonstable-pairing curve after the content-sign correction.

## Reproduction

Run **python verify_review.py** with SymPy and cypari2 available. The script reads review-inputs.json, checks pairing involutions and cross-orbit disjointness, verifies every displayed cover identity and differential, recomputes the ranks and torsion used in the arithmetic, checks every fibre and verifies local solubility. It also verifies both identities for the twist parameter symbolically and recomputes the 19 rank pairs. It uses PARI for exact arithmetic and its rank bounds. It is an independent check of the outputs and the rational-point argument. Discovery can be replayed separately from Claude's original ZIP with verify_number_fields.py, its --sweep option, and verify_galois.py.


For the extended panel, run **python verify_neighbours.py SOLVER_DIRECTORY**, with SOLVER_DIRECTORY pointing to research/elliptic-quotients/quartic-gluing-reconstruction in the extracted source. The original dependencies include python-flint. The auxiliary Magma script includes T1, T2, and d=22 in addition to the two public cases actually run. The empirical comparison is documented in magma-comparison.txt.
