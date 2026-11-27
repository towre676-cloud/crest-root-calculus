# A uniform abscissa bound without Mordell–Weil generators

Let K be a number field containing i, let E be Y² = X³−1 and let C be y² = x⁶−1. Write h for the absolute logarithmic projective Weil height and H_abs=exp(h). On E write λ(P)=2ĥ(P), the canonical quadratic height normalized by λ(P)=lim 4^(−n) h(X([2ⁿ]P)). Set X(O)=∞.

**Theorem.** If rank E(K)≤1, then every finite abscissa on C(K) satisfies

H_abs(x) ≤ (36·156816^(1/3))^(1/4) = 6.6378100886… .

The same bound holds for an individual point whose two images below are linearly dependent modulo torsion, irrespective of the rank of E(K). For a field omitting i one may apply the theorem in K(i) when its rank hypothesis is certified and restrict the result to K. The rank condition can be certified by a 2-isogeny Selmer upper bound, by decomposition into rational quadratic character spaces in a multiquadratic input, or by another independent rank proof.

## Canonical height on the fixed elliptic target

In dehomogenized coordinates the duplication morphism on the abscissa has F=X⁴+8X and G=4X³−4. Homogenize both to degree four. The two polynomial identities

9X⁷=(9X³−8)F−16XG, 36=4X²F−(X³+9)G

give homogeneous identities for 9X⁷ and 36Z⁷ by homogenization to degree seven. To obtain common coefficient 36, scale the first identity by four. At a finite place v the coefficients are integral. In the sup norm this gives ||F,G||ᵥ≥|36|ᵥ ||X,Z||ᵥ⁴. At an infinite place the sum of absolute values of all the coefficients in the first certificate is 132 and in the second it is 14. Consequently ||F,G||ᵥ≥(36/132)||X,Z||ᵥ⁴. The forward archimedean coefficient bound is 9 and the finite-place forward bound is 1. Sum with the normalized place weights and apply the product formula. This proves

−log(132) ≤ h(X(2P))−4h(X(P)) ≤ log(9).

Telescoping after dividing by powers of four gives −A ≤ λ(P)−h(X(P)) ≤ B, where A=log(132)/3 and B=log(9)/3. These estimates also hold when an abscissa is infinity.

## Two rank-one images

For an affine point with x≠0 take P=(x²,y), Q=(x⁻²,iy/x³) on E. Their naive abscissa heights both equal 2h(x). If their free components span at most one real dimension, the quadratic form gives the exact equality

λ(P+Q)λ(P−Q) = (λ(P)−λ(Q))².

It follows that for at least one sign, R=P±Q satisfies λ(R)≤A+B. Thus h(X(R))≤2A+B=log(156816)/3. The inequality includes torsion, the identity, and the equal-canonical-height case.

Let j=i for the plus sign and j=−i for the minus sign. The chord law gives X(R)=N(x)/D(x), where

N=−2j x³−3x²+2jx, D=(x²−1)(x−j)².

These degree-four homogeneous forms have no common zero, as seen in the following Bézout identities, written with Z=1:

3=(4jx³+5x²−4jx−6)N+(−8x²+6jx+3)D,

3x⁷=(−6x³+4jx²+5x−4j)N+(3x³−6jx²−8x)D.

Both coefficient sums equal 36. Homogenize to degree seven; the same place-by-place argument gives 4h(x)≤h(X(R))+log(36). Combining the bounds proves the theorem. The rational abscissa identity extends to points with a vertical tangent or pole, by equality of morphisms to the projective line. The case x=0 has height zero.

The degree-four maps here are sums of the two displayed elliptic maps. They have j=0 target. The proof does not require them to be optimal covers.

## Effective finite enumeration

For arbitrary input K of degree d satisfying the rank hypothesis, a wholly elementary terminating enumeration is as follows. If an abscissa has primitive irreducible integer minimal polynomial of degree n≤d, its Mahler measure equals H_abs(x)^n. Every coefficient at exponent j then has absolute value ≤ binom(n,j)·7^n. Enumerate the finitely many primitive integer polynomials in these boxes, ask for all their roots in K, take exact square roots of x⁶−1 in K and restore the two infinities. Testing that actual roots lie below the bound is unnecessary, since every accepted point satisfies the source equation and every genuine point occurs. The large boxes affect speed. The optional script rank_one_general.py implements this method and a smaller-coefficient Gaussian version for fields containing i.

For the practical biquadratic replay K=ℚ(i,√d), classify points with x outside ℚ(i) from the previous exact quadratic-divisor enumeration on C/ℚ(i). Only d=2,3,7 support such points. For Gaussian x choose x=α/β in coprime Gaussian integers. The multiplicative height is max(|α|,|β|), so both norms are at most 44. The replay examines all 3,165 distinct ratios in this bound. It filters by the squareclass of α⁶−β⁶ and verifies each resulting ordinate in the input field. In detail, if z=a+bi≠0 is a positive rational number times a Gaussian square, t²=a²+b² is an integer square. If t+a>0, the required rational squareclass is the squarefree part of 2(t+a), by (Re u)²=(t+a)/(2d) for z=du². If t+a=0 it is the squarefree part of |a|. Distinct positive squarefree d give distinct rational squareclasses in ℚ(i). Checking the final source identity handles the converse.

For the 121 positive squarefree radicands 2≤d≤200, the rational character decompositions prove ranks exactly one in 34 fields and zero in 30. These 64 point sets are complete. Sixty-two sets consist precisely of the four base points and (0,±i). The two others are d=2, giving 18 points, and d=3, giving 10. Remaining tests have status partial. In particular no inference about positive-rank fields beyond the theorem is made from the finite Gaussian search.

**Consequence for an independent-point obstruction.** On a field where the available finite searches have exhausted the Gaussian height box and all non-Gaussian points, every further point must yield independent P,Q in E(K) modulo torsion. For the earlier field ℚ(i,√7), every additional abscissa would lie in ℚ(i) and the previous axis-aligned coefficient search through 65 covers this height box. Indeed rationalizing α/β gives denominator at most 44 and real and imaginary numerators at most 44.

## Rank-two CM extension

If L contains i and a primitive cube root ω and rank E(L)≤2, its free part tensored with the reals, when nonzero, is one complex line with the canonical-height norm. Multiplication by ω is an isometry of order three without a nonzero fixed vector because ω²+ω+1=0 on E. Thus the six unit phases rotate through multiples of π/3. For some u=±ωʲ the angle between P and uQ is at most π/6. When 2h(x)≥A, the lower bound is

λ(P+uQ)≥(2+√3)(2h(x)−A).

The abscissa of P+uQ equals a unit times the same reduced quartic function evaluated at a unit times x. Indeed under x↦ω^(−j)x, the two target abscissas scale by ω^(−2j) and ω^(2j), giving the relative rotation ωʲ after normalization. The numerator coefficient sum is 7 and the denominator coefficient sum is 8. Therefore λ(P+uQ)≤4h(x)+log(8)+B. This proves

h(x) ≤ ((2+√3)A+log(8)+B)/(2√3).

The same upper bound holds when 2h(x)<A. Its exponential is less than 14. The sum-map height idea is related to the classical Demjanenko–Manin method. These notes make no priority claim for a rank-restricted effectivity principle.

## Replay

```sh
python rank_one_height.py --panel
python rank_one_height.py 2
```

The first command also replays the exact 36-class Jacobian calculation over ℚ(i) underlying the degree-two abscissa classification.

## A non-Galois positive-rank example

Let K=ℚ(w), w⁴−2w²+2=0, with i=w²−1. Here K=ℚ(i,√(1+i)). The extension is non-Galois over ℚ. Normality would also put √(1−i) in K, forcing √2∈K. This would identify K with ℚ(i,√2), although both 1+i and (1+i)/2 have odd (1+i)-adic valuation in ℚ(i) and neither can be a square there, ruling out that identification.

The Selmer upper bound on E(K) is one and the full torsion order is six. The point (1−3w²/2,−3w/2) lies on E and its sixfold multiple is nonzero. Thus rank E(K)=1. The coefficient-field classifier gives exactly six C(K) points, namely the two infinities, (±1,0) and (0,±i). This computation uses the defining polynomial for K and the uniform height bound. It gives a positive-rank example outside every multiquadratic extension. Run `python non_galois_rank_one.py` to replay it.

For related general results see E. Viada, *An Explicit Manin-Dem’janenko Theorem in Elliptic Curves*, Canadian Journal of Mathematics 70 (2018), 1173–1200. https://doi.org/10.4153/CJM-2017-045-0

## Improved bounds for both rank hypotheses

The preceding simpler bounds can be tightened simultaneously. Write t=2h(x), C₀=log(36), A=log(132)/3 and B=log(9)/3. Use c=2 in the collinear case and c=√3 when K contains ω and the target has rank at most two. Then H_abs(x) is at most

exp(½ max{(A+C₀+(2−c)B)/c, (A−B+√((A+B)²+4((B+C₀)/c)²))/2}).

The corresponding upper bounds are **4.278300842** for rank at most one and **4.964624902** for CM rank at most two.

Indeed choose the sign or the unit phase giving the smaller sum. For a=λ(P), b=λ(Q) one can arrange λ(P+uQ)≤a+b−c√(ab), since the phase in the CM case makes the angle at least 150 degrees. When one vector vanishes this remains true. Both a and b belong to [L,U]=[t−A,t+B]. If t<A the asserted bounds follow immediately. Otherwise the function a+b−c√(ab) is convex separately in each variable on this rectangle. Its maximum occurs at a corner. The (L,L) value is bounded by (U,U), leaving only (U,U) and (L,U).

On the other hand, the inverse quartic height inequality and the canonical comparison give 2t−A−C₀≤λ(P+uQ). The rotated map is the same quartic up to unit scalings, so this holds for every selected phase. If (U,U) supplies the upper bound, then ct≤A+C₀+(2−c)B. If (L,U) supplies it, cancellation gives c√((t−A)(t+B))≤B+C₀. Solving the quadratic inequality gives the displayed expression. The interval checks in height_constants.py certify the decimal upper bounds. The previously reported panels used the larger initial box, which still proves the same complete lists.


The executable option `rank_one_general.py --bounded-only '[1,0,0,0,1]'` runs the finite point enumeration without a class-group or rank computation. It labels the result conditional on the rank hypothesis, allowing an external rank certificate to supply that input. The default command first obtains its own certified rank bound.
