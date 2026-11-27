# Exact local 2-isogeny descent over arbitrary number fields for the fixed Mordell curve

For C given by y² = x⁶ − 1, the direct elliptic targets are E₋, E₊ with equations Y² = X³ − 1 and Y² = X³ + 1. The algorithm computes the two isogeny Selmer sizes for each target over any input number field K. Its output is a certified rank upper bound for each target. A zero bound makes the C(K) point list exhaustive by enumerating torsion and pulling back the abscissas. Positive bounds keep the point-list problem open in this procedure.

## Global squareclasses

Translate the rational root −sign to zero. The two isogenous models are Y² = X(X² + aX + b) with parameters (a,b) = (−3sign,3), (6sign,−3). An element e of a descent image is represented by the projective quartic Z² = e U⁴ + a U²V² + (b/e)V⁴. Its valuations outside S, the prime ideals over 3, must all be even. Indeed a point of the elliptic curve at a prime outside S has abscissa of even valuation if this valuation is nonzero, by inspecting the defining equation. For the local Selmer images this same argument applies prime by prime.

Write (e) = P A² with P the product of the primes over 3 at which the valuation is odd. The ideal A may be fractional. In the certified ordinary ideal class group, compute the coordinates v of [P] in a cyclic decomposition with basis gⱼ of orders hⱼ. The allowable coordinates z of [A] are exactly all solutions of 2zⱼ ≡ −vⱼ modulo hⱼ. Choose the corresponding products A of ideal basis generators. Principalize P A² and check that the returned generator defines that very ideal. Multiply its generator by one representative of every unit modulo squares, using the torsion generator and fundamental units. This enumerates every allowable squareclass exactly once. Principal changes in A multiply e by a square times a unit, while distinct parity vectors or distinct choices of [A] cannot become equal modulo squares. The class and unit computations use PARI's unconditional bnf certification.

## Exact local squareclasses

For a bad prime ideal above 2 or 3, choose a uniformizer π. Replace e by r = e / π^(2 floor(v(e)/2)). Both r and b/r are locally integral. This replacement preserves local solubility, as follows by scaling U and Z. For a locally integral unit with global denominators, compute a globally integral lift congruent modulo π^(1+2v(2)). Verify that congruence exactly. Hensel's lemma shows that the ratio is a local square. Equal local squareclasses can now share a single test.

The two primitive coordinate charts for the quartic use t in the valuation ring and t in its maximal ideal, respectively, interchanging the leading and constant coefficients. Let F(t)=r t⁴ + a t² + b/r. Enumerate residue disks about c with depth n. If F(c) is square, there is a local point. If v(F(c)) > 2v(F′(c)), Hensel gives a local root and hence a point with ordinate zero. For the other disks, Taylor expansion bounds v(F(c+h)−F(c)) from below, for v(h)≥n, by

min{n+v(F′(c)), 2n+v(6rc²+a), 3n+v(4rc), 4n+v(r)}.

If this exceeds v(F(c)) + 2v(2), every value on the disk has the same nonsquare class and the disk is discarded. Otherwise split into all residue subdisks. There is no imposed precision cutoff.

Termination follows by compactness. At a nonroot in a nested chain the value is either a nonzero square or a nonsquare. In the first case nearby centers eventually have the same square class and produce a witness. In the second case the Taylor bound eventually rejects the whole disk. At a root the derivative is nonzero, since b(a²−4b)≠0, and the Hensel test eventually succeeds. A finitely branching infinite unresolved tree would contain an infinite nested chain. Thus the local test terminates. The remaining enumeration is finite once certified class group and unit data are supplied; those data are effectively computable by classical number-field algorithms. PARI documents that bnfcertify itself may loop if provisional BNF data are incorrect. The supplied software requires successful certification before using any result, and obtained it for every field in the receipts.

At real places the first isogeny requires e positive and this suffices by evaluating at an endpoint of the projective quartic. For b<0 one of the endpoints has positive value for either sign of e. At a prime away from 2 and 3, valuation normalization gives a smooth genus-one reduction with a point by the finite-field Hasse bound, and that point lifts. Complex places impose no restriction. The survivors are therefore precisely the two isogeny Selmer groups. For their cardinalities A and B, the descent identity gives rank E(K) ≤ log₂(AB) − 2. The default squareclass enumeration is exhaustive. An optional finite global budget yields a partial result when exceeded. Rank-zero can also occur with a positive Selmer upper bound.

## Computed coverage

On the previous 145 higher-degree models the completed point lists rise from 23 to 52, with 32 of 75 cubic, 7 of 35 quartic and 13 of 35 quintic fields certified. An additional sample chosen to have nontrivial ideal class group certifies 46 of 200 cubic, 8 of 100 quartic and 17 of 100 quintic fields. The six repeated field polynomials across these two panels all remain partial. The union comprises 539 distinct defining polynomials with 123 complete classifications, of which 122 have exactly four points and one has eight. These panels illustrate coverage and have no statistical interpretation. All 123 completed lists include exact coordinates. All 539 entries include four exact local Selmer cardinalities.

For a concrete class-number-four example, let w³ − w² − 8w − 10 = 0. Both primes above 3 are nonprincipal, with classes 3 and 1 in ℤ/4. The two Selmer sizes for each of the two elliptic targets are (2,2). Thus C(K) consists of the two points (±1,0) and the two points at infinity. A separate quadratic regression compares the computed upper bounds with independently certified lower ranks of rational quadratic twists in all 243 signed squarefree fields with |d|≤200.

## Replay

The directory contains the inherited elliptic torsion and pullback routine, the new solver, the full per-field receipts, and their verifier. Install cypari2 with PARI/GP. The calculations used cypari2 2.2.4 and PARI 2.17.2.

```sh
python classgroup_descent.py '[-10,-8,-1,1]'
python verify_classgroup_descent.py --quadratic-check
```

The coefficient list is ascending. The verifier optionally accepts a positive worker count. Torsion abscissas are checked by substitution into the original genus-two equation.

## References

PARI Group, *Catalogue of functions, General number fields*, entries for bnfinit, bnfcertify, bnfisprincipal, idealprimedec, idealchinese and nfislocalpower. https://pari.math.u-bordeaux.fr/dochtml/html/General_number_fields.html

J. H. Silverman and J. T. Tate, *Rational Points on Elliptic Curves*, on 2-isogeny descent and homogeneous quartics. The global Selmer condition above uses the corresponding construction over a number field.
