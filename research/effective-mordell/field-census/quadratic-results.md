# The full quadratic field window and an infinite family of prime radicands

Let C be the smooth projective model of y² = x⁶ − 1 and S its two rational points at infinity and its two affine points (±1,0). For every signed squarefree integer d with 0 < |d| ≤ 200 and d ≠ 1 the point set C(ℚ(√d)) is now complete. There are 243 such fields. In 232 the list is S, in one it has six points, and in ten it has eight. The six fields among the last 30 that acquire points have radicands ±19, ±95, ±182. The other four eight-point fields are the previously certified d = −3, −2, −7, 7. The six-point field is ℚ(i).

For d = 19, 95, 182 let (m,n,h) be respectively (5,3,28), (8,7,39), (3,1,2). Over ℚ(√d) the extra points are (±m/n, ±h√d/n³). Over ℚ(√−d) they are (±n/m, ±h√−d/m³). Each pair of signs varies independently. These formulas follow from m⁶−n⁶ = d h².

## Rational descent

The earlier exhaustive rational divisor calculation proves that every nonrational quadratic abscissa on C is over ℚ(√−2) or ℚ(√−3). For each remaining field one needs only rational x on z² = d(x⁶−1). Taking x to 1/x interchanges the positive and negative twists for x≠0. The zero and infinite abscissas are handled in the earlier classification. Hence it suffices to take d > 0.

Write Q₊ = x²+x+1 and Q₋ = x²−x+1. Let a,b be their positive squarefree representatives. They divide the odd radical of 3d as proved in the first descent. The earlier local and rank-zero tests leave exactly two reflected ordered pairs in every outstanding positive radicand. Use x↦−x to put (a,b) in the orientation in the table below. Set k = sf(dab), where sf means the squarefree part. A nonbranch point satisfies x²−1 = k r². Set s = (x−1)/r and t = s². Then

x = (k+t)/(k−t), s⁴+3k² = a A², 3s⁴+k² = b B².

These identities result from substituting in Q₊ and Q₋ and clearing (k−t)². A nonbranch rational solution has t finite and positive, with t ≠ k. Conversely the two quartic square conditions imply all three original squareclasses and therefore the rational twist condition.

## Eisenstein squares

Write s = U/V in coprime rational integers and F = ℚ(π), where π² = −3. Use the integral basis (1,ω) with π = 1+2ω and ω²+ω+1 = 0. Let β = U² + kπV² and γ = kV² + πU². Their norms are a and b times integer squares. The ring ℤ[ω] is Euclidean, with unit squareclasses represented by ±1. A prime of F dividing both β and its conjugate, or both γ and its conjugate, lies over 2, 3 or a divisor of k, as follows by subtracting and adding and using gcd(U,V)=1. Consequently outside 6kn, where n is the prescribed norm squareclass of the factor, every valuation is even. This holds at inert primes as well because any divisor there divides the conjugate.

A squareclass e with eβ or eγ square is enumerated as follows. Choose a unit sign. The exponent of π over 3 has parity v₃(n). A prime q inert in F cannot divide n, and for q dividing 6k its generator can occur with exponent zero or one. At a split q dividing n choose one of the two conjugate prime generators. At a split q dividing 6k but not n choose either neither or both. These choices include every permitted class by the norm parity and the ideal divisibility statement. In the code each representative is an algebraic integer, so testing a square modulo an integer modulus is a valid necessary condition when its valued product is a square in F.

For each candidate pair (δ,η) require δβ, ηγ to be squares in ℤ[ω] modulo each successive modulus 16, 9, 7, 13, 19, 97. Simultaneously test a Nβ and b Nγ as integral square residues. Representatives (U,1) and (1,V), with the underlying rational prime dividing V in the second chart, exhaust primitive rational projective pairs. Scaling a pair by a rational unit scales β and γ by its square. The table lists the number of candidates initially and after each of the six moduli in this order.

| d | a | b | k | Counts |
|---:|---:|---:|---:|---|
| 19 | 1 | 19 | 1 | 32,2,1,1,1,1,1 |
| 31 | 3 | 31 | 3 | 32,2,1,0,0,0,0 |
| 43 | 1 | 43 | 1 | 32,2,1,0,0,0,0 |
| 67 | 1 | 67 | 1 | 32,2,1,0,0,0,0 |
| 79 | 3 | 79 | 3 | 32,2,1,0,0,0,0 |
| 95 | 1 | 57 | 15 | 128,8,2,1,1,1,1 |
| 103 | 3 | 103 | 3 | 32,2,1,0,0,0,0 |
| 127 | 3 | 127 | 3 | 32,2,1,0,0,0,0 |
| 129 | 1 | 129 | 1 | 32,2,1,0,0,0,0 |
| 139 | 1 | 139 | 1 | 32,2,1,1,0,0,0 |
| 151 | 3 | 151 | 3 | 32,2,1,1,1,0,0 |
| 163 | 1 | 163 | 1 | 32,2,1,0,0,0,0 |
| 182 | 7 | 13 | 2 | 64,8,2,2,2,2,1 |
| 193 | 1 | 193 | 1 | 32,2,2,0,0,0,0 |
| 199 | 3 | 199 | 3 | 32,2,1,1,0,0,0 |

## The three nonempty cases

For d = 19, 95, 182 respectively the remaining pairs (δ,η) are (1,−3−5ω), (−2,16+2ω), (−2−6ω,6+8ω). Over F the factors δ(t+kπ), ηπ(t+k/π) and t are squares. In the first two cases put h = −kδη, c = −πh/k. Then X=ct gives a point on Y² = X(X+3h)(X−h). Scaling h and c by the rational square 1/4 gives a smaller equivalent model in the case d=95.

For d=182 use instead the square factors t, conjugate(δ)(t−kπ), and ηπ(t+k/π). Their product gives a point on Y² = X(X−2cπ)(X+2c/π), where c = conjugate(δ)πη and X = ct. The other first-choice model in this case has rank two.

The supplied final_three.m proves rank bounds [1,1] and obtains a certified odd-index Mordell–Weil subgroup on each chosen model over F. Magma's rank-one elliptic Chabauty result for the rational projection X/c has index obstruction 8,4,8 respectively. Those integers are powers of two, hence coprime to the subgroup indices. With full Minkowski class-group proof bounds, the complete rational projection sets are

| d | t-values |
|---:|---|
|19|∞, 0, 1/4, 4|
|95|∞, 0, 1, 225|
|182|∞, 0, −1, 4|

The single combined script was executed in Magma V2.29-10, seed 4120333681, and all assertions passed. Excluding zero, infinity, negative t and inconsistent ordered squareclasses gives x=5/3, 8/7, −3 respectively. Reflection and the reciprocal twist give exactly the six additional lists asserted above. The squareclass arithmetic is replayed in verify_new_descent.py; the three assertions involving elliptic Chabauty are replayed using Magma.

## A uniform ray-class theorem

Set M = 64·81·7 = 36288. Let p be any rational prime admitting an Eisenstein integer g of norm p with g ≡ 1+7ω mod Mℤ[ω]. Then both C(ℚ(√p)) and C(ℚ(√−p)) equal S.

In fact p ≡ 43 mod M. For a,b dividing 3p the primitive square tests on aQ₊, bQ₋ and pab(x²−1), homogenized as binary quadratics, modulo 64 and 81 therefore exactly reproduce the p=43 tests, using the labels p and 3p in place of 43 and 129. There are just three unordered survivor types: (1,3), (p,3p), and (1,p). In the first two, sf(ab)=3 yields the fixed elliptic curve Y²=X(X²+3X+9), with X=3x². Its rank zero and all four torsion points give only x=0,±1 in the possible finite fibres. For the third, k=1. Norm-p representatives use g or its conjugate, with the same signs and possible inert factor 2 as in the case p=43. Their residues at 16, 9, 7 coincide with those for 1+7ω and its conjugate because these moduli divide M. The candidate counts are therefore identically 32, 2, 1, 0. This proves the result for every qualifying prime.

There are infinitely many such p. The prime ideals in the ray class of (1+7ω) modulo M have positive Dirichlet density by the ray class field and Chebotarev. The degree greater than one prime ideals of this imaginary quadratic field have zero Dirichlet density, so infinitely many members have degree one. For each such ideal equality of ray classes supplies a generator congruent to 1+7ω; its norm is a rational prime. The generator-norm criterion is decidable with Cornacchia's algorithm and the twelve associate/conjugate choices. Apart from p=43, explicit examples are 1,317,109,291 and 1,317,290,731, with generators (36289,36295), (1,36295) respectively in the basis (1,ω). All p here are 7 modulo 12 and lie outside the earlier restricted multiquadratic tower.

Reference for the prime distribution statement: Alessandro Cobbe, *Steinitz classes of tamely ramified Galois extensions of algebraic number fields*, Chapter 1, the ray-class corollary of Chebotarev’s theorem.

## A density theorem for prime radicands

There is an explicitly decidable family of rational primes of natural density **927/1024** such that C(ℚ(√p)) = C(ℚ(√−p)) = S. This gives a lower bound for the density in the sense that the set of primes with this point property contains a set of that natural density. It does not require that the entire set with the point property has a natural density.

First let d be any signed squarefree integer whose prime divisors are 2, 3 or primes equal to 2 modulo 3. A primitive numerator U² ± UV + V² is odd, and no odd prime q = 2 mod 3 can divide it: a root for the ratio U/V would give a nontrivial cube root of unity in the residue field (with sign changed in one polynomial), while neither U nor V can vanish in a divisible primitive form. Hence a,b belong to {1,3}, and sf(ab) is 1 or 3. If k′=sf(ab), the squareclass of Q₊Q₋ gives the rational elliptic curve Y²=X(X²+k′X+k′²), X=k′x². These two curves have rank zero over ℚ; the first has rational torsion of order 2, the second order 4. Their only possible finite square abscissas give x=0,±1. The earlier quadratic divisor enumeration handles d=−2,−3, and x=0 requires d=−1. This proves a full point formula for the entire multiplicative family and deals with all but finitely many prime radicands inert in F.

For the split primes p outside 2,3,5,7,13,17,19,97, write p=g·conjugate(g) in the Eisenstein integers. There are four classes 1,3,p,3p for each of a,b. If p occurs in both or in neither, sf(ab) is 1 or 3 and the same fixed elliptic curves exhaust the possible abscissas. Otherwise reflection puts the ordered pair in the form (a,pv) with a,v independently 1 or 3. The integer k=sf(ap²v) belongs to {1,3} and is independent of p. The two norm tests in each of these four named cases have 32 global labels. For δ, take one of {±1,±2} times π if a=3, and times 1 otherwise. For η, take one of {±1,±2} times either g or conjugate(g), and times π if v=3. Thus there are exactly 128 labels before local tests.

At each rational prime q in the support of the test, a label first requires the three rational homogeneous square conditions for (a,pv). It then requires the two quartic norm conditions and the two Eisenstein square conditions in the preceding descent. The rational moduli are 64 at 2, 81 at 3, and q at the other places. The Eisenstein moduli are 16 at 2, 9 at 3, and q elsewhere. For fixed global labels the local output depends only on the class of g in the units of ℤ[ω]⊗ℤq modulo squares. In fact replacing g by a local unit square times g multiplies the conjugate by another unit square, and multiplies p through the norm by a rational q-adic unit square. Every condition is invariant under these changes. Local unit squareclasses are obtained by enumerating units mod 8 at 2 and mod q elsewhere. At the dyadic place 1+8O consists of squares by Hensel, and at odd places square units lift from the residue fields, including at the ramified prime over 3.

For each such local squareclass, prime_density.py outputs the bit mask of the labels left possible. Intersection composes the masks across different rational primes. The Chinese remainder theorem gives every tuple of local squareclasses. The successive exact counts of all tuples and of tuples with empty intersection are as follows.

| Auxiliary prime included | Total tuples | Empty intersection |
|---:|---:|---:|
| 2 | 8 | 2 |
| 3 | 16 | 6 |
| 5 | 32 | 12 |
| 7 | 128 | 82 |
| 13 | 512 | 382 |
| 17 | 1024 | 774 |
| 19 | 4096 | 3224 |
| 97 | 16384 | 13280 |

To justify the distribution, take the ray modulus 8·3·5·7·13·17·19·97 of F. Its ideal ray class group is the product of residue-unit groups modulo the image of the global units, since F has class number one. Its quotient by squares is exactly the product of local unit squareclass groups modulo that same image. Ray class field theory and Chebotarev equidistribute prime ideals in this finite group. Choosing a global unit multiple of g permutes the global labels because the Eisenstein units have squareclasses ±1. Replacing g by its conjugate interchanges the two labels above p. Accordingly emptiness is invariant under global units and conjugation, so the ratio 13280/16384 computed with unquotiented tuples is also its proportion among split rational primes. Prime ideals of degree greater than one have density zero among prime ideals of F; each split rational prime contributes two degree-one prime ideals with the same empty/nonempty decision. Rational split and inert primes each have natural density one half. Combining the two certified subsets gives (1+13280/16384)/2 = 927/1024.

This argument uses no rank or parity assumption on the varying elliptic twists. The exact finite group computation and the per-prime decision procedure are supplied in prime_density.py. Counting primes directly is unnecessary for the density proof.

## Signed extension through 1000

For survivors without an archived Chabauty calculation, the implementation also tests the Eisenstein and rational forms modulo 49, 61 and 193. These are the same necessary square conditions with larger or additional moduli. They close the last squareclasses for both signs of 511, 527, 623 and 777. The exact signed replay finds 936 further complete lists, all with four points, among the 972 signed squarefree radicands between 201 and 1000 inclusive. Combined with the 243 earlier lists this gives 1,179 complete of 1,215. The remaining 36 are explicitly marked partial in extended-quadratic.json.
