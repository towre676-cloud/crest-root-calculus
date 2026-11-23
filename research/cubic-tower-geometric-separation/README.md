# Geometric simple factors through the sixth cubic-tower stage

The genus-2,917 Jacobian at stage six has **127 geometrically simple, pairwise geometrically nonisogenous factors**, all defined over the rationals. The previous rational trace-rank result allowed any of the larger signed factors to split after extending the ground field. Here good reductions impose incompatible dimensions on every proposed proper geometric subvariety. Powered Frobenius polynomials then separate the factors of each common dimension. [Tower and projectors](../cubic-tower-factors/README.md) · [Rational trace rank](../cubic-tower-higher-supports/README.md)

| Support in the tenth stage | Number of certified factors | Dimension of each | Geometric separation |
|---|---:|---:|---|
| Base, singletons, pairs | 111 | 1 or 4 | Earlier absolute-simplicity and separation certificates |
| All triples | 240 | 12 | All 28,680 pairs separated at 41 |
| All four-index subsets | 420 | 32 | 87,989 of 87,990 pairs separated at 41, the last at 53 |
| Five-index subsets of $\{0,1,2,3,4,5\}$ | 12 | 80 | All 66 pairs separated at 83 |
| $\{0,1,2,3,4,5\}$ | 2 | 192 | The signs separated at 193 |

These certificates identify **785 geometric simple classes** in the tenth stage, of total dimension 18,045. The signed products on the first three, four, five and six indices are complete geometric decompositions with 15, 31, 63 and 127 simple factors, respectively. Their genera are 55, 217, 811 and 2,917. A further calculation settles all 504 five-index factors and separates all 420 six-index factors. [Larger support certificates](../cubic-tower-size-five/README.md) · [Complete seventh stage](../cubic-tower-stage-seven/README.md). The integral higher-support Frobenius polynomials also remain open. The [Cartier calculation](../cubic-tower-cartier/README.md) gives an independent construction and the 240 triple simplicity proofs. The [certificate](certificate.json) here supplies the separation data and the 434 additional absolute-simplicity proofs in dimensions 32, 80 and 192.

## Frobenius on the finite fibres

Let $S$ be a nonempty support of size $k$. Write $f=-4z^3+z^2+z+1$ and keep the [simultaneous good-reduction conditions](../cubic-tower-factors/README.md) for the indices actually used. Let $d=k2^{k-1}$ be the dimension of either signed abelian factor. We give an independent algorithm for its entire numerator modulo a prime $p>\max(3,d)$. Its finite matrices use the polynomial ring in $z$ over the prime field. No enumeration of $\mathbb F_{p^d}$ is needed.

Over the finite $z$-line the coordinate ring is free over $\mathbb F_p[z]$, with equations and monomial basis

$$
 w^2=f(z),\qquad t_i^3+3(z+i)t_i=w,\qquad w^b\prod_{i\in S}t_i^{e_i},\quad
 b\in\{0,1\},\quad e_i\in\{0,1,2\}.
$$

The normalized cubic covers over the elliptic base have disjoint finite branch loci, so the affine fibre product gives the required points even when one coordinate ramifies. Quotient this free module by the sum of the smaller modules omitting one or more indices. Its basis consists of the terms with all $e_i\in\{1,2\}$. This quotient realizes the alternating subset trace, by the coordinate decomposition of the free basis. The involution negates $w$ and all $t_i$. For sign $(-1)^\delta$, the parity condition $b+\sum e_i\equiv\delta\pmod2$ leaves $2^k$ basis monomials $v_j$.

The $p$-power map is semilinear in $z$. Define $M_{ij}\in\mathbb F_p[z]$ by $v_j^p=\sum_i M_{ij}(z)v_i$ on this quotient. The [implementation](affine_frobenius.py) reduces each cubic power separately. Multiplication by $t_i$ takes $(c_0,c_1,c_2)$ to $(w c_2,c_0-3(z+i)c_2,c_1)$ in the basis $(1,t_i,t_i^2)$. Suppressing the uniquely determined power of $w$ modulo two leaves ordinary polynomials in $z$. The relation $w^2=f$ supplies every missing pair of powers of $w$.

Give $z,w,t_i$ weights $2,3,1$. Both defining reductions preserve the weight bound. If $W_i$ is the weight of basis monomial $v_i$, then $2\deg M_{ij}\le p W_j-W_i$. Construct the following prime-field matrix $A$ using one positive integer $a\le\lfloor W_i/2\rfloor$ for each row label $i$, and analogously for its columns:

$$
 A_{(i,a),(j,b)}=[z^{pb-a}]M_{ij}(z).
$$

For $q=p^n$, the finite signed fibre trace summed over $z\in\mathbb F_q$ equals $-\operatorname{Tr}(A^n)$ modulo $p$. Indeed, at a given $z$ the $q$-power action has matrix $M(z)M(z^p)\cdots M(z^{p^{n-1}})$. In a finite fibre algebra its trace counts rational geometric points modulo $p$. Frobenius is nilpotent on the nilradical, and its trace on the reduced functions is the trace of the permutation of the geometric points. The same assertion applies to Frobenius composed with the involution. Passing to the signed subset quotient therefore gives the required projected finite trace.

Summing a power $z^h$ over $\mathbb F_q$ gives $-1$ for positive multiples of $q-1$ and zero otherwise. In a product term along a cycle, write the exponents as $h_0,\ldots,h_{n-1}$. Whenever $\sum_r p^r h_r=(q-1)a_0>0$, the cyclic equations $h_r=p a_{r+1}-a_r$ have a unique solution in positive integers. Divisibility follows successively from the displayed equality modulo $p$, and positivity follows from $h_r\ge0$. The weight bound implies $p(a_{r+1}-W_{i_{r+1}}/2)\le a_r-W_{i_r}/2$. Cycling gives $a_r\le W_{i_r}/2$. Thus these terms are exactly the ones counted by $\operatorname{Tr}(A^n)$, including their multiplicities.

The normalized fibre at infinity contributes an additional trace. Put $c_i$ equal to the quadratic character of $-(12i+1)/12$ in $\mathbb F_p$. In sign $(-1)^\delta$ its $n$th trace is

$$
 b_n=\frac12\left(\prod_{i\in S}(1+c_i^n)+(-1)^\delta\prod_{i\in S}(1-c_i^n)\right).
$$

This follows from the single rational cubic sheet and the remaining pair exchanged by the involution. The nonempty subset projection cancels $H^0$ and $H^2$, leaving $\operatorname{Tr}(F^n\mid H^1(A_S^{(-1)^\delta}))\equiv\operatorname{Tr}(A^n)-b_n\pmod p$. Here $b_{2n}=r=2^{k-1}$, so the infinity eigenvalues $1,-1$ have multiplicities $(r+b_1)/2,(r-b_1)/2$. Compute $\det(1-TA)$ using exact modular similarity to an upper Hessenberg matrix and its principal-minor recurrence. Divide the resulting power series by $(1-T)^{(r+b_1)/2}(1+T)^{(r-b_1)/2}$ through degree $d$. Newton's identities determine these coefficients because $p>d$. Reciprocity forces every coefficient above $d$ to vanish modulo $p$. This gives the entire reduced Weil polynomial. Its degree equals the $p$-rank and is $d$ precisely for ordinary reduction.

## Certifying geometric isogeny information

For factors of dimension $d$, put

$$
 M_d=\prod_{\ell\ \mathrm{prime},\,\ell-1\le4d^2}
 \max\{\ell^a:\ell^{a-1}(\ell-1)\le4d^2\}.
$$

Any ratio of two Frobenius roots from such factors has degree at most $4d^2$. If the ratio is a root of unity, its order divides $M_d$, since each prime power in that order has Euler totient at most $4d^2$. Raising the roots to $M_d$ therefore identifies every pair that could coincide over a finite extension. This bound applies between two different factors as well as within one polynomial.

If two characteristic-zero factors were geometrically isogenous, their good reductions would be isogenous over a finite extension. Some power of the two Frobenius multisets would agree. Pairing the corresponding roots, the preceding bound shows that already their $M_d$th powers agree over characteristic zero. Reducing those multisets at $p$ preserves equality of their nonzero elements. The [program](geometric.py) factors the observed prime-field unit polynomial. For each irreducible of degree $m$ it raises a companion matrix to $M_d\bmod(p^m-1)$, then gathers the resulting irreducible factors with multiplicity. A difference of these multisets excludes geometric isogeny. At 41 the sole four-index collision is $(1,3,4,8)^+$ against $(2,3,4,5)^-$. Their powered polynomials differ at 53.

For absolute simplicity, suppose a proper geometric abelian subvariety of dimension $e\le d/2$ existed. Such an $e$ exists whenever a geometric decomposition has more than one constituent, by complete reducibility. After extending the ground field it specializes at every witness prime. Use a prime where the ambient reduction is ordinary. The subvariety is then ordinary, so precisely $e$ of its $2e$ Frobenius roots are units. For a finite extension on which it is defined, its characteristic polynomial is a rational submultiset of the full Frobenius roots raised to some power $nM_d$. Further powering is injective on the **distinct** characteristic-zero $M_d$th powers, by construction. It is Galois equivariant. Consequently this rational submultiset, with its multiplicities, lifts to a rational submultiset at $M_d$. Its roots are algebraic integers, and taking powers preserves which roots are units. Reducing at $p$ gives a monic divisor of degree $e$ of the powered unit polynomial. The degrees of the irreducibles, repeated with their multiplicities, therefore restrict $e$ to their subset sums. Intersect these restrictions over the listed good primes. The intersection is empty for every factor in the certificate. This proves geometric simplicity globally even when different local primes exclude different dimensions. The [Cartier verifier](../cubic-tower-cartier/geometric_dimensions.py) implements the dimension calculation using orders in finite fields. [Howe and Zhu, absolute simplicity over finite fields](https://arxiv.org/abs/math/0002205)

## Replay

With Python 3.11 or later, NumPy and SymPy installed, run these commands from this directory. On Windows, `py` can replace `python`.

```sh
python verify.py
python verify.py --recount --workers 4
python verify.py --cartier-cross
python affine_frobenius.py 41 0 1 2 3
```

The first command checks every stored good-prime and dimension witness, separates every pair in the four groups, and reevaluates a higher-support matrix in each dimension. Both commands also recompute all 110 earlier singleton and pair polynomials modulo 41 and compare two extension counts in each of dimensions 12, 32, 80 and 192 against direct fibre enumeration. A separate Cartier matrix construction checks all 32 coefficients of both signs on the first four-index support. The second command rebuilds all stored higher-support observations from the curve equations. The optional Cartier comparison independently recomputes all 420 four-index polynomials at 41 using holomorphic differentials. Their 13,860 coefficients agreed byte for byte, with SHA-256 `6090bef37b7b25bc338c1c9bf54f0ccc160937f300bc156208a9b2146d9601bf` when concatenated in subset/sign order. The standalone coefficient command prints both reduced Frobenius polynomials for one support. The certificate records coefficients in increasing powers of $T$. The callable `signed_numerators(p, support)` returns positive and negative polynomials in this order, through the middle coefficient. The coefficient algorithm requires simultaneous good reduction on that support and $p$ larger than its signed dimension.
