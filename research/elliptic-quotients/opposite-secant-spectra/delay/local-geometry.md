# Sublattice criteria for local square cancellation

For an opposite-secant pair, write $h=\theta_{M_1}+\theta_{M_2}-\theta_{M_3}-\theta_{M_4}$ for its signed parent-neighbour theta sum. There is a direct test for $h$ to vanish at every index in the square class of $\mathbb Q_p$. It pairs finite collections of explicitly constructed sublattices by integral isometries. It explains all 18 five-adic cases in the [determinant census](README.md#extension-through-child-determinant-200000) and, with the layered construction below, all 90 two-adic cases. The isometries persist in three infinite families of geometrically realized pairs with equal elliptic quotient counts in all degrees and distinct ordinary theta series.

## The local criterion

Let $P_p$ retain coefficients at square units of $\mathbb Z_p$, and let $U(p^2)$ send $a(n)$ to $a(p^2n)$. The identities $P_p h=0$ and $U(p^2)h=0$ suffice. A local square index either has valuation zero and is counted by $P_p$, or is divisible by $p^2$ and is counted by $U(p^2)$. The latter identity also makes every subsequent $U(p^2)$ iterate vanish. The [parent transfer](README.md#reduction-to-the-parent-genus) then gives equality of geometric optimal elliptic quotient counts in every degree.

### Odd primes with split rank-two reduction

Suppose $M$ is one of the neighbours and $p$ is odd. Assume the reduction of its Gram form has a one-dimensional radical and is split on the two-dimensional quotient, and that $p^2$ divides its determinant. Over $\mathbb Z/p^2\mathbb Z$ choose a basis in which $Q(x,y,z)\equiv cxy$ with $c$ a unit. To obtain this basis, lift the two nonsingular isotropic directions to vectors of norm zero modulo $p^2$ by a linear correction using the gradient. Their pairing is a unit. The orthogonal complement of the lifted pair is generated locally by the cross product of the two Gram covectors. Its norm is divisible by $p^2$ by the determinant hypothesis. The implementation checks all these conditions directly.

Let $R$ be the preimage of the radical modulo $p$. For each of the $(p-1)/2$ lines in the binary quotient on which nonzero vectors have square norm, let $A_i$ be its preimage in $M$. Their intersection is $R$. The unit-square vectors are the disjoint union of $A_i\setminus R$, so
$$
 P_p\theta_M=\sum_i\theta_{A_i}-\frac{p-1}{2}\theta_R.
$$
For $0\le a,b\le2$, let $K_{ab}$ impose divisibility by $p^a$ on $x$ and by $p^b$ on $y$ in the split basis, with no condition on $z$. The set of vectors of norm divisible by $p^2$ is $K_{20}\cup K_{11}\cup K_{02}$. Inclusion–exclusion gives
$$
 U(p^2)\theta_M
 =\theta_{K_{20}/p}+\theta_{K_{11}/p}+\theta_{K_{02}/p}
  -\theta_{K_{21}/p}-\theta_{K_{12}/p}.
$$
The intersection $K_{22}$ cancels when adding the triple intersection. Division by $p$ in each lattice scales the norm by $1/p^2$. Thus the signed sums of radical lattices, of square-line lattices, and of the five-term divisibility expression can be checked separately. Matching each of the three signed lists by isometries proves the local exclusion. For $p=5$ every case in the prior census satisfies all three matchings. Each certificate uses 16 explicit isometry matrices across the three lists.

### A two-adic criterion

Assume $Q$ has integral Gram matrix and only norms $0,1\pmod4$. Then $Q(x)\equiv\ell(x)^2\pmod4$ for the linear parity functional $\ell(x)=Q(x)\bmod2$. Evaluating $Q$ on the coordinate vectors and their pairwise sums shows $B(x,y)\equiv\ell(x)\ell(y)\pmod2$. The norm-divisible-by-four vectors thus form the sublattice $R_2=\ker\ell$, and $U(4)\theta_M=\theta_{R_2/2}$. If $\ell(x)=1$, then
$$
 Q(x+2y)-Q(x)=4(B(x,y)+Q(y))\equiv0\pmod8.
$$
The unit-square condition therefore depends only on $x\bmod2M$. If the set of such residues is an affine coset $v+S\subset M/2M$, let $H$ and $J$ lift $S$ and $S+\mathbb F_2v$ respectively. The unit-square vectors are exactly $J\setminus H$, so $P_2\theta_M=\theta_J-\theta_H$. An empty set contributes zero. The same two isometry checks on the signed neighbour lists certify the exclusion. The condition on the residue sets and all the isometries are checked for the 25 two-adic cases with $U(4)h=0$.

For a second layer, apply the same idea to the signed sum $U(4)h$ on the four rescaled lattices $R_2/2$. Require an integral Gram matrix at this step. The norm-divisible-by-four residues still depend only on parity. If they form a subspace, take its preimage. Otherwise use the preimage of zero and the disjoint nonzero parity cosets. The coset of a nonzero residue $v$ is the difference between the lattices $2M+\mathbb Zv$ and $2M$. After scaling these terms by $1/2$, their theta sums give $U(4)^2h$.

For the unit-square test on $U(4)h$, fix a parity residue $v$ of norm one modulo four. Norm modulo eight now depends on residues modulo four. Among lifts $v+2y$, the congruence $Q(v+2y)\equiv1\pmod8$ is affine linear in $y\bmod2$, since its variable part divided by four is $B(v,y)+Q(y)$. The nonempty set of lifts is a coset of an even sublattice. Negation preserves the set, so its representative has order two modulo that sublattice. Its theta series is again a difference of two lattice theta series. Matching both second-layer lists proves $P_2U(4)h=0$ and $U(4)^2h=0$. Along with $P_2h=0$, this proves all local-square coefficients vanish.

The second-layer isometries explain 12 further cases. The remaining 53 two-adic cases are settled by the layered construction in the next section. With them, every one of the 108 local-obstruction cases in the census has a geometric proof.

## All valuation layers at two

The two layers above are the beginning of a recursion that needs neither the parity hypothesis nor a bound on the number of layers. [Layer certificates](two-adic-layers.json) · [Replay](dyadic_layers.py)

**Admissible forms.** Call a positive ternary Gram matrix admissible if its diagonal is integral and twice the matrix is integral. Let \(S_0\) be the signed list of the four neighbours. For an admissible form \(G\) on \(L\) put \(a=1\) when \(G\) is integral and \(a=2\) otherwise. Then the identity
\[
 Q(r+2^aw)-Q(r)=2^{a+1}B(r,w)+4^aQ(w)
\]
shows that divisibility of the norm by four depends only on \(r\bmod2^aL\). The same identity with \(a+1\) shows that the condition \(Q\equiv1\pmod8\) depends only on \(r\bmod2^{a+1}L\).

**The divisible vectors.** The vectors of norm divisible by four form a union \(D\) of cosets that is closed under all integer multiples. It is therefore a union of sublattices of \(L\) containing \(2^aL\). Let \(\mathcal F\) be the family of intersections of the maximal such sublattices. Möbius inversion on \(\mathcal F\) gives
\[
 1_D=\sum_{F\in\mathcal F}c_F1_F,\qquad c_F=1-\sum_{F'\supsetneq F}c_{F'}.
\]
Indeed, the members of \(\mathcal F\) containing a vector \(v\in D\) are exactly those containing the least member \(F_v\), and the coefficients sum to one over them. Every \(F\) has all norms divisible by four, so \(F\) with the form \(Q/4\) is admissible. Hence
\[
 U(4)\theta_L=\sum_Fc_F\,\theta_{(F,Q/4)},
\]
and the replay checks the indicator identity on every residue.

**The unit-square vectors.** The vectors \(P\) of norm \(1\bmod8\) are not a union of lattices, but a multiple of such a vector has norm \(1\bmod8\) or divisible by four. Take a residue-defined set \(Z\) containing all even multiples of \(P\) such that \(P\cup Z\) and \(Z\) are both unions of sublattices. Then
\[
 P_2\theta_L=\sum_F c_F\theta_F-\sum_{F'}c'_{F'}\theta_{F'}
\]
by the same inversion. The default is \(Z=D\). A unit-square set is unchanged on passing to the sublattice it generates, so one may also replace \(L\) by that sublattice first.

**Cancellation and closure.** A signed list is reduced by sorting its terms into isometry classes, by exact integral isometries checked with \(T^{\mathsf t}GT=G'\) and \(\det T=\pm1\), and summing their signs. Every construction depends only on the isometry class of \(G\). So \(R_j\), the reduced class vector of \(S_j=D(S_{j-1})\), represents \(U(4)^jh\) exactly. If the unit-square class vectors of \(R_0,\dots,R_{m-1}\) all vanish and \(R_m\) is a rational combination of \(R_0,\dots,R_{m-1}\), then linearity gives
\[
 P_2U(4)^jh=0\quad\text{for every }j\ge0 .
\]
This is exactly the vanishing of \(h\) at every square of \(\mathbb Q_2\), so the parent transfer gives equal elliptic quotient counts in every degree. No coefficient bound enters.

**The difficult cases.** Five cases need a unit-square identity involving all four neighbours. The indices are 71, 84, 100, 102 and 111. For each neighbour pair, every residue class \(1\bmod8\) modulo \(64\) has an explicit differing coefficient by norm \(512\). Thus no two of the four unit-square theta series can be identified within any of those progressions.

The search uses the indicator of \(P\) on \(L/4L\). It writes that indicator as a rational combination of the sets \(P\cap K\), over all sublattices \(4L\subseteq K\subseteq L\), and matches terms by the isometry class of the lattice generated by \(P\cap K\). The replay checks the indicator identities on every residue and the cancellation in every class. The five actual solutions have coefficients zero or one, so each yields a disjoint partition into positive pieces. Each piece contains eight unit-square residues modulo four. In the first four cases, each positive neighbour contributes one piece to each of the two negative neighbours. In case 111, each connection has two pieces. The explicit matrices in the [partition certificates](dyadic-partitions.json) verify these matches without the rational search.

Indices 71 and 100 have isometric collections of unit-square hulls, which are the lattices generated by the unit-square vectors. The same applies to 84 and 102 after exchanging the two signs. Thus these five inputs give three distinct collections of hulls at the initial valuation layer. Their Gram determinants are \(23296\), \(66304\) and \(13312\), respectively.

| Additional two-adic cases | Layers including the closing relation | Unit-square proof used |
|---|---:|---|
| 48 | 3 | complement \(D\) or the unit-square hull, with the general relation at layer 0 for 100 and 102 |
| 3 | 4 | complement \(D\) or the unit-square hull, with the general relation at layer 0 for 71, 84 and 111 |
| 2 | 6 | complement \(D\), unit-square hull, and residue complements modulo 16 |

Together the 53 cases use 57 hull steps, 47 complement-\(D\) steps, six steps with the complement given by norms \(0,4\bmod16\), and five general sublattice relations. The same procedure also re-proves the 37 cases of the earlier sections. It covers all 90 two-adic cases in the finite census. These closure lengths do not give a general termination bound.

For the thirteen nonvanishing controls, the checker verifies the nonzero coefficients at the archived two-adic square indices. Those witnesses reject square exclusion directly.


## Three infinite families

For each integer $t\ge0$, use the parent forms below. The two congruences on each line define the two index-three children in that family.

| Prime and valuation layers | Parent Gram matrix | Child congruences modulo 3 | Child Gram determinant |
|---|---|---|---:|
| 2, one layer | $\begin{pmatrix}5&1&1\\1&13&-3\\1&-3&17+72t\end{pmatrix}$ | $y-z=0$, $x+2y-z=0$ | $4608(2+9t)$ |
| 2, two layers | $\begin{pmatrix}8&0&2\\0&32&4\\2&4&41+288t\end{pmatrix}$ | $y-z=0$, $x+z=0$ | $18432(5+36t)$ |
| 5, one layer | $\begin{pmatrix}8&2&2\\2&13&3\\2&3&41+1800t\end{pmatrix}$ | $y-z=0$, $2x+2y-z=0$ | $36000(1+45t)$ |

The secant endpoints are fixed and the lift bases of their 3-neighbours are fixed because the parents are constant modulo nine. With the bases written in parent coordinates, the congruences defining their local sublattices are also fixed. At two the perturbations are divisible by eight and, in the two-layer family, by 32. After the first division of norms by four the second congruence test is therefore unchanged. At five the perturbation is divisible by 25. These congruences hold in 3-neighbour bases because their denominators are powers of three. The stored matching matrices preserve both the base norm and the square of the parent $z$-coordinate, an additional matrix identity checked separately for every isometry. They therefore preserve all the displayed deformations. The families have no distinguishing quotient degree.

In the one-layer two-adic family the children have two versus zero ordinary representations of norm five, and the signed neighbour coefficient at norm five is $-4$. In the two-layer family the children have two versus zero ordinary representations of norm eight and the signed neighbour coefficient at norm eight is $-4$. Thus $U(4)h$ has a nonzero coefficient at norm two in the latter family. In the five-adic family the children have two versus zero ordinary representations of norm eight. The signed neighbour coefficient is $-2$ at norm five when $t=0$, and $-4$ at norm eight when $t\ge1$. The first two parameters are checked by exact representation counts. For every $t\ge1$ the displayed norms see only $z=0$ because the perturbation alone contributes at least $72/9$, $288/9$ or $1800/9$ at any nonzero $z$ in a neighbour. This also applies to the children. The child lattices are inequivalent. The signed neighbour sums are nonzero, which rules out a pairwise isometry of the original neighbour collections as the explanation for these families.

All three parents have norms $0,1\pmod4$ throughout. The children do not represent one at $t=0$ and the deformations only increase norms. At two the parent vector $(0,0,1)$ has norm one modulo eight in all three families. For the two-adic families the binary restriction $z=0$ has determinant 64 and 256, respectively, so it is nondegenerate at every odd prime. In the five-adic family the restriction has determinant 100, and at five $(0,0,1)$ has norm one modulo five. A nondegenerate binary form over an odd residue field represents one, and the usual smooth lifting applies. At three each child contains a unit-square line because it is a secant. The [primitive CM realization criterion](README.md#exact-enumeration-and-local-tests) now supplies genus-two Jacobians for every parameter in all three families.

## Replay

~~~sh
study=research/elliptic-quotients/opposite-secant-spectra/delay
python -B "$study/local_geometry.py" --verify "$study/local-geometry.json"
python -B "$study/dyadic_layers.py" --verify "$study/two-adic-layers.json"
python -B "$study/dyadic_layers.py" --all-dyadic --controls
~~~

The replay constructs each 3-neighbour and each sublattice directly from the registered parent and secant endpoints. It verifies every isometry by an integral matrix identity and a determinant check. For all three infinite families it verifies the simultaneous matrix identities and the bounded ordinary representation checks. The second command compares the 53-case receipt. The third replays all 90 two-adic cases and runs the thirteen nonvanishing controls and the compact partition checks. The cancellation proofs do not use modular-form coefficient cutoffs. [Exact matrix certificates](local-geometry.json)
