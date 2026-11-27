# Sublattice criteria for local square cancellation

For an opposite-secant pair, write $h=\theta_{M_1}+\theta_{M_2}-\theta_{M_3}-\theta_{M_4}$ for its signed parent-neighbour theta sum. There is a direct test for $h$ to vanish at every index in the square class of $\mathbb Q_p$. It pairs finite collections of explicitly constructed sublattices by integral isometries. It explains all 18 five-adic cases in the [determinant census](README.md#extension-through-child-determinant-200000) and 37 two-adic cases, including all 25 in which $U(4)h=0$. The isometries persist in three infinite families of geometrically realized pairs with equal elliptic quotient counts in all degrees and distinct ordinary theta series.

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

The second-layer isometries explain 12 further cases. The remaining 53 two-adic cases already have modular-form certificates. Extending the isometry construction to them remains open.

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
~~~

The replay constructs each 3-neighbour and each sublattice directly from the registered parent and secant endpoints. It verifies every isometry by an integral matrix identity and a determinant check. For all three infinite families it verifies the simultaneous matrix identities and the bounded ordinary representation checks. The conclusions do not use modular-form coefficient cutoffs. [Exact matrix certificates](local-geometry.json)
