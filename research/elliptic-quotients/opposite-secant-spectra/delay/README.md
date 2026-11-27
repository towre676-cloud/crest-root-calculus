# First differences and equal elliptic quotient spectra

The level-three opposite-secant construction can produce **equal quotient counts in every degree**, even when the refined Humbert forms are inequivalent. A family below proves this for infinitely many pairs of genus-two curves. An exhaustive sweep through child Gram determinant 200,000 finds 280 such pairs and measures the first difference for every pair in the sweep whose spectra differ.

Three further infinite families have identical full spectra through explicit sublattice cancellations at 2 and 5. The same criterion explains 55 cases in the determinant census by exact isometries. [Local constructions and matrix certificates](local-geometry.md)

A separate effective family has finite first differences beyond any fixed multiple of the third successive-minimum length, even while the three parent minima become asymptotically equal. Along that family $\delta\gg\lambda_3(G)^{4/7}$ for the **squared** third minimum. [Construction and exact replay](balanced-delay.md)

Here $r_F(n)$ is the number of vectors of norm $n$ in a ternary lattice $F$, with both signs counted. Write $s_F(m)=r_F(m^2)$. Divisor inversion passes between these ordinary counts and primitive square counts. In particular the first differing degree and its signed excess are the same for either convention. On a refined Humbert lattice the primitive counts give geometric optimal elliptic quotient classes. [Kani, elliptic subcovers](https://doi.org/10.1007/s40316-018-0105-6)

## Reduction to the parent genus

Let $G$ be an integral positive ternary lattice with $\det G\equiv1\bmod3$. Its isotropic conic over $\mathbb F_3$ has four points. Partition them into opposite secants $e,e'$, and let $L_e,L_{e'}$ be the inverse images of their two planes in $G$. They have identical counts in every degree prime to three by [the opposite-secant argument](../README.md#the-level-three-classification).

For each isotropic point $v$, choose an integral lift $w$ with $Q(w)\equiv0\bmod9$. This is possible by varying the lift by $3y$ because the form is nondegenerate modulo three. Define
$$
 H_v=\{y\in G:B(y,w)\equiv0\bmod3\},\qquad
 M_v=H_v+\mathbb Z(w/3).
$$
Here $B$ is the integral Gram pairing. The lattices $M_v$ are the four 3-neighbours of $G$. They have the same determinant and agree locally with $G$ away from three. At three they are unimodular of the same rank and determinant square class, so they too lie in its genus. Write
$$
 A(m)=\sum_{v\in e}s_{M_v}(m)-\sum_{v\in e'}s_{M_v}(m),\qquad
 B_0(m)=\sum_{v\in e}s_{H_v}(m)-\sum_{v\in e'}s_{H_v}(m).
$$

For $a\not\equiv0\bmod3$, an ambient vector $aw+3y$ has norm divisible by nine exactly when $B(y,w)\equiv0\bmod3$. Thus the vectors of norm $9m^2$ with nonzero residue on $v$ are in bijection under division by three with the vectors of norm $m^2$ in $M_v\setminus H_v$. The zero residue makes a common contribution. On the tangent plane $v^\perp$, the quotient by its radical has norm square class $-\det G\equiv2\bmod3$. It follows that $H_v$ represents no square prime to three, and the only vectors reducing to an isotropic residue in $H_v$ have residue on $v$ or zero. Consequently
$$
 s_{L_e}(3m)-s_{L_{e'}}(3m)=A(m)-B_0(m)=B_0(3m),\qquad
 B_0(m)=0\quad (3\nmid m).
$$
In the middle identity the vectors of zero residue in each $H_v$ contribute $s_G(m)$ and cancel. Induction gives the first-difference formula
$$
 \boxed{\quad\delta(L_e,L_{e'})=3\min\{m\ge1:A(m)\ne0\}\quad}
$$
with both sides infinite if the indicated sets are empty. At the first difference the excess is $A(m)$ itself. In particular matching the two unordered pairs of neighbours by isometry proves equality of the entire quotient spectra.

The transfer also gives the complete difference sequence. Let $D_{mathrm{prim}}(m)$ be the signed difference of the two primitive square counts and let $A_{mathrm{prim}}(m)$ be the corresponding signed sum of the four parent-neighbour counts. Then
$$
 D_{mathrm{prim}}(m)=
 sum_{j=1}^{v_3(m)}(-1)^{j-1}A_{mathrm{prim}}(m/3^j),
 qquad
 A_{mathrm{prim}}(m)=D_{mathrm{prim}}(3m)+D_{mathrm{prim}}(m).
$$
Indeed $B_0(m)=0$ if $3\nmid m$ and $B_0(3m)$ is the ordinary child difference at $3m$. The ordinary relation is therefore $D(3m)=A(m)-D(m)$. With $T_3 f(m)=f(m/3)$ when $3\mid m$ and zero otherwise, this is $(1+T_3)D=T_3A$. Divisor inversion commutes with $T_3$, so the relation holds for primitive counts and can be expanded at each $m$ as the displayed finite sum.

When neighbour isometries do not settle the question, apply the $t=1$ Shimura lift to the signed sum of the four neighbour theta series. It is cuspidal in weight $3/2$ because all four neighbours lie in one genus. At theta level $N_G$, its lift is holomorphic of weight two, trivial character and level $N_G/2$. Holomorphic weight two suffices for Sturm's test, including when a unary component has a noncuspidal image. Its coefficient of index $n$ is $A(n)$ plus a linear combination of $A(d)$ for proper divisors $d\mid n$. Thus the least nonzero $A(m)$ is at most $[\mathrm{SL}_2(\mathbb Z):\Gamma_0(N_G/2)]/6$, unless all $A(m)$ vanish. This moves the Sturm test to the **parent level**. Compared with the child theta level, which has an extra factor of nine at three, the resulting first-degree bound is about three times smaller. [Shimura–Niwa level and holomorphy, in Purkait, Theorem 3.1](https://doi.org/10.5802/jtnb.865)

## An infinite family with identical full spectra

For $t\ge0$ integral put
$$
 G_t=\begin{pmatrix}5+36t&1-36t&2\\1-36t&5+36t&2\\2&2&8\end{pmatrix},
 \qquad Q_t=Q_0+36t(x-y)^2.
$$
Take $L_e$ with $x+2y+2z\equiv0\bmod3$ and $L_{e'}$ with $x+z\equiv0\bmod3$. Their endpoints, in order, are
$$
 ((0,1,2),(1,1,0)),\qquad ((1,0,2),(1,1,2)).
$$
For every $t\not\equiv3\bmod5$, the two lattices are refined Humbert forms of geometrically nonisomorphic genus-two curves. Their optimal elliptic quotient counts agree in every degree. Their ordinary theta series first differ at norm $5+36t$, and their child Gram determinant is $1440(1+18t)$.

**Neighbour isometries.** In endpoint order, the following matrices are three times integral bases of the four neighbours, where columns are the basis vectors of the scaled lattices.
$$
 W_1=\begin{pmatrix}9&0&0\\0&3&2\\0&0&1\end{pmatrix},\quad
 W_2=\begin{pmatrix}3&0&1\\0&3&1\\0&0&3\end{pmatrix},\quad
 W_3=\begin{pmatrix}3&0&2\\0&9&0\\0&0&1\end{pmatrix},\quad
 W_4=\begin{pmatrix}9&6&2\\0&3&2\\0&0&1\end{pmatrix}.
$$
Equivalently each is a basis of the lattice generated by $3H_v$ and the lift $w$ of the given endpoint. Since $G_t\equiv G_0\bmod9$, the bases work for the whole family. The first isometry is $(x,y,z)\mapsto(-y,-x,-z)$. The second has matrix
$$
 R=\frac19\begin{pmatrix}2&-7&-8\\-7&2&-8\\6&6&3\end{pmatrix}.
$$
Both preserve $Q_0$ and $(x-y)^2$, and carry $W_3\mathbb Z^3$ to $W_1\mathbb Z^3$ and $W_4\mathbb Z^3$ to $W_2\mathbb Z^3$, respectively. For the latter, $W_2^{-1}R W_4$ has rows $(0,-1,-1),(-3,-2,-1),(2,2,1)$ and determinant $\pm1$. The first-difference formula proves equality at all squares.

**Inequivalence.** On the plane $x=y$, the reflection $(x,x,z)\mapsto(x,x,-x-z)$ preserves the norm and exchanges the two restrictions. Put $s=x+y$ and $r=x-y$. Then
$$
 Q_0=2r^2+\tfrac52s^2+8(z+s/4)^2.
$$
Off the plane $r=0$, the smallest possible value is five, attained exactly by the signed first and second coordinate vectors. This follows because $|r|\ge2$ gives norm at least eight, and if $|r|=1$ then $s$ is odd. At norm $5+36t$ the second coordinate vector and its negative belong only to $L_{e'}$. The first coordinate vector belongs to neither lattice. There is no other difference below that norm because the plane contributions agree. Hence the forms are inequivalent for every $t$.

**Realization.** The form is $(x+y)^2\bmod4$ and neither child represents one. For an explicit unbounded subsequence, take a positive $r\equiv7\bmod(36\cdot5\cdot47)$ and $t=(r^2-49)/36$. The common vector $(2,1,1)$ has norm $r^2$ and the root is coprime to $1440(1+18t)=720(r^2-47)$. This proves the realization criterion directly on that subsequence. For all $t\not\equiv3\bmod5$, the local test below proves the same criterion. At two, $(2,1,1)$ has norm $1\bmod8$ when $t$ is even, and $(1,0,0)$ does when $t$ is odd. At three both secants contain a line of unit square norm. At every odd prime other than three and five, the restriction $x=y$ is nondegenerate and represents one modulo the prime. At five that restriction has rank one, but the coefficient of $(x-y)^2$ is a unit when $t\not\equiv3\bmod5$, so the full lattice has rank at least two. The realization theorem now gives the claimed Jacobians. The primitive binary form in the genus description has discriminant $-360(1+18t)$, so taking $1+18t$ to be primes congruent to one modulo 18 yields infinitely many CM fields. Taking $1+18t=s^2$ with $s\equiv1\bmod18$ prime to five keeps the field $\mathbb Q(\sqrt{-10})$ fixed. [Kani, CM realization and binary genus description](https://www.queensu.ca/academia/kani/sites/acadkaniwww/files/uploaded_files/Preprints/refHum5.pdf)

## Exhaustive sweep through child determinant 100,000

The sweep ranges over **integral Gram** parents of determinant prime to three with $\det G\equiv1\bmod3$ and norms $0,1\bmod4$. The latter condition is necessary for the children since their index in $G$ is odd. The primitive realization criterion and the norm-one test select the Jacobian pairs within that range. The sweep concerns level-three opposite secants of these parents.

| Quantity | Count |
|---|---:|
| Reduced Gram matrices before isometry removal | 50,332 |
| Parent lattice isometry classes | 40,180 |
| Unordered pairs of inequivalent children | 90,027 |
| Pairs with a norm-one vector in both children | 2,689 |
| Other pairs obstructed by a local unit-square test | 16,600 |
| Pairs realized by genus-two Jacobians under the primitive criterion | 70,738 |
| Realized pairs with a finite first difference | 70,608 |
| Realized pairs with exactly equal full spectra | 130 |

Of the 130 equalities, **74** follow from neighbour isometries and **56** from Shimura–Sturm on the parent. In each of those 56 cases the signed sum of the four ordinary neighbour theta series is nonzero, even though its square coefficients all vanish. The first nonzero ordinary norm is 5 in 35 cases, 8 in 15, 12 in two, 13 in three and 20 in one. The smallest equal-spectrum pair has determinant 1,440 and is $t=0$ in the family above. All 31 pairs in the earlier census are recovered by exact child isometry tests. The determinant cutoff concerns lattices. These counts do not give the number of curve isomorphism classes realizing each lattice.

| Extremum over pairs with differing spectra | First difference $\delta$ | Child determinant | Parent squared minima |
|---|---:|---:|---|
| Largest $\delta$ | 66 | 83,088 | $(8,29,49)$ |
| Largest $\delta/\sqrt{\det L}$, about 0.360 | 60 | 27,792 | $(8,13,32)$ |
| Largest $\delta/\sqrt{\lambda_3(G)}$, about 11.14 | 60 | 85,680 | $(20,24,29)$ |

Here $\lambda_3(G)$ denotes the **squared** third successive minimum. The table is evidence about finite first differences. It puts no upper bound on the asymptotic exponent. Write $\Delta(X)$ for the largest finite first difference over pairs of determinant at most $X$. The [delayed-difference family](../README.md#an-infinite-family-with-an-unbounded-first-difference) gives $\Delta(X)\gg\sqrt X$. The Shimura–Sturm argument gives $\Delta(X)\ll X\log\log X$ at large $X$, because the theta level of an integral Gram lattice is at most four times its determinant. A power bound between these remains open.

## Extension through child determinant 200,000

The extended replay recovers all 31 archived pairs and agrees with the original totals in each earlier determinant band. It enumerates 150,300 parent isometry classes and 364,337 pairs of inequivalent children. Of these, 289,274 pairs meet the primitive Jacobian criterion. Their spectra differ for 288,994 pairs and agree in every degree for 280. Neighbour isometries account for 159 of the equalities. The other 121 have a nonzero ordinary neighbour theta difference by norm 28 and zero square coefficients by Sturm.

Of those 121 nonzero sums, 90 are supported entirely on nonsquares over $\mathbb Q_2$ and 18 over $\mathbb Q_5$. The other 13 have explicit nonzero indices that are squares at every prime of the theta level. The Hecke decomposition accounts for all 13. Their central blocks belong to eight weight-two newform orbits with exact twisted central zeros. Four of the 13 sums also have blocks supported on nonsquares over $\mathbb Q_2$. [Local-square criterion](local-square.md) · [Hecke packets and exact replay](residual-packets.md)

| Extremum over pairs with differing spectra | First difference $\delta$ | Child determinant | Parent squared minima |
|---|---:|---:|---|
| Largest $\delta$ | 102 | 189,360 | $(9,9,264)$ |
| Largest $\delta/\sqrt{\det L}$, about 0.360 | 60 | 27,792 | $(8,13,32)$ |
| Largest $\delta/\sqrt{\lambda_3(G)}$, about 14.10 | 81 | 182,988 | $(28,28,33)$ |

The last two extrema measure different effects. The largest finite degree occurs on an elongated parent, while the greatest normalized delay relative to the third minimum occurs when all three squared minima are close together. The [nearly round construction](balanced-delay.md) proves that the latter normalized delay is unbounded. The optimal exponent in the child determinant remains open. [Exact extended receipt](census-200000.json)

### Exact enumeration and local tests

For a Minkowski-reduced ternary Gram matrix with diagonal $a\le b\le c$ and off-diagonal entries $u,v,w$, the inequalities are $2|u|,2|v|\le a$, $2|w|\le b$, and $2(\pm u\pm v\pm w)\ge-a-b$ for the four sign patterns of positive product. Its diagonal entries are the squared minima, and $abc\le2\det G$. The generator enumerates these bounded integer entries, imposing the mod-four condition via diagonal residues and off-diagonal parities. Exact isometries remove duplicates at fixed determinant and minima. Each split child has a cyclic order-nine part of its discriminant group at three, so its unique integral index-three overlattice is its parent. Exact parent automorphisms then remove duplicate unordered pairs. [Minkowski inequalities and determinant bound](https://doi.org/10.1017/S1446788700024964)

For the realization test it suffices to represent one over $\mathbb Z_p$ for $p$ dividing twice the child determinant. At two the test is representation of $1\bmod8$, and at an odd prime it is representation of $1\bmod p$. An odd-prime solution has nonzero gradient and lifts. At all other odd primes the lattice is unimodular of rank three and the test is automatic. Hasse–Minkowski and weak approximation applied to $Q(\mathbf x)=m^2$ then give a rational point with affine coordinates $\mathbf x/m$ integral at the finitely many tested primes. Clearing denominators gives an actual square with root prime to the determinant. This also implies primitivity of the form. Checking for norm one then distinguishes the product cases. [Weak approximation for smooth quadrics](https://www.aimath.org/WWN/qptsurface2/articles/html/9a/) · [Primitive CM realization](https://www.queensu.ca/academia/kani/sites/acadkaniwww/files/uploaded_files/Preprints/refHum5.pdf)

### Reproduction

From the repository root, with a C++17 compiler and the [opposite-secant replay dependencies](../README.md#reproduction-and-prior-work) installed, run

~~~sh
study=research/elliptic-quotients/opposite-secant-spectra/delay
python -B "$study/census.py" --verify "$study/census.json"
python -B "$study/census.py" --determinant-limit 200000 --verify "$study/census-200000.json"
~~~

The census uses integer square-shell enumeration and exact PARI isometries. PARI independently checks each finite extremum and the primitive transfer identity through degree 81, including the corresponding counts on the parent neighbours. The family replay tests five parameters, including the ordinary theta difference and the exact neighbour isometries. The two proofs above apply to the entire stated families.
