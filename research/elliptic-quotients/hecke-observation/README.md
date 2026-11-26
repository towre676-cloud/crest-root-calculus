# Prime-power separation and the oldforms retained by square counts

At the 36 inert self-gluing levels in the [replay](prime-powers.json), the numbers of optimal elliptic quotients in degrees $1,p,p^2,\ldots$ distinguish the 7,580 refined Humbert forms counted across the levels. The auxiliary prime $p$ is fixed at each level. The levels cover all nine maximal imaginary quadratic orders of class number one and go through gluing prime $103$. At most ten steps in the $p$-neighbour graph suffice. The comparison includes pairs in different genera.

The choice of $p$ matters. At $D=-3$, $\ell=5$, the prime $p=17$ leaves two pairs tied at every power of $17$. Both pairs separate in degree two, and $p=7$ separates every pair at the same level.

## From neighbours to prime-power counts

Fix $D,\ell$ and an auxiliary prime $p\nmid 2D\ell$. Join every integral ternary lattice in the $p$-neighbour closure to its $p+1$ neighbours with multiplicity. The closure is finite because each neighbour stays in the genus. Let $A$ be the matrix with row $i$ recording the neighbours of class $i$. The column $R_j$ gives the ordinary counts of norm $p^{2j}$ on all classes, including both vector signs. Then
$$
R_1=(A-\chi I)R_0,\qquad
R_{j+1}=A R_j-pR_{j-1}\quad(j\ge1),\qquad
\chi=\left(\frac{D}{p}\right).
$$
Here $R_0$ records norm-one counts. A gluing lattice has zero norm-one count, while its neighbours can have norm-one vectors. The geometric quotient multiplicity at degree $p^j$, for $j\ge1$, is the corresponding entry of $R_j-R_{j-1}$. Geometric multiplicities count both complementary elliptic directions, as in the [inert gluing theorem](../inert-gluing/README.md).

For the recurrence, a $p$-neighbour belongs to an isotropic line in $L/pL$. A vector of norm $n$ outside $L$ in one of these neighbours is $y/p$ for a unique $y\in L\setminus pL$ of norm $p^2n$. The line of $y$ determines the neighbour. These vectors contribute $r_L(p^2n)-r_L(n)$ to the sum of the counts over neighbours. A vector $x\in L\setminus pL$ belongs to the neighbours whose lines are perpendicular to its reduction. There are $1+\left(\frac{-\det(L)n}{p}\right)$ such lines, including the isotropic case, where the number is one. A vector in $pL$ belongs to every neighbour. Summing gives
$$
\sum_{L'\text{ a }p\text{-neighbour}}r_{L'}(n)
=r_L(p^2n)+\left(\frac{-\det(L)n}{p}\right)r_L(n)+p\,r_L(n/p^2),
$$
where a nonintegral index contributes zero. Here $\det(L)=4|D|\ell^2$, so the formula at $n=p^{2j}$ gives the claimed recurrence. This is the square-index specialization of the ternary Eichler commutation relation. [Classical neighbour and theta calculation](https://doi.org/10.1017/S0027763000019218)

Each $R_j$ is a monic degree-$j$ polynomial in $A$ applied to $R_0$. Thus two profiles agree at every power of $p$ exactly when their coordinates agree on the columns $R_0,AR_0,\ldots,A^{h-1}R_0$, where $h$ is the size of the relevant closed graph. Cayley–Hamilton gives the reverse implication. The main replay distinguishes profiles modulo $65\,521$, which suffices for exact inequality over the integers. It also checks the first three degrees $1,p,p^2$ against direct primitive enumeration. The all-powers collisions below use an exact invariant subspace over the integers.

### An exact obstruction for $p=17$

At $D=-3,\ell=5$, write $\omega=(1+\sqrt{-3})/2$. The positive genus has scalar parameters $0,1+\omega$, in that order. Its $17$-neighbour matrix with the gluing forms first, and its norm-one column, are
$$
A_+=\begin{pmatrix}0&12&6\\2&10&6\\2&12&4\end{pmatrix},
\qquad R_0=(0,0,2)^{\mathsf T}.
$$
The first two rows have the same sum on the first two columns and the same third entry, so $A_+$ preserves equality of the first two coordinates. Thus every degree $17^j$ gives the same geometric multiplicity on the pair, although degree two gives $6$ versus $2$. The negative genus has a second pair, $\omega,1+2\omega$, with the same invariant-subspace argument and degree-two counts $0$ versus $2$.

## A square-null cuspidal eigensystem that occurs in just one genus

At $D=-7,\ell=5$, exhaustive integral enumeration gives genus sizes $5$ and $6$. The $T(3^2)$ eigenvalue $1$ occurs once, in the negative genus. Its theta image is a nonzero cusp form $h$, with $[q^5]h=8$ in the recorded basis. The first $240$ squares give zero, and the weight-two Sturm bound at $M=700$ certifies $\mathcal S_1 h=0$ for the Shimura lift $\mathcal S_t$ in square class $t$.

The same bound applied to $\mathcal S_5 h$ proves
$$
\mathcal S_5 h=8\bigl(F(z)+F(2z)\bigr),
$$
where $F$ is the normalized weight-two newform of the curve $[0,1,1,-131,-650]$ of conductor $35$. The exact modular symbol gives $L(F,1)/\Omega^+=1$. The twist by $-7$ has conductor $245$ and root number $-1$, so $L(F\otimes\chi_{-7},1)=0$. This identifies an actual twisted central-value vanishing for one square-null eigenline. The all-square vanishing was also verified independently by the coefficient calculation. [Curve 35.a1](https://www.lmfdb.org/EllipticCurve/Q/35/a/1)

## A shared eigensystem with three independent square lifts

At $D=-11,\ell=7$, the full genera have $15$ and $17$ classes. The $T(3^2)$ eigenspace of eigenvalue $-1$ has one dimension in the positive genus and two in the negative genus. Let $F$ be the normalized weight-two newform of the curve $[0,-1,1,-10,-20]$ of conductor $11$, and set $G(z)=F(z)+F(2z)$. The exact $\mathcal S_1$ images of the three basis vectors have the following coefficients in the oldform directions.

| Genus | $G(z)$ | $G(7z)$ | $G(49z)$ |
| --- | ---: | ---: | ---: |
| Positive | 8 | 58 | 98 |
| Negative | 4 | 40 | 0 |
| Negative | 2 | 2 | 0 |

The determinant is $-7056$. Each lift identity is checked through index $672$, the weight-two Sturm bound at $M=2156$. All three lifts have the same eigenvalues at primes away from $2,7,11$. Equivariance and injectivity of $\mathcal S_1$ on this eigenspace give the same good-prime eigensystem on the source. There is no square-null direction in this packet, even after mixing the two genera. The square counts retain independent information in the oldform directions at $7$. [Isogeny class 11.a](https://www.lmfdb.org/EllipticCurve/Q/11/a/)

In general, within a Hecke eigensystem with source spaces $V_+,V_-$ in the two genera, the number of new cross-genus relations after quotienting out the two within-genus kernels is
$$
\dim\bigl(\mathcal S_1(V_+)\cap\mathcal S_1(V_-)\bigr).
$$
Indeed, a relation gives a common image by taking opposite signs on the two components. The kernel of this map is exactly the sum of the two within-genus kernels. Matching good-prime eigenvalues identifies an ambient oldspace. The intersection depends on the actual lifted subspaces inside it.

## Exact replay

The [prime-power verifier](observe.py) constructs each closed neighbour graph from the scalar-part congruence forms. The [framework controls](controls.py) certify the two full genus enumerations used for the eigensystems, the collision at $p=17$ and both modular-form calculations.

For the full genus certificates, set $N=\det(L)$ and take an integral Gram matrix in a size-reduced LLL basis with Lovász parameter $3/4$. If the squared Gram–Schmidt lengths are $a,d_2,d_3$, then $d_2\ge a/2$, $d_3\ge d_2/2$ and $ad_2d_3=N$. Hence $a^3\le8N$ and $d_2\le\sqrt{2N/a}$. Size reduction bounds the second diagonal entry by $d_2+a/4$ and both inner products with the first basis vector by $a/2$ in absolute value. It also bounds the last cross term once those entries are fixed. The determinant equation determines the remaining diagonal entry. The verifier applies the quadratic-form test modulo $4$ and compares the primary discriminant forms of the even lattice $2L$. Signature and discriminant form characterize the genus. Every matching form is exactly isometric to a member of the neighbour closure. [Genus criterion](https://home.hiroshima-u.ac.jp/ichiro-shimada/ronsetsu/shimada_kobe.pdf)

~~~sh
python -m pip install cypari2 numpy
study=research/elliptic-quotients/hecke-observation
python -B "$study/observe.py" --verify "$study/prime-powers.json"
python -B "$study/controls.py" | diff - "$study/controls.json"
~~~

The original ten pair-module tests in the supplied commit `9d0fff6fdb6408edd4aedab484bb9be04552b7a8` were also reproduced. The [joint square-kernel study](../square-spectrum-kernel/README.md) gives the all-degree counts of relations. For a uniform prime-power result, the open issue is whether every inert level admits a good auxiliary prime whose norm-one profile separates all of its gluing classes.
