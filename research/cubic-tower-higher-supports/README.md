# Higher factors in the cubic tower

The genus-393,661 tenth-stage Jacobian has at least **2,047 distinct simple isogeny classes over $\mathbb Q$**. One can choose a different simple constituent in each of its 2,047 signed subset factors, including the elliptic base. The previous computation identified the 21 elliptic and 90 four-dimensional factors and proved their geometric simplicity. The new bound gives at least **1,936 additional rational simple classes**, with one representative in each higher signed support. Their dimensions, multiplicities and possible identifications over $\overline{\mathbb Q}$ remain open. [Earlier factor construction](../cubic-tower-factors/README.md) · [Rank certificate](rank.json)

A separate calculation gives the full genus-55 third-stage zeta numerator modulo 3 at 11 and both degree-24 Frobenius polynomials modulo 3 for all 120 three-index supports at 41. The underlying congruence applies at every stage and every simultaneous good prime $p>3$. [Coefficients and point counts](mod-three.json)

Subsequent computations determine the geometric simple factors through the [seventh tower stage](../cubic-tower-stage-seven/README.md) and give new integral coefficients for genus 55. [Cartier matrices and zeta data](../cubic-tower-cartier/README.md) · [Geometric separation](../cubic-tower-geometric-separation/README.md) · [Five- and six-index computations](../cubic-tower-size-five/README.md)

## Frobenius rank over the rationals

Take the factors $B_j$ to be the base and the two factors $A_S^\pm$ for each nonempty subset of $\{0,\ldots,9\}$. At a simultaneous good prime $p$, put $M_{pj}=\operatorname{Tr}(F_p\mid H^1(B_j))$. The 2,055 rows come from the first simultaneous good primes in increasing order, from 41 through 18,329. The first column is the elliptic base. For each increasing nonzero binary mask, with bit $i$ representing index $i$, the positive and negative traces follow, in that order. The matrix has rank **2,047 modulo 127**. Hence it also has 2,047 independent columns over $\mathbb Q$. [Portable rank verifier](prime_rank.py)

For the implication about isogenies, write $B_j\sim_{\mathbb Q}\prod_b C_b^{m_{bj}}$ with the $C_b$ pairwise nonisogenous and simple over $\mathbb Q$. Each $C_b$ has good reduction at the witness primes because its rational Tate module is a direct summand of that of a good Jacobian. Additivity of the trace factors $M$ through the integer multiplicity matrix $(m_{bj})$. The latter therefore has full column rank. Its row count bounds the number of distinct simple classes from below. A nonzero maximal minor also has a nonzero permutation term. This selects a different constituent with $m_{bj}>0$ in every column. The original 111 simple columns use their own 111 rows in every such selection. Restricting to the columns of any $k$-index stage inside stage ten gives the bound $2^{k+1}-1$ there. These arguments concern isogenies defined over $\mathbb Q$. A rationally simple class can split or become isogenous to another class after extending the base field.

The good-prime sieve tests the discriminant of $f=-4z^3+z^2+z+1$, the leading coefficients and discriminants of the quadratic branch polynomials $D_i$, and all required pairwise resultants. For each prime, the counter computes the ordinary and sign-twisted fibre sums proved in the [tower construction](../cubic-tower-factors/README.md). Splitting the ten coordinates into two groups of five lets one matrix product sum the 1,024 subset products. The absolute sum for each dot product is at most $(2p+1)2^{10}$, so the binary64 operations are exact on these integers. Rank reduction takes Schur complements modulo 127 in blocks of 32, with integer intermediates bounded by $32\cdot126^2$. An independent counter checks all 2,047 traces at each of the first four primes. A second implementation in C++ uses integer arithmetic, individual subset products and ordinary modular Gaussian elimination. It reproduces the matrix digest and the rank independently. [C++ verifier](trace_rank.cpp)

## Quadratic resolvents modulo 3

Let $P_{S,+}$ and $P_{S,-}$ denote the local numerators of the signed factors for nonempty $S$, with the convention $P(T)=\det(1-TF_p)$. Write

$$
D_i=(12i+1)z^2+(12i^2+1)z+4i^3+1,\qquad h_U=\prod_{i\in U}(-27D_i).
$$

For nonempty $U$ let $H_{U,b}$ be the zeta numerator of the smooth projective model of $y^2=f^b h_U$, with $b=0$ or $1$. Its genus is $|U|-1+2b$. Put $P_E=H_{\varnothing,1}$. Then in $\mathbb F_3(T)$,

$$
P_{S,+}(T)=\frac{\prod_{\varnothing\ne U\subseteq S} H_{U,\,|U|\bmod2}(T)}{(1-T)(1-pT)},\qquad
P_{S,-}(T)=P_E(T)\prod_{\varnothing\ne U\subseteq S} H_{U,\,1-(|U|\bmod2)}(T).
$$

The quotient on the left is a polynomial after reduction. These formulas replace the dimension $|S|2^{|S|-1}$ of each signed factor by auxiliary curves of genus at most $|S|+1$.

**Proof.** Over $E:w^2=f(z)$ the common splitting field on support $S$ has geometric group $S_3^{|S|}$. Extending the involution $w\mapsto-w$ by sending every labeled root of each cubic to its negative gives a central involution $\sigma$ over the $z$-line. Negating the roots preserves the defining symmetric relations with $w$ negated. The disjoint splitting fields give the independent labels in each cubic. The resulting group over $\mathbb P^1_z$ is $S_3^{|S|}\times\langle\sigma\rangle$. Write $\eta$ for the character of the second factor, corresponding to $\sqrt f$.

The permutation representation on the sheets of $X_S$ is the regular representation of $\langle\sigma\rangle$ tensored with the three-letter permutation representations of each $S_3$. Forgetting an index projects its permutation representation onto the trivial summand. The complementary summand is the standard two-dimensional representation $V_i$. Thus the Artin $L$-function of $(\bigotimes_{i\in S}V_i)\otimes\eta^\delta$ is the signed numerator, for $\delta=0$ on the positive sign and $\delta=1$ on the negative sign. Indeed, the trivial and permutation representations give the base and cover zeta functions by counting Frobenius orbits on the inertia orbits of sheets at every closed point. Decomposing by subsets and the central involution gives the assertion. The degree-zero and degree-two terms occur only for the empty support and the trivial character.

Choose the augmentation lattice in $\mathbb Z_3^3$ for $V_i$. Modulo 3 the diagonal line lies in the augmentation plane and is fixed. Its one-dimensional quotient has the sign character, as follows from the determinant. A tensor-product filtration therefore has one character $\eta^\delta\prod_{i\in U}\operatorname{sgn}_i$ for every $U\subseteq S$. Inertia has order one or two at every point. At each $D_i$ branch value a single cubic ramifies. Above the roots of $f$ and infinity the splitting extension over $E$ is unramified and only $E\to\mathbb P^1$ ramifies. Averaging over inertia is consequently an integral projector at 3. Reduction preserves its invariant summand and taking inertia invariants preserves the filtration. Frobenius preserves the resulting filtration at each closed point, so the Euler determinant on the original lattice is the product of the character determinants modulo 3. This proves the coefficientwise congruence for entire Euler products, including terms of degree divisible by 3.

The Vandermonde $v_i$ of the three cubic roots satisfies $v_i^2=-27D_i$. Cubic permutations act on $v_i$ by their sign, while $\sigma$ negates all three root differences and sends $v_i$ to $-v_i$. Therefore $wv_i$ realizes the sign character trivial on $\sigma$, with square $f(-27D_i)$. The product character for $U$ has square class $h_U f^{(|U|+\delta)\bmod2}$. Its $L$-function is the numerator of the corresponding double cover since its product with the trivial $L$-function is the double-cover zeta function. For $U=\varnothing,\delta=0$ the $L$-function itself is $1/((1-T)(1-pT))$. The other empty term is $P_E$. This gives both formulas.

## Étale cyclic covers realizing the signed factors

The same splitting field gives geometric models over $\mathbb Q$ for the square of every signed factor. Fix a nonempty $S$ of size $k$, let $Y$ be the smooth projective common splitting cover over the $z$-line, and put $G=S_3^k\times\langle\sigma\rangle$. All its deck transformations are defined over $\mathbb Q$. Indeed the independent cubic permutations and simultaneous negation act on the rational splitting-field presentation above. The geometric group already realizes the full degree $2\cdot6^k$, so the splitting field has constant field $\mathbb Q$.

Put $H=A_3^k\times\langle\sigma\rangle$. The curve $K=Y/H$ has function field obtained by adjoining the $k$ square roots of $f(-27D_i)$ to $\mathbb Q(z)$. These square classes are independent even geometrically, since each has two branch values unique to $D_i$. Its cover of the line has degree $2^k$ and inertia of order two at $2k+4$ values, the roots of the $D_i$ and $f$, and infinity. Riemann–Hurwitz gives $g(K)=1+k2^{k-1}$. Every nontrivial inertia group of $Y$ has order two and maps nontrivially to $G/H$. At a root of $D_i$ the $i$th sign is ramified. At the roots of $f$ and infinity every $f(-27D_i)$ has odd valuation. Thus $Y\to K$ is everywhere étale.

Choose a nontrivial cubic character $\chi_i$ of each $A_3$ and define $\psi_\delta=(\prod_i\chi_i)\eta^\delta$ on $H$, for $\delta=0,1$. The quotient $D_\delta=Y/\ker\psi_\delta$ is a connected étale cyclic cover of $K$ of degree three or six, respectively. The kernels are subgroups of the rational deck group, so the curves and their deck actions descend over $\mathbb Q$ even though the character values use roots of unity. The degree-three quotient of $D_1$ is $D_0$, and its degree-two quotient is $L=K(\sqrt f)$. Indeed the kernels of the order-three and order-two components intersect in $\ker\psi_1$. Define the faithful part of $J(D_1)$ as the abelian image of the rational group-algebra projector onto the two order-six characters, equivalently the identity component of the common kernel of the norm maps to $D_0$ and $L$. Then

$$
\operatorname{Prym}(D_0/K)\sim_{\mathbb Q}(A_S^+)^2,\qquad J(D_1)_{\mathrm{faithful}}\sim_{\mathbb Q}(A_S^-)^2.
$$

To see the multiplicity two, let $R_\delta$ be the faithful part of the rational regular representation of $H/\ker\psi_\delta$. Over a splitting field its two characters are $\psi_\delta$ and its inverse. Each induces to $W_\delta=(\bigotimes_{i\in S}V_i)\otimes\eta^\delta$ in $G$, because either nontrivial character of $A_3$ induces to the standard representation of $S_3$. Thus $\operatorname{Ind}_H^G R_\delta\simeq W_\delta^{\oplus2}$ over $\mathbb Q$, by equality of rational representation characters. For a subgroup $U$ and $V=H^1(Y_{\overline{\mathbb Q}},\mathbb Q_\ell)$, pullback gives $H^1(Y/U)=V^U=\operatorname{Hom}_G(\mathbb Q_\ell[G/U],V)$, equivariantly for the absolute Galois group of $\mathbb Q$. Applying the faithful projector on $D_\delta$ replaces the permutation module by $\operatorname{Ind}_H^G R_\delta$. Applying the signed subset projector on $X_S$ selects $W_\delta$. This identifies the faithful rational cohomology with two copies of that of $A_S^\pm$. The isogenies now follow from the isogeny theorem over $\mathbb Q$. [Faltings, *Endlichkeitssätze für abelsche Varietäten über Zahlkörpern*](https://doi.org/10.1007/BF01388432)

For $k=3$ the curves $K,D_0,D_1$ have genera 13, 37 and 73. The two displayed abelian images have dimension 24 each. At a simultaneous good prime the exact numerator identities are $P_{D_0}=P_K P_{S,+}^2$ and $P_{D_1}P_K=P_{D_0}P_L P_{S,-}^2$. The second equality removes the contributions from both proper nontrivial cyclic quotients, whose common contribution is $J(K)$.

## The third stage at 11

For $S=\{0,1,2\}$ define $a=T^2+T-1$, $b=T^2-T-1$, $c=T^2+1$ and $d=T^6-T^5-T-1$ in $\mathbb F_3[T]$. The two new degree-24 factors and the full degree-110 numerator are

$$
P_{S,+}=(T^2-1)^5 a b^6,\qquad P_{S,-}=(T^2-1)^4 b c^4 d,\qquad
P_{X_S}=-(T^2-1)^{26}a^2 b^{15} c^6 d^2.
$$

Six exact direct counts on the cubic tower check the two signs and the full numerator modulo 3.

| Extension degree | Points on the tower | Fixed points of sign times Frobenius | Trace on $A_S^+$ | Trace on $A_S^-$ |
|---:|---:|---:|---:|---:|
| 1 | 9 | 13 | −2 | −2 |
| 2 | 293 | 65 | −2 | −2 |
| 3 | 1,314 | 1,276 | −65 | −59 |
| 4 | 15,957 | 15,841 | 26 | 210 |
| 5 | 161,429 | 159,333 | 638 | 178 |
| 6 | 1,783,142 | 1,782,164 | −7,865 | −761 |

At 41 the 240 signed triple polynomials supply 5,760 additional numerator degrees modulo 3. Together with the 762 degrees already known integrally, this gives 6,522 of 787,322 numerator degrees at stage ten modulo 3. The verifier also compares every coefficient of the seven archived exact factors at 11 and 111 at 41 with the congruence. To count the auxiliary hyperelliptic curves it groups field elements by their quadratic-character signs and takes a Walsh transform of these counts. Zeros of a branch polynomial are counted separately. Newton identities over the integers and reciprocity give the auxiliary numerators before reduction modulo 3. [Congruence verifier](mod_three.py)

## Replay

From this directory, using Python 3.11 or later with NumPy and SymPy installed, run

```sh
python prime_rank.py
python mod_three.py
```

On the reference machine the first run took about 12 seconds and the second about a minute. Each `--compute` option prints its fully rebuilt receipt. The default replay uses Python throughout. The optional independent C++17 check below prints `2047`. The prime list is also accepted on standard input by Windows builds of the program.

```sh
c++ -O3 -std=c++17 trace_rank.cpp -o trace_rank
python prime_rank.py --primes | ./trace_rank
```
