# Local square classes and full-spectrum collisions

The determinant-200,000 sweep found 121 opposite-secant pairs whose quotient spectra coincide and whose signed parent-neighbour theta sum is nonzero. For 90 of them, every nonzero ordinary coefficient is at a nonsquare in $\mathbb Q_2$. For 18 the same statement holds in $\mathbb Q_5$. These are exact support exclusions, including every valuation layer. The remaining 13 have a nonzero coefficient at an index that is a square at every prime dividing their theta level.

Two cases in that remaining group have an exact central-value explanation. Their signed sums belong to the weight-two packet of the elliptic curve 232.b. The discriminant −232 twist has central value zero and root number +1. A twist in a square class indistinguishable at the bad primes has nonzero central value.

The [sublattice criterion](local-geometry.md) gives an independent geometric proof for all 18 five-adic cases and 37 two-adic cases. Its isometries persist in three infinite families.

## A finite test for exclusion of local squares

Let $h=\sum_{n\ge1}a(n)q^n$ be the signed sum of the four parent-neighbour theta series. It is a weight-$3/2$ cusp form since the four lattices are in one genus. If $G$ is their parent Gram matrix, its character is $\chi=(\det G/\cdot)$ and its level $N$ is the least positive integer for which $N(2G)^{-1}$ is integral with even diagonal. Each theta coefficient counts both signs. [Half-integral-weight conventions](https://sites.lsa.umich.edu/kartikp/wp-content/uploads/sites/1326/2024/07/ssw.pdf)

For $p$ prime, write $U=U(p^2)$, so $[q^n]U^j h=a(p^{2j}n)$. Use level $L=\operatorname{lcm}(N,4p)$ for its cyclic span. The half-integral $U$ level formula shows that this level is stable and that the character is unchanged on units. Let $P_p$ retain the unit-square classes. For two it keeps just $n\equiv1\pmod8$. For odd $p$ it keeps the nonzero quadratic residues modulo $p$. Then
$$
 a(n)=0\text{ for every positive }n\in\mathbb Q_p^{\times 2}
 \quad\Longleftrightarrow\quad
 P_pU^j h=0\quad\text{for all }j\ge0.
$$
Indeed a local square has even valuation, and removing the corresponding power of $p^2$ leaves a unit square. [The half-integral $U$ and $V$ level formulas, Bringmann–Han–Heim–Kane, Lemma 2.1](https://doi.org/10.4064/aa250119-3-7)

If $m=8$ for $p=2$ and $m=p$ otherwise, the projected form has level dividing $L'=\operatorname{lcm}(L,m^2,mf_\chi)$, where $f_\chi$ is the conductor. Here is a direct proof for these projections. The sieve is a linear combination of translates $h(z+r/m)$, with coefficients invariant under multiplying $r$ by a unit square modulo $m$. Let $\gamma=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)\in\Gamma_0(L')$. Write $\tau_t(z)=z+t$, choose $s\equiv rd^2\pmod m$, and factor $\tau_{r/m}\gamma=\gamma_r\tau_{s/m}$. The matrix $\gamma_r$ is integral and has lower row $(c,d-sc/m)$. Its character and theta multiplier equal those of $\gamma$ after translating the argument. In fact $d-sc/m\equiv d\pmod{f_\chi}$, and the conductor of the quadratic character $(c/\cdot)$, as well as four, divides $c/m$ when $c\ne0$. For odd primes that divide $m$, this uses $m^2\mid c$. At two it follows from $4\mid c$ for odd $m$, and from $64\mid c$ for $m=8$. The case $c=0$ is unchanged. The translate permutation preserves the sieve weights. Rational translates also preserve holomorphy at cusps after passage to a common congruence subgroup. [Theta multiplier convention](https://sites.lsa.umich.edu/kartikp/wp-content/uploads/sites/1326/2024/07/ssw.pdf)

At either level $K$, the coefficient bound
$$
 B(K)=\left\lfloor\frac{[\mathrm{SL}_2(\mathbb Z):\Gamma_0(K)]}{6}\right\rfloor
$$
suffices for exact equality of two weight-$3/2$ forms with the same character. Multiply their difference by the weight-$1/2$ theta series with constant one, and use the weight-two order bound. Since $U^j h$ all lie in a single finite-dimensional space, compute consecutive powers until reaching the first linear relation of degree $r$. Verifying it through $B(L)$ proves that the first $r$ powers span the whole cyclic module. For those $r$ powers, check the projected coefficients through $B(L')$. If all projected forms vanish, this certifies the entire local exclusion. Every integer square is a local square, so the parent transfer then proves equality of the children in every quotient degree. [Order bound](https://www.wstein.org/msri06/refs/stein-book-on-modular-forms.pdf) · [Parent transfer](README.md#reduction-to-the-parent-genus)

## What the exact tests separate

| Parent-neighbour mechanism | Count |
|---|---:|
| Pairwise isometries of the two neighbour collections | 159 |
| Nonzero theta sum supported outside the square class over $\mathbb Q_2$ | 90 |
| Nonzero theta sum supported outside the square class over $\mathbb Q_5$ | 18 |
| Nonzero at an index that is simultaneously square over all level primes | 13 |
| Full-spectrum collisions in the sweep | 280 |

At two the maximal cyclic dimension is five. Every observed $U(4)$ orbit in this group is killed by $U(4)^4(U(4)+1)$. At five all 18 orbits are killed by $U(25)$ itself. These orbit bounds are statements about this census. The finite criterion applies at arbitrary levels with the stated theta hypotheses.

The last 13 all have explicit witnesses in the receipt. For example, at parent determinant 928 the coefficient at norm 33 is eight. The norm is a unit square at both 2 and 29, the full set of level primes. Thus these 13 cannot be explained by support entirely outside the local square class at a level prime, or by excluding the joint square class at all level primes. Their full-square equality is certified independently by the parent Shimura–Sturm test. [Exact data](local-square.json)

## A central-value explanation at parent determinants 928 and 3712

The parent Gram matrix is
$$
 \begin{pmatrix}8&0&2\\0&8&4\\2&4&17\end{pmatrix}.
$$
Its signed neighbour sum $h$ has theta level 928 and character of discriminant 232. Exact half-integral Hecke checks through norm 240 give its $T(3^2),T(5^2),T(17^2)$ eigenvalues as $(1,1,0)$. The $T(3^2)$ eigenvalue excludes unary theta components. Among the weight-two cuspidal newforms at levels dividing 464, the three eigenvalues isolate the packet of
$$
 E:\quad y^2=x^3+x^2-80x-304,
$$
of conductor 232. The square-class-33 Shimura lift of $h$ is exactly eight times the normalized newform of $E$, checked through the level-464 bound of 120. At parent determinant $3712=4\cdot928$, the theta level is 1856 and the same fingerprint isolates the same packet among newform levels dividing 928. Its square-class-33 lift is $8F(z)+8F(2z)$ for the normalized newform $F$ of $E$, verified through the level-928 bound of 240. [The elliptic curve 232.b1](https://www.lmfdb.org/EllipticCurve/Q/232/b/)

Write $\chi_0=\chi(-1/\cdot)$. The Waldspurger central-value factor for squarefree class $t$ uses the quadratic discriminant $\operatorname{disc}(\mathbb Q(\sqrt{-928t}))$. The two determinant square classes are equal. The discriminant is −232 for $t=1$ and −7656 for $t=33$. The twisted sums from the odd modular symbol of $E$ are exactly 0 and 128, respectively, in the replay's orientation. Consequently the first central value vanishes and the second is nonzero. Both twists have sign +1, so the first zero has order at least two. To apply Waldspurger with its auxiliary hypothesis at two, use the order-four character $\psi$ modulo 16 with $\psi(3)=i$ and $\psi(-1)=1$. Its square is the quadratic character of discriminant 8, and $E^{(8)}$ has conductor 1856, divisible by 16. For each $j\ge0$, twist $U(4)^j h$ in its Fourier coefficients by $\psi$. The half-integral form has character $\chi\psi^2$ and, if nonzero, belongs to the packet of $E^{(8)}$. The character transformation follows by using the Gauss-sum weights $\overline\psi(r)$ in the translate argument above. At good primes the Hecke formula gives the new eigenvalues $\psi(p)^2a_p(E)$. The trivial integral character meets (H1) and the new conductor meets (H2)(c). The two twists cancel in the central-value factor:
$$
 (E\otimes\psi^2)\otimes(\chi_0\psi^2)^{-1}=E\otimes\chi_0^{-1}.
$$
The factor therefore vanishes for squarefree part one, and Waldspurger's coefficient formula gives zero at every odd square of $(U(4)^j h)\otimes\psi$. The value of $\psi$ at any odd square is nonzero, so this proves the vanishing for $h$ at every square, with all powers of four included. The check uses rational modular symbols and does not assert an exact analytic order. [Waldspurger's hypotheses and factorization, Purkait, Theorem 4.3.4](https://wrap.warwick.ac.uk/id/eprint/50236/1/WRAP_THESIS_Purkait_2012.pdf) · [PARI exact twisted-symbol formula](https://pari.math.u-bordeaux.fr/dochtml/html/Modular_symbols.html)

There is also an exact relation between the two residual cases at parent determinant 7168. If $h_1,h_2$ are their signed sums in receipt order and $S_{8,a}$ keeps indices congruent to $a$ modulo eight, then
$$
 h_2=(S_{8,5}-S_{8,1})h_1.
$$
The ordinary coefficients through the weight-two order bound certify the relation. This expresses the second case in terms of the first. The subsequent Hecke decomposition completes the arithmetic description of the remaining eleven cases. [Residual packets and exact replay](residual-packets.md)

## Replay

With the original census dependencies installed, run from the repository root.

~~~sh
study=research/elliptic-quotients/opposite-secant-spectra/delay
python -B "$study/local_square.py" --verify "$study/local-square.json"
python -B "$study/census.py" --determinant-limit 200000 \
  --verify "$study/census-200000.json" --verify-square-null "$study/local-square.json"
~~~

The first command reconstructs the neighbour forms and verifies the local recurrences, projection bounds, simultaneous local-square witnesses, the progression identity and exact modular symbols. The second checks that the 121 inputs are precisely the original non-isometry cases in the exhaustive sweep.
