# The zeta factors of the nonabelian cubic tower

The genus-13 second stage of the Crest–Root tower has a Jacobian isogenous over the rationals to five elliptic curves and two absolutely simple abelian fourfolds. These seven factors are pairwise nonisogenous even over an algebraic closure. We compute its full zeta function at 11. In the genus-393,661 tenth stage, the same construction exhibits **21 elliptic curves and 90 absolutely simple fourfolds**, again pairwise geometrically nonisogenous. Their 111 Frobenius polynomials at 41 occupy 762 of the 787,322 numerator degrees. The remaining supports and their dimensions are determined by the theorem below. Further simple constituents and their multiplicities inside those supports remain open. [Exact coefficients and witnesses](receipt.json)

A subsequent calculation separates all 2,047 displayed factors in the rational isogeny group. Their simple constituents can be chosen pairwise distinct over the rationals, one per displayed factor. It also computes triple-support polynomials modulo 3. [Higher supports and trace-rank proof](../cubic-tower-higher-supports/README.md)

## The tower and its factors

Write $E:w^2=f(z)=-4z^3+z^2+z+1$. For a finite set $I$ of nonnegative integers, let $X_I$ be the smooth projective normalization obtained by adjoining one root $t_i$ of each of the cubics

$$t_i^3+3(z+i)t_i=w,\qquad i\in I.$$

The empty set gives $E$. The original Part XXI, §§91–92, supplied prime-field counts through ten levels. Part XXIII, Theorem 100.1, proves geometric disjointness of the cubic closures with group $S_3^{|I|}$ and the genus formula $g_I=1+2|I|3^{|I|-1}$ for the tower. Its proof works for arbitrary finite subsets, because the discriminant supports are pairwise disjoint. These are the inputs from the historical monograph. The factorization here uses the classical correspondence-idempotent method on this specific tower. [Crest–Root research archive guide](https://github.com/towre676-cloud/crest-root-calculus#whats-in-the-research-archive) · [Kani and Rosen, *Idempotent relations and factors of Jacobians*](https://doi.org/10.1007/BF01442878)

**Theorem.** Over $\mathbb Q$ there are abelian varieties $A_S^+, A_S^-$ for every nonempty subset $S\subset I$ (including $S=I$), for which

$$
J(X_I)\ \sim\ E\times
\prod_{\varnothing\ne S\subseteq I}(A_S^+\times A_S^-),
\qquad
\dim A_S^+=\dim A_S^-=|S|2^{|S|-1}.
$$

Here $\sim$ denotes isogeny. Each displayed factor is defined up to isogeny and is pulled up from $X_S$. The general formula leaves room for further decomposition of the factors. At a simultaneous good prime it also decomposes the first étale cohomology and factors the zeta numerator.

**Subset projectors.** On $J(X_I)$ define $e_i=\pi_i^*(\pi_i)_*/3$, where $\pi_i:X_I\to X_{I\setminus\{i\}}$ forgets the $i$-th cubic root. These are idempotents since pushforward after pullback is multiplication by three. They commute. On the common unramified locus they average independently over the three choices in each forgotten coordinate. Their compositions are finite correspondences, hence agree globally by taking closures. Products of the $e_i$ and $1-e_i$ give orthogonal rational idempotents summing to the identity. The projector with $1-e_i$ precisely for $i\in S$ cuts out the contribution first appearing on $X_S$.

**The two signs.** The involution $\sigma(z,w,t_i)=(z,-w,-t_i)$ commutes with the subset projectors. Over each of the three finite fixed points on $E$ its only fixed sheet is $t_i=0$ for every $i$. These roots are simple because $f(-i)\ne0$. Over the origin of $E$ one normalized sheet in each cubic is fixed and the other two are exchanged. Indeed, with $\tau=z/w$, the leading constants of $t_i\tau$ are $-1,1/2,1/2$. On the last two sheets $t_i=1/(2\tau)+b+O(\tau)$ with $b^2=-(12i+1)/12$, and $\sigma$ sends $b$ to $-b$. Thus there are exactly four fixed points on $X_I$. Riemann–Hurwitz gives the quotient genus $|I|3^{|I|-1}$ and anti-invariant dimension $1+|I|3^{|I|-1}$. For the empty set the plus part is zero and the minus part is $E$. Möbius inversion on subsets gives, for either sign and $k=|S|\ge1$,

$$
\sum_{j=1}^k(-1)^{k-j}\binom{k}{j}j3^{j-1}=k2^{k-1}.
$$

This identifies the dimensions of the images of the subset projectors multiplied by $(1\pm\sigma)/2$. Multiplication by a common denominator makes these rational endomorphisms honest homomorphisms with abelian images.

## The integral splitting away from 2 and 3

For $n=|I|$ put $D=2\cdot 3^n$. For every nonzero signed subset projector $P_j$ above, $DP_j$ is an integral endomorphism of $J=J(X_I)$. In fact each $e_i$ is a push-pull endomorphism divided by 3, each $1-e_i$ has the same denominator bound, and the two sign projectors have denominator 2. This also covers the base projector, whose nonzero image is isogenous to $E$.

Take the actual abelian images $B_j=\operatorname{im}(DP_j)\subset J$. Define $a:\prod_j B_j\to J$ by addition and $b:J\to\prod_j B_j$ by $x\mapsto(DP_jx)_j$. Then $ab=[D]$ since $\sum_jP_j=1$. On $B_i$ the endomorphism $DP_j$ is zero if $j\ne i$ and multiplication by $D$ if $j=i$. This follows by precomposing with the surjection $DP_i:J\to B_i$ and using $P_jP_i=\delta_{ij}P_i$. Consequently $ba=[D]$ as well. The addition isogeny therefore has kernel annihilated by $D$. At every prime $\ell>3$ it gives an isomorphism of integral Galois modules

$$
\bigoplus_j T_\ell(B_j)\ \xrightarrow{\ a\ }\ T_\ell(J)
$$

with inverse $D^{-1}T_\ell(b)$. This advances the integral-projector question **F11** on this tower. It supplies explicit integral multiples and a bound on the exponent of the gluing kernel. The minimal denominators, the precise 2- and 3-primary indices and the induced polarization types remain open.

## From each fibre to a local polynomial

For $p>3$ check that $f$ is squarefree, each $D_i=(12i+1)z^2+(12i^2+1)z+4i^3+1$ has degree two and is squarefree and prime to $f$, and the $D_i$ are pairwise coprime modulo $p$. The argument of Theorem 100.1 then applies in the residue characteristic: the normalized infinity is unramified, and the branch sections coming from the four points per cubic are disjoint and tame. These conditions give smooth normalized covers of the same genus and the same four involution-fixed points. The routines test them afresh for every support used at a witness prime.

Put $q=p^r$. Use $s=1$ for the ordinary fixed points and choose a nonsquare $s$ in $\mathbb F_q$ for the fixed points of $\sigma\operatorname{Frob}_q$. Writing $w=\eta w'$, $t_i=\eta t_i'$ with $\eta^2=s$ gives

$$
s(w')^2=f(z),\qquad s(t_i')^3+3(z+i)t_i'=w'.
$$

For a point $P$ of this base let $N_i(P)$ be the number of solutions of the second equation. At infinity the number is $2+\chi_q(s)\chi_q(-(12i+1)/12)$. For each nonempty support $S$, form

$$
b_S(s)=-\sum_P\prod_{i\in S}(N_i(P)-1).
$$

To see the formula, the fixed point count of each intermediate $X_U$ is the sum over base points of $\prod_{i\in U}N_i(P)$. This includes the tame branch points, where one ramified root contributes one geometric point. Disjoint branch supports ensure the normalized fibre product has the product fibres. Taking the alternating sum over $U\subseteq S$ cancels the degree-zero and degree-two cohomology terms in the trace formula and expands the product of differences. For a nonsquare $s$ this reasoning concerns the fixed points of the sign composed with Frobenius at the specified extension. Choosing a nonsquare again in an even extension is essential.

The traces of the $r$-th Frobenius power on the positive and negative factors are $(b_S(1)+b_S(s))/2$ and $(b_S(1)-b_S(s))/2$. For a factor of abelian dimension $d$, Newton identities recover the first $d$ coefficients of $\det(1-TF)$ from $r=1,\ldots,d$ and the functional equation gives the other $d$. This is a finite algorithm for a whole stage under the good-reduction hypotheses. The largest needed extension exponent is $|I|2^{|I|-1}$, which is 5,120 at stage ten. Enumerating the corresponding field is beyond the included implementation.

To count a depressed cubic efficiently, precompute the distributions of $x^3$, $x^3+x$ and $x^3+gx$ for a generator $g$ of $\mathbb F_q^\times$. For $x^3+ax=b$ with $a\ne0$, write $a=r^2c$ for $c\in\{1,g\}$. Then the desired count is the frequency of $b/r^3$ in the table for $x^3+cx$. The tables require linear storage in $q$. The field constructor verifies that its multiplicative orbit contains every nonzero residue.

## The genus-13 stage at 11

Indices 0 and 1 satisfy the simultaneous good-reduction tests at 11. Five ordinary point counts on $X_{\{0,1\}}$ were $7,203,1300,15971,163677$. The corresponding signed fixed point counts were $5,89,1370,14561,161295$. The first four extensions determine every factor; the fifth checks the prediction afterward. If $Z(T)=P(T)/((1-T)(1-11T))$, then

$$
P(T)=(1+3T+11T^2)^3(1-5T+11T^2)(1-3T+11T^2)P_+(T)P_-(T).
$$

| Coefficients by increasing power | $T^0$ | $T^1$ | $T^2$ | $T^3$ | $T^4$ | $T^5$ | $T^6$ | $T^7$ | $T^8$ |
|---|---:|---:|---:|---:|---:|---:|---:|---:|---:|
| $P_+$ | 1 | −4 | 15 | −27 | 124 | −297 | 1815 | −5324 | 14641 |
| $P_-$ | 1 | −2 | 11 | −7 | 156 | −77 | 1331 | −2662 | 14641 |

The cubic power of the same elliptic numerator at 11 illustrates what one reduction can conceal. Different witness primes separate the characteristic-zero elliptic factors, even after arbitrary scalar extension.

## Absolute simplicity and geometric separation

For any of these Frobenius polynomials, a ratio of two roots inside a single polynomial, or between two polynomials of degree at most eight, has at most 64 algebraic conjugates. If it is a root of unity of order $m$, then $\varphi(m)\le64$. Also $\varphi(m)\ge\sqrt{m/2}$ by the prime-power formula. The verifier computes $M=\operatorname{lcm}\{m\le8192:\varphi(m)\le64\}$. Any collision between powers of two roots over any finite extension forces a collision between their $M$-th powers.

For each fourfold, choose a good reduction at which its degree-eight Frobenius polynomial is irreducible over $\mathbb Q$. The checker also proves that the characteristic polynomial of its $M$-th power is squarefree modulo the auxiliary coefficient prime 1,000,003. It powers the companion matrix there and applies Newton identities to the matrix traces. Since its eight roots have distinct powers over every finite extension, their powers remain one Galois orbit there. The reduction is absolutely simple and therefore so is the characteristic-zero factor. A nontrivial isogeny factor in characteristic zero would specialize after a finite extension at a good prime. The analogous calculation for a pair finds a good reduction at which the two powered polynomials are coprime modulo 1,000,003. Their reductions cannot become isogenous after any extension. This proves geometric nonisogeny in characteristic zero, including possible twists. An elliptic curve is already absolutely simple. [Howe and Zhu on absolute simplicity over finite fields](https://arxiv.org/abs/math/0002205)

At stage ten, every singleton contributes two elliptic factors and every pair two fourfolds. Including $E$ gives 21 elliptic curves and 90 fourfolds. Of the fourfolds, 86 have their simplicity witness at 41. Prime 11 handles $(0,1)^+$ and $(0,2)^+$, and 13 handles $(5,9)^-$ and $(8,9)^-$. Of the 4,215 pairs with equal dimensions, 4,179 separate at 41. The remaining 36 use other primes and are individually listed in the receipt. Witnesses need good reduction only on the supporting smaller stages: the corresponding factors there pull back up to isogeny into stage ten over $\mathbb Q$. Distinct dimensions then separate the remaining pairs. This calculation imposes no simplicity assumption on the larger supports.

| Extension degree at 41 | Points on stage ten | Fixed points of sign composed with Frobenius |
|---:|---:|---:|
| 1 | 81 | 729 |
| 2 | 59,993 | 2,357 |
| 3 | 85,560 | 87,804 |
| 4 | 6,563,741 | 2,856,185 |

## Replay

From this directory, run the following commands with Python 3.11 or later. On Windows, `py` can replace `python`.

```sh
python -m pip install numpy sympy
python verify.py
```

The verifier rebuilds and compares the [receipt](receipt.json), including every published polynomial and the exceptional witness primes. `python verify.py --compute` prints the complete recomputation. Literal enumeration over the small fields checks the ordinary and twisted counts independently of the three cubic histograms. Separate checks reproduce all four stage-ten prime-field counts from §92 and reject the bad stage-two reduction at 5. The largest table in the default replay has $41^4=2{,}825{,}761$ entries. The reference implementation uses 64-bit integer arrays with checked accumulator bounds. Its default field-order budget is three million.
