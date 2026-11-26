# Density of prime-power separators

A **separating prime** $p$ at an inert gluing level $(D,\ell)$ is a prime $p\nmid 2D\ell$ for which the quotient degrees $1,p,p^2,\ldots$ distinguish all the [refined Humbert forms](../inert-gluing/README.md) at that level. The set of separating primes has a rational natural density. If one separating prime exists, the density is strictly positive. This gives unconditional positive density at all [36 levels with a finite separation certificate](../hecke-observation/prime-powers.json).

The stronger certificates here give **density one at 13 levels** and **density exactly one-half at $D=-3,\ell=11$**. At the latter level every good auxiliary prime split in $\mathbb Q(\sqrt{-11})$ fails on the same pair. Almost every inert auxiliary prime separates all forms there.

## Observation by weight-two eigenvalues

Let $\theta_i$ be the ordinary theta series of a gluing form, and let $H_i=\mathcal S_1\theta_i$ be its square-class-$1$ Shimura lift, including the Eisenstein part. The half-integral theta level is $M=4|D|\ell^2$ and the theta character on its units is $(|D|/\cdot)$. The full Shimura lift is holomorphic of weight two with trivial character at level $N=M/2$. Its positive coefficient at $m$ is
$$
b_i(m)=\sum_{\substack{d\mid m\\(d,M)=1}}\left(\frac Dd\right)r_i((m/d)^2).
$$
At any fixed good $p$, the triangular relation and content subtraction identify equality of $b_i(p^j)$ for all $j\ge0$ with equality of all prime-power quotient counts. The appropriate neighbour formula for general $n$ uses the term $+\left(\frac{Dn}{p}\right)r_i(n)$, as derived in the [prime-power recurrence](../hecke-observation/README.md#from-neighbours-to-prime-power-counts). [Full Shimura lift in weight $3/2$](https://doi.org/10.1017/S0027763000020468)

Decompose the weight-two space into cusp newforms and Eisenstein series with their dilations. Away from the level, positive coefficients of a dilation of index greater than one vanish. Thus there are finitely many normalized eigensystems $f$ and algebraic coefficients $c_{if}$ with
$$
b_i(p^j)=\sum_f c_{if} Q_j(a_p(f),p),\qquad
Q_0=1,\quad Q_1=T,\quad Q_{j+1}=T Q_j-pQ_{j-1}.
$$
Terms with the same good-prime eigensystem are combined. Each $Q_j(T,p)$ is monic of degree $j$ in $T$. Vandermonde inversion says that at $p$ a pair $i,j$ is tied at all powers exactly when the sum of $c_{if}-c_{jf}$ over each group of equal $a_p(f)$ is zero.

## From one prime to positive density

For every cusp newform use its two-dimensional Galois representation, and for an Eisenstein series associated to $(\xi,\xi^{-1})$ use $\xi\oplus \xi^{-1}\varepsilon$, where $\varepsilon$ is the cyclotomic character. These have Frobenius trace $a_p(f)$ and determinant $p$ at good primes away from the auxiliary representation characteristic. Work over a finite extension of a local coefficient field containing all eigenvalues and expansion coefficients. [Galois representations of newforms](https://math.stanford.edu/~conrad/modseminar/pdf/L01.pdf)

Let $s$ be the number of eigensystems and $G$ the algebraic monodromy group of their direct sum. For each pair of forms the equations
$$
\sum_f (c_{if}-c_{jf})Q_k(\operatorname{tr}g_f,\det g_f)=0
\qquad(0\le k<s)
$$
define a closed conjugacy-stable locus. Take the union $Z$ over all pairs. The Vandermonde argument proves that a good Frobenius lies in $Z$ exactly when its prime fails to separate. Algebraic Chebotarev gives the separating density as the fraction of connected components of $G$ **that are not contained in $Z$**. If $p_0$ separates, choose the representation characteristic away from $p_0$ so that its Frobenius exhibits a point of one such component. This proves strictly positive density, with no twist or CM hypotheses. [Algebraic Chebotarev, Theorem 5](https://vijaympatankar.weebly.com/uploads/6/7/7/1/6771298/loc-pot-eq-gal-repr-_-ff-ell-curves-rajan-manisha-vijay-.pdf)

One can also certify a finite Chebotarev class of separating primes. Choose a residual place away from the level, $p_0$, the denominators of the expansion coefficients, and at least one nonzero coefficient witness for each pair at $p_0$. Matching joint residual Frobenius with $p_0$ preserves the witnesses modulo that place.

## Restricting quadratic twists by the conductor

Suppose two trivial-character cusp newforms of levels dividing $N$ are twists by a nontrivial character $\psi$. Their determinants force $\psi^2=1$. At every prime $q$ the conductor exponents satisfy
$$
a(\pi\otimes\psi)\le \max(a(\pi),2a(\psi)),\qquad
a(\pi)<2a(\psi)\ \Longrightarrow\ a(\pi\otimes\psi)=2a(\psi)
$$
for a ramified quadratic $\psi$. To see this, work with the local Weil–Deligne representation. Its determinant is unramified. On each inertia ramification group the invariant subspace of the finite inertia representation can have dimension only zero or two, since a one-dimensional invariant line forces the quotient to be trivial by the determinant, and the finite inertia action is semisimple. A ramified representation therefore has just one inertia break and Artin conductor twice one plus that break. The monodromy correction is zero in the ramified case. An unramified inertia representation has conductor zero or one after the monodromy correction. Twisting cannot increase the break past the maximum, and it attains the character's break when that break is strictly higher. This proves the inequalities, including tame characters.

Consequently $\operatorname{cond}(\psi)^2\mid N$. Ramification away from $2,\ell$ is impossible because $D$ is fundamental. The only possible nontrivial character for $D\ne -8$ in the class-number-one list is $\chi_{\ell^*}$, with $\ell^*=(-1)^{(\ell-1)/2}\ell$. At $D=-8$, the level also allows $\chi_{-4}$ and its product with $\chi_{\ell^*}$. The same conductor restriction applies to the quadratic self-twist of a CM newform. In particular no cusp newform in the lift can have CM by the original gluing field, since $|D|^2$ does not divide $N$.

The proof also applies to collisions among the Eisenstein systems. The weights $a_p$ that occur in the lifted theta series are real because the ternary neighbour matrices are self-adjoint for the reciprocal automorphism weights. For an Eisenstein system, real $a_p=\xi(p)+p\xi(p)^{-1}$ forces $\xi(p)\in\{\pm1\}$, so $\xi$ is quadratic. Its conductor square divides $N$. Coincidences among such Eisenstein systems are governed by the signs of their ratios. The bound $|a_p(f)|\le2\sqrt p$ for cusp forms separates them from those Eisenstein systems.

The test [in code](certify.py) uses selected good-prime neighbour matrices $A_p$ on a common finite closure. They act through the corresponding $T_p$ on the square lifts. For a proposed character $\psi$, a twist would give a common eigenvalue of $X=\sum_p w_p A_p$ and $Y=\sum_p w_p\psi(p)A_p$ for any integer weights. Coprime characteristic polynomials exclude the character, even if the closure contains classes outside the gluing family. A common eigenvalue leaves the question open. To exclude CM by a possible imaginary quadratic character $\eta$, the test checks whether $\sum_{\eta(p)=-1}A_p^2$ is invertible. A CM eigensystem would lie in its zero eigenspace.

## Promoting a finite separation test to a density statement

Outside a density-zero set, two distinct non-CM cusp systems coincide at a good prime only when a true quadratic twist relating them has value $+1$ there. This follows from the positive-density trace theorem of Patankar and Rajan and the density-zero vanishing of a non-CM trace. Their connected-monodromy hypothesis is essential for this use, so the algorithm handles possible CM contributions separately. [Positive-density trace theorem and trace vanishing](https://doi.org/10.1016/j.jnt.2017.03.001)

On one sign sector of the surviving quadratic characters, the generic partition of non-CM and Eisenstein systems by eigenvalue is fixed. If a finite prime with those signs separates all pairs, that finer generic partition also separates them. If CM has not been excluded, instead test the Krylov columns $A_p^k Bv_0$, where $v_0$ gives the norm-one counts and
$$
B=\prod_{\eta\text{ still possible CM}}\left(\sum_{\eta(q)=-1}A_q^2\right).
$$
The filter kills all actual CM systems. On each generic group of equal non-CM or Eisenstein eigenvalues its eigenvalue is constant, because every factor involves squares of eigenvalues. If the filtered observation distinguishes a pair, some non-CM or Eisenstein group has nonzero total unfiltered weight for that pair. Almost surely no CM system coincides with any such group. For non-CM cusp systems this follows from the same positive-density trace theorem, since a twist of a non-CM system cannot have CM. Hence that group certifies separation of the original pair at almost every prime in the sign sector. The sectors have equal density by the finite quadratic character group.

| Claimed density | Exact levels |
| --- | --- |
| $1$ | $(-3,5),(-4,7),(-7,5),(-8,5),(-8,7),(-11,7),(-4,11),(-8,13),(-7,13),(-3,17),(-7,17),(-4,19),(-4,23)$ |
| $1/2$ | $(-3,11)$ |

For a fixed sign sector all computations use exact integer matrices and exact separation certificates. At $(-4,11)$ and $(-4,19)$ the filtered certificate works in both signs without a claim about the existence of CM in the lift. At $D=-8$ the exact tests exclude the extra characters at $\ell=5,7,13,23,47$. The last two tests use coprime characteristic polynomials modulo $65\,521$. A common factor over the rationals would have a monic integral common factor whose degree survives reduction. See the [density receipt](receipt.json) and [larger twist checks](two-adic.json).

## A genuine half-density obstruction

At $D=-3,\ell=11$, put $\omega=(1+\sqrt{-3})/2$. The scalar parameters $\omega,2+2\omega$ give distinct forms in the same genus. Their square-lift difference is $F(z)+F(2z)$ for an exact $F\in S_2(\Gamma_0(363))$ satisfying
$$
F\otimes\chi_{-11}=-F.
$$
To check the lift identity, compute square coefficients and reconstruct $F$ in an exact cusp space. Check equality through the weight-two Sturm bound $264$ at level $726$. To check the twist identity, compare the coefficients in the same level $363$. Indeed the primitive $\chi_{-11}$ twist of any newform with $11$-exponent at most two has exponent at most two by the local bound. Depleting the twisted form at $11$ preserves that bound. For exponent zero use $f-a_{11}f(11z)+11f(121z)$, for exponent one use $f-a_{11}f(11z)$, and for exponent two the coefficients at multiples of $11$ already vanish. This proves the Fourier twist acts within level $363$ on the full cusp space there. The replay checks the twist coefficients through index $264$, exceeding the bound $88$ at level $363$. [Sturm bound](https://www.lmfdb.org/knowledge/show/cmf.sturm_bound)

If $\chi_{-11}(p)=1$ and $p$ is good, the coefficients of $F(z)+F(2z)$ at every $p^j$ vanish. The two geometric prime-power spectra are therefore equal. On the negative sign sector the verified separating prime $7$ certifies separation generically. The density is exactly $1/2$.

A collection of coincidence characters also has a simple consistency test. There is an element on which every character has value $-1$ exactly when there is no odd-cardinality product of the characters equal to $1$. The smallest obstruction is a dependent triple $\psi_1,\psi_2,\psi_1\psi_2$.

## Reproduction

~~~sh
python -m pip install cypari2 numpy python-flint
study=research/elliptic-quotients/prime-density
python -B "$study/certify.py" --verify "$study/receipt.json"
python -B "$study/two_adic.py" | diff - "$study/two-adic.json"
~~~

The integer Krylov replay also checks the all-powers collisions at $p=17$ for $(-3,5)$ and at $p=17,29$ for $(-8,7)$. It finds $2$, $11$, and $3$ tied pairs respectively.

The level-wide existence problem remains open. The results above make its consequence unconditional whenever one prime is found, and give an actual level where half of all auxiliary primes fail for the same arithmetic reason.
