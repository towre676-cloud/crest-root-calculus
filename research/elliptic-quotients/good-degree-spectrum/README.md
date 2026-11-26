# Coprime quotient degrees and separating-prime density

For the nine class-number-one imaginary quadratic fields, the degrees coprime to $2D\ell$ distinguish every refined Humbert form at every odd inert gluing prime $\ell\le199$. This is an exact census of **210 levels and 554,482 forms counted across the levels**. The largest last separating degree is **499**, at $D=-163,\ell=193$. The scalar parameters $34+29\omega$ and $48$ have zero and four geometric quotients of that degree.

The density of auxiliary primes whose powers distinguish every form is **one at 204 levels** and **one-half at six levels**. No prime-neighbour matrices at these larger levels are needed for the density-one certificates. They use coefficients derived from the same integer representation counts as the coprime-degree census.

## The coefficient test

Let $P_a(m)$ be the geometric quotient count and $r_a(n)$ the ordinary representation count of the Humbert lattice $L_a$. Put $M=4|D|\ell^2$, $N=M/2$. The $t=1$ Shimura lift $H_a$ of the theta series of $L_a$ lies in $M_2(\Gamma_0(N))$. For $(m,2D\ell)=1$ its coefficient is
$$
 b_a(m)=\sum_{d\mid m}\chi_D(d) r_a((m/d)^2)
       =\sum_{k\mid m}P_a(k)\sum_{d\mid m/k}\chi_D(d).
$$
The second equality groups the ordinary counts by the content of the representing vector. The equation is triangular with diagonal one, so separation by these $b_a(m)$ is equivalent to separation by the same entire set of coprime quotient degrees. [The full lift and its coefficient formula](../prime-density/README.md#observation-by-weight-two-eigenvalues)

There is a finite decision bound even for a pair that stays tied in a search. For each prime $r\mid N$ the operator $1-V_rU_r$ deletes coefficients of indices divisible by $r$, where $V_r f(z)=f(rz)$. Applied successively, these operators send forms of level $N$ to forms of level $N\operatorname{rad}(N)$. At every stage $r$ already divides the current level, so $U_r$ preserves that level and $V_r$ raises it by at most $r$. The constant coefficient is killed as well. Checking the remaining coefficients through the weight-two Sturm bound at $N\operatorname{rad}(N)$ settles equality on all coprime indices. The census uses much smaller explicit witnesses. [Sturm bound](https://www.lmfdb.org/knowledge/show/cmf.sturm_bound)

## Odd discriminants force Steinberg at the ramified prime

Every odd discriminant in the family is $D=-q$ with $q$ prime and $q\equiv3\pmod4$. Write the ambient norm as
$$
 u^2+(2x+y)^2+qy^2.
$$
If $q^2$ divides this norm, reduction mod $q$ first gives $q\mid u,2x+y$, and division by $q$ then gives $q\mid y$. The congruence defining $L_a$ is unchanged by division by $q$ because $q\ne\ell$. Consequently
$$
 r_a(q^2n)=r_a(n)\quad(n\ge1),\qquad U_q H_a=H_a.
$$
For the second identity use the lift coefficient formula at all indices, with the character extended by zero on nonunits of $M$.

At $q$ the weight-two level has exponent one. In the $q$-old part for a cusp newform, the $U_q$ polynomial is $X^2-a_qX+q$. The eigenvalue $1$ would require $a_q=q+1$, which exceeds the cuspidal bound $2\sqrt q$. Thus each cuspidal system occurring in $H_a$ is $q$-new and has $a_q=1$. With trivial central character, the local representation of conductor exponent one is Steinberg up to an unramified twist. A CM cusp form is automorphically induced from a character, with zero local monodromy, so none of the cusp constituents here has CM. [Oldform theory](https://andreeamocanu.github.io/atkin-lehner-theory.pdf)

The [conductor restriction](../prime-density/README.md#restricting-quadratic-twists-by-the-conductor) leaves only $\chi_{\ell^*}$ as a possible nontrivial twist at odd $D$. Therefore **at every odd discriminant in the family**, coprime-degree injectivity, existence of a separating auxiliary prime and separating-prime density at least one-half are equivalent. For the forward implication, on the $\chi_{\ell^*}(p)=-1$ sector the non-CM systems have pairwise distinct eigenvalues outside a density-zero set. The two reverse implications follow directly from the definitions. The argument includes Eisenstein systems, whose real good-prime eigenvalues are $(p+1)\xi(p)$ for quadratic $\xi$ of conductor square dividing $N$. [Positive-density trace theorem](https://doi.org/10.1016/j.jnt.2017.03.001)

In addition, if $\ell\equiv1\pmod4$, quadratic reciprocity and inertness give $\chi_{\ell^*}(q)=-1$. A twist would reverse the sign of $a_q$ and is excluded. The analogous Eisenstein restriction follows from $U_q=1$. At exponent one the two $q$-old roots for an Eisenstein system are $\xi(q)$ and $q\xi(q)$, so an occurring system has $\xi(q)=1$. Its twist by a character with value $-1$ at $q$ cannot also occur. In this odd-discriminant, $\ell\equiv1\pmod4$ case, coprime-degree injectivity implies density one.

## A sector certificate that allows extra twists and CM

At even $D$, the twist restriction alone leaves a second possibility. Several CM systems with the same CM field can all have eigenvalue zero at every good prime inert in that field. Coprime-degree data can distinguish their contributions even though a single inert prime cannot. The density certificate below removes this problem before testing a sector.

For $D=-4$ the only candidate is $\chi_{\ell^*}$. For $D=-8$ the conductor also allows $\chi_{-4}$ and its product with $\chi_{\ell^*}$, except where [already excluded by exact characteristic polynomials](../prime-density/two-adic.json). The only possible CM field in either even-discriminant case has character $\chi_{\ell^*}$, and only if it is imaginary. In fact, an imaginary candidate ramified at $2$ in the $D=-8$ case has quadratic discriminant exponent two there. A CM inducing character with target level exponent at most four then has conductor exponent at most two at the ramified prime above $2$, by the CM level formula. It is trivial on $-1\in 1+\mathfrak p^2$. Its restriction to $\mathbb Q_2^\times$ would have to equal the local quadratic character on units, which sends $-1$ to $-1$. This excludes both imaginary candidates ramified at $2$. [CM level and character formula](https://dermotmccarthymath.github.io/publications/CM_Qsqrt%28-3%29.pdf)

A **sign sector** specifies the values of the allowed characters at the auxiliary prime. The sector test uses only good indices $m$ satisfying both conditions below.

1. For each candidate twist $\psi$ whose specified sign at the auxiliary prime is $+1$, require $\psi(m)=1$.
2. For each possible CM character $\eta$, require a prime factor $r$ of $m$ with odd valuation and $\eta(r)=-1$.

For a pair $a,a'$, write the difference of lifted coefficients at good indices as a linear combination of normalized newform and Eisenstein coefficients. The second condition makes every CM term zero, since a CM newform has zero odd-power coefficients at an inert prime. Outside a density-zero set in the prescribed sector, the generic groups of equal good-prime eigenvalues among the remaining systems are exactly the appropriate true twist groups and Eisenstein groups. A non-CM system and a CM system coincide only on a density-zero set. The first condition makes $a_m(f)$ constant on each generic group, by the twist formula $a_m(f\otimes\psi)=\psi(m)a_m(f)$. A nonzero difference of $b$ at one admissible $m$ therefore forces a nonzero original weight sum in some non-CM or Eisenstein group. The prime-power Vandermonde criterion then separates the pair at almost every prime in that sector. A finite list of indices separating all pairs certifies the entire sector up to density zero. [Vandermonde criterion](../prime-density/README.md#observation-by-weight-two-eigenvalues) · [Trace coincidence theorem](https://doi.org/10.1016/j.jnt.2017.03.001)

This test does not require identifying the true twists among the candidates. Equally sized sign sectors partition the good primes up to finitely many primes. Certifying all sectors gives density one.

## Where exactly half the primes work

The six levels are
$$
(-3,11),\quad(-7,3),\quad(-19,3),\quad(-43,3),\quad(-67,3),\quad(-163,3).
$$
The sector certificates give the lower bound one-half. The [level $(-3,11)$ identity](../prime-density/README.md#a-genuine-half-density-obstruction) gives the matching upper bound there.

For each of the five new levels at $\ell=3$, the forms with scalar parameters $1$ and $1+\omega$ are in the same genus and their square-lift difference is
$$
 H_1(z)-H_{1+\omega}(z)=F(z)+c F(2z),\qquad
 F\in S_2(\Gamma_0(9|D|)),\qquad F\otimes\chi_{-3}=-F.
$$
The coefficient $c$ is $-1$ at $D=-7$ and $1$ at the other four levels. The replay reconstructs $F$ in the exact cusp space and checks the lift through the weight-two Sturm bound $6(|D|+1)$ at level $18|D|$. It also checks $F$ is anti-invariant under the Fourier twist through that bound. Twisting a trivial-character cusp form by $\chi_{-3}$ preserves level $9|D|$ by the [local conductor and depletion argument](../prime-density/README.md#a-genuine-half-density-obstruction). Thus the twist check is conclusive. At every good prime split in $\mathbb Q(\sqrt{-3})$, the difference has zero coefficient at every power of that prime. This rules out the entire split sector.

## Exact replay

The [coprime census](coprime.py) enumerates primitive norm-square vectors in the common ambient lattice. For a degree prime to $\ell$, an ambient primitive vector in $L_a$ is primitive there and conversely. The defining congruence becomes an affine line in the two trace coordinates of $a$ whenever the last two vector coordinates are nonzero as a pair modulo $\ell$. Accumulating these lines counts the contributions for all ambiguous classes together. Both vector signs are counted. The [sector replay](sectors.py) uses the exact divisor sums above and refines the list of possible ties separately in each sign sector. The [five half-density replays](half_density.py) use PARI only for the exact modular form and representation computations.

~~~sh
python -m pip install cypari2
study=research/elliptic-quotients/good-degree-spectrum
python -B "$study/coprime.py" --verify "$study/coprime.json"
python -B "$study/sectors.py" --verify "$study/sectors.json"
python -B "$study/half_density.py"
~~~

See the [coprime witnesses](coprime.json), the [sector certificates](sectors.json) and the [twist identities](half-density.json). The supplied comparison is commit `35d459c517b13fedfdf594962a7e292093000e43`.

Coprime-degree injectivity beyond the census remains open. At odd discriminants it is now sufficient for the uniform single-prime problem. The sector criterion also applies when a level might have more than one quadratic twist.

The [incidence transfer](../forced-square-kernel/README.md) computes the forced relations in the full split genus and identifies its residual kernel with the nonsplit one. It also gives uniform separation with residue-orbit labels at good degrees of order at most $\ell^2$ for fixed $D$.
