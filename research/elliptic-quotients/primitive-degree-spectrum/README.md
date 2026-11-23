# Counting every primitive elliptic map on an infinite family of quartics

For the [smooth plane quartics $C_a$](../primitive-degree-four/README.md) with integer $a\equiv4\pmod {15}$, the primitive elliptic maps over an algebraic closure can be listed in every degree, and every class has a representative over the rationals. A map is primitive here if there is no intermediate curve with both map degrees greater than one. Classes are taken up to target isomorphism and translation; the source is fixed.

**Theorem.** Write $P,Q$ for the two maps to $E_a:y^2=x(x-1)(x-a)$ in the construction. There is one further double cover $R$ to a second elliptic curve. The full list of primitive classes is $R$ together with
$$
[m]P+[n]Q,\qquad (m,n)\in\mathbb Z^2,\quad \gcd(m,n)=1,
$$
where simultaneous sign gives the same class. The degree of the displayed map is $2(m^2+n^2)$. Each single source in this progression therefore has primitive maps of unbounded degree, with explicit recipes over $\mathbb Q$. For example, there are three primitive classes of degree two, two of degree four, none of degree six, and four of degree 26. The arithmetic progression contains infinitely many geometrically distinct sources because the invariant of their degree-two quotient varies and a fixed smooth genus-three curve has finitely many involutions.

## The three elliptic factors

Over the $r$-line write $y_1^2=r(r-1)(r-a)$ and $y_2^2=(r+1)r(r+1-a)$. Then $P=(r,y_1)$ and $Q=(r+1,y_2)$. The invariant $z=y_1y_2/r$ gives the third sign quotient
$$
R:C_a\longrightarrow E'_a,\qquad z^2=(r-1)(r-a)(r+1)(r+1-a).
$$
Its two points over infinity are rational. The three maps have different characters under the two sign involutions. Hence each cross composition of a pushforward and a different pullback on the Jacobians is zero, while composition along the same map is multiplication by two. The product of the three elliptic pullbacks therefore gives an isogeny $E_a^2\times E'_a\to J(C_a)$, by independence and equality of dimensions.

Put $j(\lambda)=256(1-\lambda+\lambda^2)^3/[\lambda^2(1-\lambda)^2]$. The two invariants are $j(a)$ and $j(\lambda')$ with $\lambda'=a(a-1)/[(a-2)(a+1)]$, a cross ratio of the four roots in $E'_a$. At $a\equiv4\pmod {15}$ the first invariant has negative 3-adic valuation, so $E_a$ has no complex multiplication. At 5 the first invariant is integral. The cross ratio $\lambda'$ has negative 5-adic valuation, and hence $v_5(j(\lambda'))=2v_5(\lambda')<0$. The two curves cannot be geometrically isogenous. A cyclic prime isogeny between them would carry an integral invariant to an integral invariant, by the monic integral modular polynomial; any isogeny over the algebraic closure can be factored into prime isogenies. The second curve is also non-CM. At $a=4$ the two invariants are $35152/9$ and $1906624/225$. [Classical modular polynomials](https://arxiv.org/abs/1202.3985)

## Which directions give primitive maps?

Put $f_{m,n}=[m]P+[n]Q$. Orthogonality gives $f_*f^*=[2(m^2+n^2)]$. Also $P_*f^*=[2m]$ and $Q_*f^*=[2n]$. For coprime $m,n$ this forces $\ker f^*\subseteq E_a[2]$. When both are odd, the restriction is $(P+Q)^*$; the three squareclass tests in the [degree-four proof](../primitive-degree-four/README.md) exclude all three nonzero elements. When just one integer is odd, the restriction is $P^*$ or $Q^*$. The function for an individual 2-torsion class then has exactly two odd places on the $r$-line, one of them infinity. Each of the three quadratic subfields of the degree-four source has four odd places. This excludes all three classes for the individual covers as well. Thus $f^*$ has trivial kernel and the dual quotient $f_*$ has connected kernel. The same holds for the third double cover $R$, since a factorization through an elliptic isogeny would leave a degree-one map from a genus-three curve to a genus-one curve.

A higher-genus intermediate would have to have genus two and receive an unramified map of degree two from the source, by Riemann–Hurwitz. A deck transformation of such a map would be fixed-point-free. Every involution on a smooth plane quartic acts on the canonical plane and fixes a line meeting the curve. Genus zero cannot map nontrivially to an elliptic curve. This establishes full geometric primitivity.

Conversely, a map to any elliptic curve induces a nonzero homomorphism from $J(C_a)$, so its target is geometrically isogenous to $E_a$ or $E'_a$. In the first case, compose with an isogeny to $E_a$. The rational Hom space is two-dimensional with basis $P_*,Q_*$, since $\operatorname{End}^0(E_a)=\mathbb Q$ and the other elliptic factor is nonisogenous. Clearing denominators gives a unique rational direction represented by a coprime pair up to sign. The identity component of the kernel is unchanged by postcomposition with an isogeny or by multiplication by a nonzero integer. The optimal quotient in this direction is exactly $f_{m,n}$, whose kernel on Jacobians is connected. Every other map in the direction factors through it up to translation. The second Hom space gives the single direction $R$. All exhibited maps are rational, and a rational point at infinity on the source allows normalization by translation over $\mathbb Q$.

For $d=2N$, the number of primitive classes is
$$
\mathbf 1_{N=1}+\tfrac12\#\{(m,n):\gcd(m,n)=1,\ m^2+n^2=N\}.
$$
It is zero for odd $d$. In the Gaussian integers, a primitive norm is divisible by no prime $3\pmod4$ and by at most one factor of 2. Each split prime dividing $N$ contributes two choices, since all its powers have to be placed in one of its two Gaussian factors. The four units and identification by simultaneous sign leave $2^{1+s}$ directions, where $s$ is the number of distinct split primes, whenever those divisibility conditions hold.

Every geometric map factors through its unique primitive quotient. A primitive elliptic quotient has exactly $\sigma_1(h)=\sum_{e\mid h}e$ target kernels of order $h$. Thus if $c(d)$ is the primitive class count, the full geometric count in degree $D$ is $\sum_{d\mid D}c(d)\sigma_1(D/d)$. For example, there are eleven geometric degree-four classes. Nine come from the three double covers followed by one of three 2-isogenies, and two are primitive.

## All rational map degrees in a smaller progression

For every $a\equiv4\pmod {4{,}849{,}845}$, including $a=4$, the degrees of **all** rational elliptic maps are exactly
$$
\{2(m^2+n^2):(m,n)\in\mathbb Z^2\setminus\{(0,0)\}\}.
$$
Over the algebraic closure every positive even degree occurs via isogenies after a double cover. Every geometric primitive degree is even, so no odd degree occurs.

The modulus is $15\cdot7\cdot11\cdot13\cdot17\cdot19$. At those last five good primes, the two targets have the following Frobenius traces throughout the progression.

| Prime | $E_a$ | $E'_a$ |
|---:|---:|---:|
| 7 | 0 | 0 |
| 11 | −4 | 0 |
| 13 | −2 | −2 |
| 17 | 2 | 6 |
| 19 | 4 | −4 |

For each target, these traces give a nonsquare $a_p^2-4p$ modulo every allowed odd prime isogeny degree in Mazur’s theorem, using a different residue characteristic from that candidate. Both targets have rational 2-torsion. For the quartic this follows from the rational branch roots. A rational isogeny kernel with elementary divisors $b\mid c$ contains the full $b$-torsion; its cyclic remainder can only have 2-power order. The rational isogeny degrees are therefore squares or twice squares. Postcomposing a primitive map by these degrees preserves the form twice a sum of two squares, by the two-square product identity. Conversely, the displayed maps work for all integer pairs. [Mazur, Theorem 1](https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf) · [Exact residue witnesses](receipt.json)

## Replay and use

```sh
python verify.py
python degrees.py 4 26 --directions
```

The second command returns $(2,3),(2,-3),(3,2),(3,-2)$. Each pair is an explicit group-law recipe using $P$ and $Q$ above. The direction list is optional. The count alone uses integer factorization. Outside the smaller progression the program returns a rational verdict only when the results above decide it. The [verifier](verify.py) checks the branch parities, independent residue counts, odd-isogeny exclusions, and degree counts through 500 against direct lattice enumeration.
