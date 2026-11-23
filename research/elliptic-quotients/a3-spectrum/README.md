# All elliptic maps on the primitive degree-four example

For the [plane quartic](../primitive-degree-four/ciani.json) with parameter $a=3$, the complete set of degrees of nonconstant maps to varying elliptic curves over $\mathbb Q$ is

$$
\{2(m^2+n^2):m,n\in\mathbb Z,\ (m,n)\ne(0,0)\}.
$$

The degrees of rational maps that admit no nontrivial intermediate curve even over $\overline{\mathbb Q}$ are given by the same formula with $\gcd(m,n)=1$. Examples of primitive degrees are $2,4,10,20,26,34,50$. Degree six is absent over the rationals. These results concern the one fixed curve at $a=3$. The geometric classification below uses three integer coordinates. In particular there are four geometric degree-three classes and all four require the quadratic field $\mathbb Q(\sqrt{-2})$. The construction of primitive maps $mP+nQ$ itself works throughout the integer family. [Quartic and maps](../primitive-degree-four/README.md) · [Arithmetic checks](verify.py) · [Output](receipt.json)

## Primitive maps of unbounded degree on one source

Let $P,Q:C\to E$ be the two rational double covers constructed in the main study. The corresponding homomorphisms on Jacobians satisfy $P_*P^*=Q_*Q^*=[2]$ and $P_*Q^*=Q_*P^*=0$. For coprime integers $m,n$, put $f=mP+nQ$. Then $f_*f^*=[2(m^2+n^2)]$. Also $P_*f^*=[2m]$ and $Q_*f^*=[2n]$, so every element of $\ker f^*$ is killed by two.

When $m,n$ are both odd, $f^*$ on $E[2]$ equals $(P+Q)^*$, which the three branch squareclass tests prove is injective. When exactly one is odd it equals $P^*$ or $Q^*$. Each of these is injective: a nonzero element in the kernel would define an étale double isogeny through which the ramified degree-two cover factors. A degree-one map from the source to that elliptic cover is impossible. Thus $f^*$ is injective and $f_*$ has connected kernel, by duality of abelian varieties. Up to translation this quotient cannot factor through a nontrivial isogeny.

An intermediate curve of genus two would require a degree-two map from $C$ to that curve by Riemann–Hurwitz. Such a map would be unramified. Every involution of the smooth plane quartic acts on its canonical plane and fixes a line meeting $C$, excluding an unramified double quotient. Riemann–Hurwitz excludes higher genus intermediates with either nontrivial incoming degree or genus exceeding that of $C$. A rational curve cannot map nontrivially to $E$. This proves geometric primitivity of every displayed coprime map.

## The full lattice over the rationals

Write $E:y^2=r(r-1)(r-3)$. The third sign quotient $R$ maps to the curve $E'$ with $y^2=g(r)=(r-1)(r-3)(r+1)(r-2)$. It has a rational point at infinity. The identity

$$
X(X-1)(X-3)(r-1)^4=-2g(r),\qquad X=2(r-2)/(r-1)
$$

identifies $E'$ with a quadratic twist of $E$ by $-2$. At 5 the curves have respective traces $2$ and $-2$, so they are nonisogenous over the rationals. The invariant of either is $21952/9$. A complex multiplication invariant is an algebraic integer, so both geometric endomorphism rings are $\mathbb Z$.

The three elliptic images $P^*E,Q^*E,R^*E'$ in the Jacobian have pairwise zero cross compositions, by their different characters under the two sign changes. For the $R$ quotient the same computation gives $R_*R^*=[2]$. Their addition is an isogeny to the three-dimensional Jacobian, since the three pushforwards separate the images. Comparison on each image gives

$$
P^*P_*+Q^*Q_*+R^*R_*=[2].
$$

Normalize maps to respect a rational base point, or work with their induced homomorphisms. For any rational homomorphism $h:J_C\to E$, the identity gives $2h=aP_*+bQ_*$ with integers $a,b$. The $R$ term vanishes by rational nonisogeny. Subtract integral multiples of $P_*,Q_*$ from $h$ so the two coefficients on the right are zero or one. If exactly one is one, the right side has map degree two while the left degree is divisible by four. If both are one, the right side has degree four, requiring a degree-one map from $C$ to $E$. Thus both are zero and $\operatorname{Hom}_{\mathbb Q}(J_C,E)=\mathbb ZP_*\oplus\mathbb ZQ_*$. The same reasoning gives $\operatorname{Hom}_{\mathbb Q}(J_C,E')=\mathbb ZR_*$.

Every rational elliptic target admitting a map is isogenous over the rationals to $E$ or $E'$. Compose the map with an isogeny back to that representative. The result is a multiple of a primitive direction $mP+nQ$ or of $R$. Composition with the isogeny preserves the identity component of the kernel of the induced Jacobian map. Since the primitive direction has connected kernel, the original map factors through its quotient, over the rationals and up to translation. The curve has rational points at infinity, so the normalization of maps by translation is defined over this base field.

## The remaining isogenies

Both elliptic targets have full rational two-torsion. The explicit prime-isogeny sieve excludes every other possible prime degree from Mazur’s list. The squared traces of the two quadratic twists coincide. At primes $5,7,11,17$ their discriminants supply a nonsquare modulo every candidate odd prime. The proof uses the fact that a Galois-stable cyclic subgroup of prime order would give a Frobenius-stable line in the torsion representation. [Mazur, Theorem 1](https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf) · [Exact witnesses](receipt.json)

A rational isogeny kernel with elementary divisors $a\mid b$ contains all $a$-torsion. After multiplying by $a$, the residual cyclic kernel of order $b/a$ is Galois-stable and each of its prime factors must be two. Thus the isogeny degree is a square or twice a square. Each square and twice-square occurs. Multiplication scales an integral sum-of-two-squares vector. The identity $2(m^2+n^2)=(m+n)^2+(m-n)^2$ absorbs an additional factor two as well. It follows that every rational elliptic map has one of the degrees stated, and $mP+nQ$ realizes each value. For primitive maps, the intermediate isogeny must have degree one.

## Three directions over the algebraic closure

Fix an isomorphism $E'\to E$ over $K=\mathbb Q(\sqrt{-2})$ using the fractional linear $X$ above and $Y=\sqrt{-2}\, y/(r-1)^2$. Let $S$ denote its composition with the third quotient map $R$, working up to translation. Its geometric endomorphism coordinates are now in $\mathbb Z$ along with those of $P,Q$. For any $h:J_C\to E$, the same projector identity gives $2h=aP_*+bQ_*+cS_*$ with integers $a,b,c$. Subtract integer multiples to reduce the coefficients to zero or one. A residue vector of weight $j$ has degree $2j$ by orthogonality. For weight one or three this degree is indivisible by four. Weight two would require a degree-one map from $C$ to $E$. Thus all coefficients are even, proving $\operatorname{Hom}_{\overline{\mathbb Q}}(J_C,E)=\mathbb Z^3$ in this basis.

For a primitive integer vector $v=(a,b,c)$, put $h_v=aP+bQ+cS$. The three adjoint compositions imply $\ker h_v^*\subseteq E[2]$, since $\gcd(a,b,c)=1$. Its exact order is determined by the three functions $x-e$, $e=0,1,3$, on $E$. Their pullbacks for the sum of the maps represent the products of the individual pullback squareclasses. For $P$ the odd places of $r-e$ are $e,\infty$. For $Q$ they are $e-1,\infty$. For $S$ they are $X^{-1}(e),1$, where $X^{-1}(0)=2$, $X^{-1}(1)=3$, $X^{-1}(3)=-1$. The sign change or translation by the chosen base point has no effect on the induced homomorphisms on $\operatorname{Pic}^0$.

On the base line, the odd support of a product is the symmetric difference of these supports, taken with multiplicities $a,b,c$ modulo two. It is a square in the function field of $C$ precisely when the result is one of $\varnothing,B_1,B_2,B_1\mathbin{\triangle}B_2$, where $B_1=\{0,1,3,\infty\}$ and $B_2=\{-1,0,2,\infty\}$. This is the squareclass criterion for the three quadratic subfields of the biquadratic cover. The verifier independently factors the actual rational functions.

| $(a,b,c)\bmod 2$ | Nonzero point of $E[2]$ killed by $h_v^*$ |
|---|---|
| 001, 010, 100, 101, 110 | None |
| 011 | $(1,0)$ |
| 111 | $(3,0)$ |

Consequently the primitive geometric quotient attached to $v$ has degree
$$
d(v)=\frac{2(a^2+b^2+c^2)}{\delta(v)},\qquad
\delta(v)=\begin{cases}2& b,c\text{ both odd},\\1&\text{otherwise}.\end{cases}
$$
Indeed, quotienting $J_C$ by the connected kernel of the induced map gives an elliptic curve $E_v$. The remaining isogeny $E_v\to E$ is dual to the isogeny with kernel $\ker h_v^*$, so its degree is $\delta(v)$. The map to $E_v$ cannot have an elliptic intermediate or a higher-genus intermediate by the genus argument above.

These exhaust the primitive geometric classes, taken up to target isomorphism and translation. Any elliptic target receiving a map is geometrically isogenous to $E$. Composing with an isogeny to $E$ produces a nonzero lattice vector. Its rational line specifies the same connected kernel. Conversely the same connected kernel gives the same rational line because $\operatorname{Hom}(E_v,E)$ has rank one for a non-CM curve. Thus the classes are in bijection with primitive integer triples modulo simultaneous sign. Conjugation over $K/\mathbb Q$ acts on these lines by $(a,b,c)\mapsto(a,b,-c)$. Exactly the lines with $c=0$ or $a=b=0$ are rational, where rational optimal maps were exhibited above. The remaining lines have exact class field $K$.

For any positive degree $D$, every geometric class is obtained by composing its primitive quotient with an isogeny. An elliptic curve over an algebraically closed field of characteristic zero has $\sigma_1(N)=\sum_{e\mid N}e$ subgroups of order $N$, by the index-$N$ sublattice count. The full number of classes of degree $D$ is
$$
\sum_{v\,\mathrm{primitive}/\{\pm1\},\ d(v)\mid D}\sigma_1\big(D/d(v)\big).
$$
The sum is finite because $a^2+b^2+c^2\le d(v)\le D$. The [counting program](map_classes.py) evaluates it for any supplied positive degree.

| Degree | Primitive geometric classes | Primitive classes defined over $\mathbb Q$ | All geometric classes including isogeny compositions |
|---:|---:|---:|---:|
| 2 | 5 | 3 | 5 |
| 3 | 4 | 0 | 4 |
| 4 | 4 | 2 | 19 |

In degree four, fifteen classes come from following each of five double covers by one of three 2-isogenies. The other four are primitive. The degree-three classes are exactly the four sign lines of $(\pm1,\pm1,\pm1)$.

All even degrees occur because the curve has a double cover and geometric isogenies exist in every degree. An odd primitive degree requires $a,b,c$ all odd and equals $a^2+b^2+c^2\equiv3\pmod 8$. Conversely every positive $d\equiv3\pmod 8$ is a sum of three odd squares by the Gauss–Legendre three-square theorem. Divide the representing triple by its common divisor $g$ to get a primitive triple of degree $d/g^2$. A subsequent isogeny of degree $g^2$ realizes $d$. The complete geometric degree set therefore comprises all positive even integers and the odd multiples of positive integers congruent to three modulo eight. [Three-square theorem and algorithmic proof](https://doi.org/10.1090/mcom/3349)

## Independent cubic-count upper bound

The original Part XXX geometric solver at 11 selects the residue field of size 121. It finds eight special classes. The upper bounds at precisions $11$, $11^2$, $11^3$ are respectively $8,8,4$. The bounded solver itself reports a lower bound of zero, since it recovers no characteristic-zero formulas in this run. The four distinct sign directions constructed above supply the matching lower bound of four. This verifies the geometric cubic count by a second route and shows explicitly why the first residue-field scan overcounts it. [Original engine receipt](cubic-lifting.json) · [Replay](verify_cubic_lifting.py)

## Replay

Run `python verify.py` with SymPy installed to check the twist identity, the nonisogeny witness, the prime-isogeny sieve, all seven parity rows and the small-degree counts. Run `python map_classes.py 4 --directions` for a list of all contributing primitive directions with their fields and isogeny multiplicities.
