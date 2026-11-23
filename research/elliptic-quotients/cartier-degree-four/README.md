# Recovering every optimal elliptic quotient on a fresh positive quartic

The [quartic equation](quartic.json) alone yields two geometrically primitive degree-four maps and one double quotient. These are all its optimal elliptic quotient classes over the algebraic closure, up to translation and target isomorphism. A second input mixes the projective coordinates and has all fifteen monomials. The calculation recovers the same three classes from that equation independently. Each run finds its source point, differential directions and target coefficients. The final certificate bounds all geometric directions regardless of coefficient height or field of definition. [Mixed equation](mixed.json) · [Full rational map formulas and finite-field receipts](receipt.json)

The source was constructed as a positive benchmark. In its first set of coordinates its equation is

$$
29880U^4-42336U^3V+21360U^2V^2-4416UV^3+320V^4
+(-489U^2+348UV-52V^2)W^2+2W^4=0.
$$

| Canonical differential line | Optimal degree | Target invariant $j$ |
|---|---:|---:|
| $W$ | 2 | $1728$ |
| $U$ | 4 | $17779581/32166277$ |
| $2U-V$ | 4 | $11664/13$ |

The unique double quotient distinguishes this curve from the previous three-involution family and the five-involution benchmark at $a=3$. For an integer $d>0$, all geometric map classes of degree $d$, including isogeny compositions, number
$$
\sigma_1(d/2)+2\sigma_1(d/4),
$$
where $\sigma_1(n)$ is the sum of the positive divisors of $n$ and is zero at a noninteger argument. Indeed every elliptic map factors through its unique connected-kernel quotient followed by an isogeny. An elliptic curve over the algebraic closure in characteristic zero has exactly $\sigma_1(n)$ finite subgroups of order $n$. In particular the quartic has five degree-four classes in total and admits precisely the positive even map degrees.

## How the equation determines the maps

The Cartier operator preserves the pulled-back differential line of every elliptic map over a good prime where the pullback is nonzero, by functoriality of Cartier. In a plane quartic this operator on the three canonical coordinate differentials can be computed from the coefficients of $F^{p-1}$. Put $e_1=(2,1,1),e_2=(1,2,1),e_3=(1,1,2)$. With images as columns its matrix over $\mathbb F_p$ has entries
$$
C_{ij}=[X^{p e_i-e_j}]F^{p-1}.
$$
The algorithm enumerates primitive integral coordinate triples up to simultaneous sign and retains the common projective eigendirections modulo its good primes. In this example height four suffices in both coordinate presentations. This finite search supplies candidates. Completeness comes from the independent specialization argument below. [Costa, Harvey and Sutherland, *Counting points on smooth plane quartics*](https://doi.org/10.1007/s40993-022-00397-8)

Given a candidate line $L$, the four coefficients of the restriction of $D_PF$ to $L=0$ are linear in the possible involution center $P$. Their matrix has rank at least two. Otherwise its three binary cubic gradient restrictions would be multiples of one binary cubic, producing a common zero and hence a singular point on the line. Any candidate center is consequently unique. The program tests $L(P)\ne0$ and the reflection identity on the full quartic in a frame for $L$ and $P$. A successful identity gives a double quotient. For a degree-four map the same test identifies precisely the imprimitive cases, by the fixed-line criterion. [Involution proof](../involution-atlas/README.md) · [Degree-four criterion](../degree-four-jets/README.md)

For the remaining lines select a rational point where the line is nonzero and a smooth affine chart $F(t,s)=0$ with $F_s(0)\ne0$. Write the normalized differential as $h(t)dt$ with $h(0)=1$, and send the selected point to the elliptic origin by translation. It remains to find the two coefficients $A,B$ of $y^2=x^3+Ax+B$. At each prime above 27 the program enumerates **every** nonsingular normalized pair in $\mathbb F_p^2$. It computes the first 27 coefficients of $R=t^2x$ from
$$
(tR'-2R)^2=4h^2(R^3+At^4R+Bt^6),\qquad R(0)=1.
$$
Only the scalars $4(n+1)$ with $1\le n\le26$ must be inverted in this recurrence. Primes that divide chart denominators or lose a required fixed matrix rank are skipped.

A degree-four map gives at least two linearly independent cubic numerator-denominator pairs $N,D$ satisfying $t^2N-RD=0$. The search eliminates the ten numerator coefficients by ordinary linear algebra. It also requires $D$ to vanish through first order, leaving eight denominator coefficients. The resulting truncated matrix is $15\times8$, with entries linear in the 27 coefficients of $R$. Rank at most six is necessary. A batch test of one $7\times7$ minor rejects most targets, and full finite-field elimination checks every surviving candidate. When one smooth pair survives at a prime, its two residues enter a Chinese remainder computation. Two consecutive equal rational reconstructions trigger a characteristic-zero check. A modular coincidence alone never certifies a map.

In the original coordinates each direction used 16 primes from 29 through 97, covering 64,240 parameter pairs per direction including the singular pairs discarded directly. In the mixed coordinates the run used 22 primes through 131, covering 149,398 pairs per direction over the accepted primes. On the original input the point is $[2:11:16]$ with $U=2+t,V=11+s,W=16$. The two characteristic-zero target pairs are
$$
(464/10673289,5806720/941480149401),\qquad
(38416/10673289,-75295360/941480149401).
$$
The 57-term verifier gives rank 18 for each, reconstructs $x$ as a quotient of cubics, and verifies the global elliptic identity by exact polynomial division. The other coordinate is $y=dx/(2h(t)dt)$, interpreted using the rational differential on the entire quartic. Both lines pass the primitivity check. Under $(U,V,W)=(X+Y,Y+Z,X+Y+Z)$ the search independently finds $[5:-3:14]$ and the three predicted projective lines. [Finite-field search](recover.py) · [End-to-end program](cartier_search.py)

## Geometric exhaustion at 17

The verified reflection yields an even model $aW^4+b(U,V)W^2+c(U,V)=0$. The code checks invertibility of the projective change modulo 17. In every fibre over $[U:V]$ it solves the quadratic in $W^2$ and sums its actual square-root counts, including multiple roots and the point at infinity on the base. The counts over the first three extensions are $30,312,4878$. The homogeneous Jacobian ideal certifies smoothness there. The real Weil polynomial is $s(s+4)(s+8)$, where $s=T+17/T$. The trace-zero elliptic block is supersingular. The other two blocks are ordinary with distinct geometric endomorphism fields $\mathbb Q(\sqrt{-13})$ and $\mathbb Q(i)$. Equivalently the replayed root-ratio calculation separates all three elliptic blocks over the algebraic closure of the residue field and gives geometric multiplicity one in each.

An optimal elliptic quotient determines a Rosati-symmetric rank-two projector on the Jacobian. Distinct optimal directions give different projectors. Specialization at good reduction is injective on rational endomorphisms after a field extension defining them. A multiplicity-one geometric elliptic block has just one such projector. Therefore at most three optimal directions can occur, over any number field and in any degree. The three verified directions attain this upper bound. A smooth plane quartic in characteristic zero has no unramified involution, so the connected-kernel elliptic quotients are also primitive as curves. [Specialization and projector proof](../polarized-reduction/README.md)

## Provenance

One construction starts with $D:y^2=x^3+3x$ and the function $g=x^2+x+3-2y$. At $k=2$ the normalization of the double cover $t^2=g(16g^2-129g+363)$ supplies the benchmark. Its canonical coordinates can be taken as $U=-3x^2/2+xy+y/2-3/2$, $V=x^3-3x^2/4+3x/2+5y/4+3/4$, and $W=t$. The discriminant of the quartic for $x$ over $\mathbb Q(g)$ is $-256g^2(g-1)(16g^2-129g+363)$. The two nonzero finite branch values of the elliptic double cover each have fibre type $(2,1,1)$ on $D$, while $0$ has type $(2,2)$ and infinity has type $(4)$. This explains how to choose a positive example. The discovery program reads only the submitted plane equation.

## An infinite family with all optimal degrees determined

For every integer $k\equiv2\pmod {17}$ put
$$
B=k^2-k+1,\quad D_k:y^2=x^3+(k-2)x^2+Bx,\quad
g=x^2+(k-1)x+B-2y,
$$
and let $C_k$ be the smooth projective normalization of
$$
t^2=gH(g),\qquad H(z)=16z^2-(24k^2-4k+41)z+9k^4+21k^3+4k^2+5k+25
$$
over $D_k$. Every $C_k$ is a nonhyperelliptic genus-three curve. Its three optimal elliptic quotients have degrees **2, 4, 4** over $\overline{\mathbb Q}$ and are defined over $\mathbb Q$ as Jacobian quotients. The full geometric map count $\sigma_1(d/2)+2\sigma_1(d/4)$ therefore holds throughout this progression.

To check the family at 17, the quartic defining $x$ over the $g$-line is
$$Q=(x^2+(k-1)x+B-g)^2-4(x^3+(k-2)x^2+Bx).$$
Its discriminant in $x$ is $-256g^2(g-(k-1)^2)H(g)$, and $Q(x,0)=(x^2+(k-3)x+B)^2$. Modulo 17, $D_k$ is smooth and the latter quadratic has distinct roots. The two roots of $H$ are simple, nonzero, and different from $(k-1)^2$. Each corresponds to one simple ramification and two unramified points of the degree-four function $g$. Infinity is its unique pole divisor $4O$. These statements hold on both fibres over the local ring at 17, since the corresponding discriminants and resultants are units there.

More explicitly the divisor of $gH(g)$ is $2A+R-12O$, with $A$ the two points above $g=0$ and the double point above each root of $H$. The four remaining unramified points above the roots of $H$ form $R$. These disjoint sets extend as étale divisors after an unramified base extension that labels the points. This follows from the squarefree factors at $g=0$ and Hensel lifting at each simple discriminant zero. The double cover defined by the square root, using the line bundle $6O-A$, consequently has smooth reduction, ramified exactly at the four points of $R$. Riemann–Hurwitz gives genus three. Its special fibre is the worked quartic. Indeed the three displayed functions $U,V,W$ specialize to its equation identically. Their ratio $U/V$ is nonconstant with pole degree at most six on $D_2$ and at most twelve on the cover, so the resulting map to the smooth plane quartic is separable even at 17. Two genus-three curves admit no separable map of degree greater than one. The generic fibre is nonhyperelliptic as well. The multiplication map from $\operatorname{Sym}^2 H^0(K)$ to $H^0(2K)$ has rank six on the special plane quartic, hence also on the generic fibre, by base change for the canonical powers.

The degree-two map to $D_k$ and degree-four map to $E_k:t^2=zH(z)$ exist over $\mathbb Q$. The second map has degree four since the source double extension of $D_k$ is ramified at $R$. At $k=2$, the invariant of $E_2$ identifies it with the first primitive degree-four direction up to geometric isomorphism because both targets have invariant $17779581/32166277$, and the other two elliptic blocks are geometrically nonisogenous. The three special optimal degrees are exactly $2,4,4$. Indeed the prime-to-17 part of each special degree agrees with its corresponding characteristic-zero degree at $k=2$, and each local projector multiplier is a power of two at 17. For all $k$ in the progression $D_k$ and $E_k$ specialize to these same distinct blocks. The optimal degree underlying the exhibited map to $E_k$ divides four and has prime-to-17 part four, so it is four.

For the remaining elliptic subvariety put $e_F=1-e_D-e_E$, where $e_D,e_E$ are the symmetric projectors onto the first two distinct blocks. The two projectors are orthogonal since the blocks are geometrically nonisogenous. Their complement has rank two on first homology and defines an elliptic subvariety over $\mathbb Q$. The endomorphism $4e_F=4-4e_D-4e_E$ is integral because the first two polarizing degrees are 2 and 4. The complementary optimal degree divides four by the projector denominator criterion. Its prime-to-17 part is four on the remaining special block. It follows that the degree is exactly four. The common local bound of three optimal directions proves exhaustion. The line above infinity has two rational lifts on $C_k$ since the leading coefficient of $gH(g)$ at its pole is 16, allowing the maps themselves to be based over $\mathbb Q$.

These curves include infinitely many geometric isomorphism classes. The $j$-invariant of $E_k$ is a nonconstant rational function of $k$, as is already seen from its different values at 2 and 19. Each fixed source in this family has only three optimal target classes. [Symbolic identities and specialization checks](verify.py)

## Replay

From this folder, with NumPy and SymPy installed, run

```sh
python cartier_search.py quartic.json
python cartier_search.py mixed.json
python verify.py
```

The numerical limits control the candidate search. The exact output field `all_optimal_directions_exhausted` requires the certified count of distinct maps and involutions to attain the independent local bound. A smaller height, insufficient good primes, nonrational directions, missing rational chart points or repeated geometric blocks can leave that field false on other inputs. The general characteristic-zero elimination of all four jet parameters is a separate method. The successful benchmark here uses arithmetic reconstruction and a geometric upper bound. [Verification of both presentations and a cutoff control](verify.py)
