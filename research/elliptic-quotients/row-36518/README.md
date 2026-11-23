# All elliptic map degrees on row 36518

Every elliptic map from Sutherland row 36518 factors, up to translation on its target, through a single rational double cover to $E:Y^2=X^3+16$. Over an algebraic closure the possible degrees are all positive even integers. Over the rationals they are exactly $2n^2$ and $6n^2$ for $n\geq1$. Thus a degree-six elliptic map exists over the rationals on this curve. The same statement holds for an infinite family of geometrically distinct smooth quartics through the catalogue row. [Original data](https://math.mit.edu/~drew/genus3curves.html) · [Source equation](quartic.json) · [Exact checks](verify.py)

## The double cover

Set $(x,y,z)=(u-2v,u,v+t)$ in the original equation and put

$$
A=4u^2-8uv+2v^2,\qquad B=-2u^4+2u^3v-4u^2v^2+8uv^3-3v^4.
$$

The equation becomes $t^4+At^2+B=0$. More generally set $r=t^2+5kuv$ and $C_k:r^2+Ar+B=0$ for any integer $k$. Every member has the same smooth reduction modulo 5 and the rational point $[u:v:t]=[0:1:1]$. The fixed-point quotient of the involution $t\mapsto -t$ has coordinate $w=2r+A$ with $w^2=A^2-4B$. The identity

$$
A^2-4B=(2u)^3(u-v)+16(u-v)^4
$$

gives $X=2u/(u-v)$ and $Y=w/(u-v)^2$. This is a rational identification of the quotient with $E$ since $X$ is a fractional linear coordinate on the $[u:v]$-line. The nontrivial sign change witnesses the degree two of the map on smooth projective curves.

## The complementary surface

The counts at 5 for extension degrees $1,2,3$ are $(8,40,140)$ on the source and $(6,36,126)$ on the target. The verifier counts in each fibre over the finite projective line, including infinity. Newton identities give the characteristic polynomial $(T^2+5)h(T)$, with the elliptic factor $T^2+5$ and

$$
h(T)=T^4+2T^3+4T^2+10T+25.
$$

The exact root-ratio resultant is $625(U-1)^4q_1(U)^2q_2(U)$, with

$$
q_1=5U^4+6U^3+6U^2+6U+5,\qquad
q_2=25U^4+20U^3+26U^2+20U+25.
$$

All three quartics $h,q_1,q_2$ are primitive and irreducible modulo 11 with degree preserved. Neither $q_1$ nor $q_2$ is monic, so neither can have a root of unity. The four roots of $h$ have pairwise distinct powers at every exponent, and those powers form a single Galois orbit. The complementary abelian surface is therefore geometrically simple over the residue field. Any characteristic-zero elliptic factor would specialize after a finite extension at that good prime, contradicting simplicity.

Consequently any elliptic map kills the complementary surface at the level of Jacobians. Its difference from its composition with the involution is constant. Evaluation at any geometric fixed point of the sign change makes that constant zero. The map descends through the double cover, over its own field of definition. This yields exactly one geometric optimal elliptic quotient. Over the algebraic closure, cyclic isogenies of every positive degree exist from the fixed target, proving the geometric degree assertion. [Detailed specialization and factorization argument](../row-23640/README.md)

## Rational isogenies

The point $(0,4)$ has order three since its horizontal tangent meets the cubic with multiplicity three. A rational subgroup of order two would give a rational root of $X^3+16$, which is irreducible. Mazur’s prime-isogeny theorem leaves the odd primes $3,5,7,11,13,17,19,37,43,67,163$. For each candidate other than three the table gives a prime $p$ where the Frobenius discriminant $a_p^2-4p$ is a nonsquare modulo that candidate. This excludes a rational line in the corresponding torsion. [Mazur, Theorem 1](https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf)

| $p$ | $a_p$ | Newly excluded prime degrees |
|---:|---:|---|
| 5 | 0 | 11, 13, 17, 19, 37 |
| 7 | −1 | 5 |
| 11 | 0 | 7, 43 |
| 17 | 0 | 67 |
| 41 | 0 | 163 |

For a Galois-stable kernel with elementary divisors $a\mid b$, its subgroup killed by $a$ is $E[a]$. The image under multiplication by $a$ is cyclic of order $b/a$. Every prime dividing $b/a$ must yield a rational prime isogeny from $E$. Thus every rational isogeny degree is $a^2 3^j$, either a square or three times a square. Multiplication maps and their compositions with the rational 3-isogeny realize both sets. The source double cover multiplies each degree by two.

For the family, $E$ is fixed and the branch divisors of the double covers are the zeros of $w-A-10kuv$. These give distinct sections up to scalar of the same degree-four line bundle on $E$. Any isomorphism of the source curves must identify their unique degree-two quotient maps up to an automorphism of $E$. The stabilizer of the line bundle has at most $16\cdot6=96$ geometric elements, from $E[4]$ and the six pointed automorphisms at $j=0$. Each geometric isomorphism class can therefore occur for at most 96 integer parameters.

## Replay

With Python 3.11 or later and SymPy installed, run `python verify.py` in this folder. The program reuses the small finite-field arithmetic from the adjacent row 23640 package. It verifies the original-coordinate identity, smooth reduction, point counts, factorization, exact root-ratio certificate and isogeny exclusions. [Recorded result](receipt.json)
