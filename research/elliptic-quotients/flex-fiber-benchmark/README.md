# A plane quartic with exactly one primitive elliptic quotient, of degree four

The [input](quartic.json) is a smooth plane quartic over the rationals in mixed coordinates. It has no geometric involutions. From its coefficients alone, the [solver](branch_fibres.py) finds a degree-four map, identifies its elliptic target, and proves that every elliptic map over the algebraic closure factors through this one. Thus its geometric degrees are exactly $4n$ for positive integers $n$. Its rational degrees are exactly $4n^2$ and $12n^2$. There is a single rational class in each of those degrees, up to target isomorphism and translation. The same conclusions hold on an infinite arithmetic family described below.

## Running both positive benchmarks

```sh
python branch_fibres.py quartic.json
python verify.py
python benchmark.py --seconds 90
```

The fast run takes the quartic and a good prime (5 by default). No projection, branch points, elliptic invariant, differential or map formulas are supplied. It solves a two-variable polynomial system exactly and verifies the map by substitution in the source function field. A geometric upper bound from the special fibre certifies exhaustion of all elliptic directions regardless of the search strategy. If the certificates do not meet, the program returns an unresolved status. The reconstruction searches rational projection centres on the given line at infinity and a particular trigonal branching pattern. That search alone cannot exclude maps on other inputs. [Exact output](receipt.json)

The timed benchmark also runs the earlier unrestricted incidence on the chart containing the positive differential. In mixed coordinates its line is $36Z-24X-24Y=0$, so the last nonzero coefficient is the third. The [timing record](benchmark.json) gives the runtime, solver versions and outcome of both commands. The wall-clock cap applies only to the raw elimination. Linux also limits that subprocess to one GiB of address space. The [general elimination theorem](../degree-four-incidence/README.md) applies on all three charts; this experiment times the relevant positive chart.

## Recovering the map

Write $X,Y,Z$ for the supplied coordinates. On $Z=1$ the algorithm finds a rational point at infinity and projects from it. The resulting coordinates are $t=X+Y$ and $s=Y$. The source equation is cubic in $s$. Its discriminant as a polynomial in $t$ has two double roots $a=1/4, b=1$, one simple quadratic factor and one simple quartic factor. The program uses the rational factors and tries each degree-four product $H(t)$ of the simple factors. In this example there is exactly one possibility.

Consider a degree-four function $g$ on the $t$-line with a simple zero at $a$ and a simple pole at $b$. When its remaining zero and pole have order three, $g=A/B_0$ with $A=c(t-a)L(t)^3$ and $B_0=d(t-b)M(t)^3$, where $L,M$ are affine linear and may be constant. Take the four roots of $H$ to be the fibre $g=1$. Normalize $L(b)=M(a)=1$. The equation $A-B_0=H$ fixes $c=H(b)/(b-a)$ and $d=H(a)/(b-a)$. Writing $L=1+l(t-b)$ and $M=1+m(t-a)$, the difference $A-B_0-H$ is divisible by $(t-a)(t-b)$. The coefficients of the remaining quadratic give three equations in $l,m$. The complete lexicographic basis on this input is

$$l-\tfrac12,\qquad m.$$

The two remaining simple branch points have a common value $\beta=1/9$ under $g$. Seek a target cubic in $x,z$ with $z=g$,
$$x^3-Dzx-Dz(z+B)=0.$$
Its discriminant as a cubic in $x$ has the other simple roots $1,\beta$ when $B^2=\beta$ and $D=27(1+\beta+2B)/4$. The solver tries both rational square roots. Depress the source cubic to $v^3+p(t)v+q(t)=0$. An identification $x=kv$ must satisfy $k^2p=-Dg$ and $k^3q=-Dg(g+B)$. Their quotient determines $k=(g+B)p/q$. The successful sign gives $D=3$, $B=-1/3$, and

$$
z=\frac{(t+1)^3(4t-1)}{9(t-1)},\qquad
x=-\frac{(t+1)(2st-t^2-2t+1)}{3(t-1)} .
$$

The target is smooth. The cubic identity reduces to zero modulo the original projection equation. The $t$-pencil has degree three, $g$ has degree four and projection to $z$ on the target has degree three. The map therefore has degree four. Calculating $dz/(3x^2-Dz)$ on the source also verifies its canonical differential line independently of the branch matching.

The [involution atlas](../involution-atlas/README.md) finds zero involutions over the algebraic closure, so this map has no degree-two intermediate. Over $\mathbb F_5$ the first three point counts are 8, 38 and 110. The Frobenius polynomial is
$$
(T^2+3T+5)(T^4-T^3+6T^2-5T+25).
$$
The root-ratio certificate finds exactly one geometric elliptic isogeny block, of multiplicity one. Distinct optimal characteristic-zero elliptic quotients give distinct projectors under specialization. This bounds the total number of primitive elliptic directions by one, in all degrees. The remaining abelian surface factor is geometrically simple. This also determines the radical of the *unrestricted* projected incidence on this input. In the affine differential chart $[e_0:e_1:1]$ it is $(3e_0+2,3e_1+2)$. The other two charts have the unit ideal. This conclusion follows because the positive line is $[-2:-2:3]$ and the specialization certificate excludes every other geometric line. It determines the geometric solutions even though the raw Gröbner run reached its time limit. [Polarization and specialization proof](../polarized-reduction/README.md)

## Infinitely many positive inputs

For every integer $r\equiv1\pmod5$, take
$$
a_1=3r^2,\qquad a_3=-1-2r^3,\qquad
E_r:\ y^2+a_1xy+a_3y=x^3,\qquad g(t)=\frac{t^3(4t-5)}{t-2}.
$$
Form the fibre product over $y=g(t)$. An integral plane quartic model is obtained by setting $b=(t-1)s-r$, $A=a_1t(4t-5)(t-2)$, and $B_1=(4t-5)(t-2)(t^3(4t-5)+a_3(t-2))$. Its affine equation and map are
$$
F_r(t,s)=\frac{b^3-Ab-B_1}{(t-1)^2}=0,\qquad
(x,y)=\left(\frac{tb}{t-2},g(t)\right).
$$
The quotient defining $F_r$ is a polynomial with integral coefficients and plane degree four. At $r=1$, replacing $t,s$ by $X+Y+1,X-Y+2$ gives the benchmark. The verifier checks these identities symbolically before evaluating the input. The common smooth reduction at 5 ensures geometric smoothness of every member. The degree-three and degree-four function fields over the base line are linearly disjoint by their coprime degrees.

Over the algebraic closure the rational cover $g$ has ramification types $(3,1)$ at 0 and infinity, and $(2,1,1)$ at 1 and $625/4$. Its derivative is $6t^2(t-1)(2t-5)/(t-2)^2$. Its transitive monodromy group has order divisible by four, contains a 3-cycle and contains a transposition. This forces $S_4$. On $E_r$ the discriminant equation of the degree-three projection gives $-27y^2(y-1)(y-a_3^2)$. The value $a_3^2$ is a branch point outside the branch locus of $g$ for every integer in the progression. Indeed it is an integer square different from zero and one, and $625/4$ is nonintegral. The cubic extension for $E_r$ therefore cannot be contained in the splitting field for $g$. After base change the latter still has group $S_4$ because a cubic minimal polynomial could become reducible over the normal splitting field only by acquiring a root there, which would contain a conjugate of the cubic field. The natural four-letter action is primitive, proving primitivity of the exhibited map on every member.

Every member has the same multiplicity-one elliptic block at 5. The exhibited map is therefore its unique primitive geometric elliptic quotient. Its compositions with isogenies give precisely $\sigma_1(n)$ geometric classes of degree $4n$, since an elliptic curve in characteristic zero has $\sigma_1(n)$ subgroups of order $n$. Here $\sigma_1$ is the sum of positive divisors. There are no other positive degrees. The family has infinitely many geometric source isomorphism classes. Its unique elliptic target has nonconstant invariant $j=(a_1^4-24a_1a_3)^3/[a_3^3(a_1^3-27a_3)]$, so infinitely many of these targets and their sources differ.

## The rational spectrum

Restrict further to $r\equiv1\pmod{7735}$, where $7735=5\cdot7\cdot13\cdot17$. The elliptic traces at these primes are $-3,-4,-1,-3$. For every candidate prime isogeny degree in Mazur's theorem apart from 3, at least one of these Frobenius polynomials $T^2-a_pT+p$ has no root modulo that candidate, at a different residue characteristic. This excludes a rational cyclic isogeny of that prime degree. The verifier tests the roots explicitly, including the prime 2. The flex $(0,0)$ on $E_r$ gives a rational point of order three. At 7 there is no eigenvalue modulo 9, excluding a rational cyclic kernel of order nine. There is exactly one Galois-stable order-three subgroup. Two stable lines would diagonalize Frobenius at 7, whose repeated eigenvalue modulo 3 is 1. This would give nine rational 3-torsion points, while $\#E_r(\mathbb F_7)=12$. [Mazur, *Rational isogenies of prime degree*](https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf)

Every finite rational isogeny kernel with elementary divisors $n\mid m$ contains $E_r[n]$. After division by this full subgroup the remaining cyclic part has order one or three. Both possibilities occur and each subgroup is unique. The plane model has the rational point $[0:1:0]$ before the coordinate change, so maps can be normalized by rational translation. The unique primitive quotient descends over the rationals. Every rational map factors through it over the rationals by the universal property of the connected Jacobian quotient. Consequently there is exactly one rational class in each degree $4n^2$ or $12n^2$, and none in other degrees. For example the benchmark has geometric maps in degree eight and no rational map of that degree.
