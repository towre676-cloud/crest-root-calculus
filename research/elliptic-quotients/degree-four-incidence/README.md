# A projective incidence for primitive degree-four elliptic maps

For any smooth plane quartic over the rationals, the three elimination ideals in [incidence.py](incidence.py) cut out exactly the geometric primitive degree-four classes, up to target isomorphism and translation, in the projective plane of differentials. Their coefficients retain Galois conjugacy. The supplied [certificate verifier](verify_projective.py) checks explicit map formulas and their primitivity on the [previous degree-four example](../primitive-degree-four/README.md). The nontrivial unspecialized Gröbner systems were generated. Their full elimination remains unbenchmarked.

## Why three plane quartics suffice

Work over an algebraic closure of characteristic zero. Write $H=K_C=\mathcal O_C(1)$. For a degree-four map $f:C\to E$, choose any elliptic origin $O$ and put $M=4H-f^*(3O)$. Its degree is four. Riemann–Roch gives $h^0(M)\ge2$. Its fixed divisor has degree at most one. A larger fixed divisor would leave a pencil of degree at most two on a nonhyperelliptic genus-three curve. After removing the fixed part, a general section has simple zeros and avoids the support of the fixed part. For example, choose a basepoint-free pencil in the moving system and a regular value of its separable map to the line. Thus some section $s$ of $M$ has four distinct zeros.

Multiply the three coordinates of the Weierstrass embedding $E\subset\mathbf P^2$ by $s$. The result is three sections of $4H$ with precisely those common zeros. Every section lifts to a quartic form by the plane restriction sequence. The kernel on degree-four forms is spanned by $F$. The program omits the leading monomial of $F$ from its fourteen-coordinate basis for each form. This construction accommodates arbitrary source base points, target translations and ramification. [Riemann–Roch](https://stacks.math.columbia.edu/tag/0B5B)

## The incidence equations

Write $A,B,C_0$ for the three forms and work on the ordinary short cubic $y^2=x^3+ax+b$ with affine functions $x=A/C_0, y=B/C_0$. The cubic equation modulo $F$ is

$$
 B^2C_0-A^3-aAC_0^2-bC_0^3=0.
$$

For affine source coordinates $u,v$ in $W=1$, let $D=F_v\partial_u-F_u\partial_v$. A nonzero canonical form has expression $\ell(u,v,1)\,du/F_v$. The pullback of $dx/(2y)$ gives the second polynomial identity

$$
 C_0 D(A)-A D(C_0) = 2\lambda B C_0\,\ell(u,v,1) \pmod F .
$$

Here $\lambda\ne0$. In homogeneous coordinates the right side gains a factor $W$ and the two sides have degree ten. The program computes all coefficients of both remainders by exact rational monomial reduction. The charts use the last nonzero coefficient of $\ell$ as pivot.

Four distinct marked source points satisfy $A=B=C_0=0$. Four further distinct points have $C_0\ne0$ and map to the same affine point $(x_t,y_t)$ with $y_t\ne0$. All eight satisfy $F=0$ and have projective representatives normalized by a common auxiliary linear form. Another auxiliary linear form separates points within each set of four. Scalar inverses impose the nonzero conditions, including the discriminant $4a^3+27b^2\ne0$. These auxiliary choices exist for every map by choosing a general unramified fibre away from the images of the base divisor.

Conversely the differential identity excludes constant maps because $\lambda B C_0\ell\ne0$ in the function field. The image lies on the given smooth cubic. Removing the common base divisor of length $h$ on the smooth source gives $3\deg f=16-h$. The first set of marks gives $h\ge4$, hence $\deg f\le4$. The four distinct points of one fibre give $\deg f\ge4$. In particular the incidence handles an initially unknown base divisor without assuming a fixed affine chart for its support. [Degrees on curves](https://stacks.math.columbia.edu/tag/0AYQ)

## Primitivity and finite projection

All intermediate curves of a composite degree-four map would be double quotients. A rational intermediate cannot map to an elliptic curve. A genus-two intermediate would require an étale double quotient by Riemann–Hurwitz. Every involution fixes a line in the canonical plane, and that line meets the quartic, excluding this case. A genus-one intermediate gives a bielliptic involution. Its invariant differential line is represented by its projective axis. Conversely if the pullback differential is invariant under such an involution, then the difference of the two elliptic maps obtained before and after the involution has zero differential and is constant. Evaluating at a fixed point makes this constant zero. Thus the map descends to the bielliptic quotient.

The [involution atlas](../involution-atlas/README.md) computes the full set of these axes. For each orbit it returns an irreducible polynomial $q(T)$ and an axis vector $a(T)$. The incidence inverts

$$
 \operatorname{Res}_T\big(q(T),(1,s,s^2)\cdot(\ell\mathbin\times a(T))\big).
$$

If the candidate equals any conjugate axis, this resultant vanishes identically in the auxiliary scalar $s$. If all axes differ, a single choice of $s$ avoids the finitely many forbidden roots for all orbits. This removes imprimitive maps before elimination, including those associated with nonrational involutions.

A primitive class has connected kernel for its induced map $f_*:J_C\to E$. Its differential line determines the kernel's Lie algebra and hence the connected complex subtorus, so the line specifies the class uniquely. The set of lines is finite. In the endomorphism lattice of $J_C$, the symmetric integral endomorphism $T=f^*f_*$ satisfies $T^2=4T$, has rank two on rational first homology and has Rosati squared norm $\operatorname{Tr}(T^2)=32$. Positivity gives finitely many such lattice points, each determining its kernel Lie algebra. Consequently each elimination ideal gives the actual image as a set of geometric points. In particular closure of the projection adds no spurious differential directions. For any reported line the specialized original incidence has a point over a finite extension. Its form coefficients recover a representative. One can find such a point by enumerating tuples of algebraic numbers with exact equality tests, or by finite-type algebraic sampling. A line's field of definition determines its class. A specified elliptic origin and map over that field require their own descent checks. [Rosati involution and positivity, Theorem 14.3](https://www.jmilne.org/math/CourseNotes/AVc.pdf)

## Exact degree check on an explicit map

The [verifier](verify_projective.py) computes the rank of the $112\times66$ coefficient matrix of all degree-six multiples of $F,A,B,C_0$ in degree ten. Choose one form $G$ that is nonzero modulo $F$. The complete intersection $F=G=0$ has length sixteen, and its Hilbert function is sixteen from degree six on by the $(4,4)$ resolution. A linear form avoiding its support trivializes all twists there. Consequently degree-six multipliers restrict to all functions on this finite scheme. In degree ten the images of the other two forms then generate their entire ideal on the scheme. The dimension $66-\mathrm{rank}$ is the scheme-theoretic common base length, with multiplicities. This check works for explicitly supplied presentations even when their base points coalesce. Cubic identity, nonsingularity and differential are checked separately.

Run the specialized controls and produce one generic elimination job with Python and SymPy:

```sh
python verify_projective.py
python incidence.py ../primitive-degree-four/ciani.json 2 --singular quartic-chart2.sing
```

Use pivot arguments 0, 1 and 2 to cover the differential plane. Run an exported job in Singular, or supply `--solve` in place of `--singular FILE` for unrestricted exact SymPy elimination. The zero ideal on the pivot-zero chart means the unique direction is present. The unit ideal means the chart is empty. On the example, the exporter produced 143 equations and 101 auxiliary variables on each of the two positive-dimensional coefficient charts. The SymPy solver returned the unit ideal on the pivot-zero chart. Singular was unavailable for this replay.

For an arbitrary supplied rational map, run `python verify_projective.py --check SOURCE.json MAP.json`. The source uses the repository's quartic monomial format. The map JSON contains rational short-Weierstrass coefficients `a` and `b`, and a `forms` array of three quartics in target order. Each quartic is an array of `[i,j,k,coefficient]` terms. Rational coefficients can be supplied as strings.

The [recorded controls](receipt.json) verify quartic presentations for all four primitive classes on this example. Their degree is four and their common base length is four. The rational differential directions are $[1:1:0]$ and $[1:-1:0]$. The other two are $[1:0:\pm\sqrt{-2}]$ and are conjugate over $\mathbf Q(\sqrt{-2})$. The independent [endomorphism lattice calculation](../a3-spectrum/README.md) proves these four classes exhaust the list. The double quotient followed by a 2-isogeny has the same degree and base length but matches the involution axis $[1:0:0]$. The double quotient itself has degree two, base length ten and the same axis.
