# Reconstructing optimal elliptic degrees from a bielliptic equation

Start with the integer coefficients of a squarefree quintic or sextic \(f\) defining \(C:y^2=f(x)\). This study recovers the full **geometric** optimal elliptic degree spectrum when \(C\) admits a geometric reduced elliptic involution and its two elliptic factors are geometrically isogenous and non-CM. The implementation certifies non-CM by rational singular-moduli exclusion, nonintegrality or distinct ordinary reductions. It computes the isogeny and the relevant action on 2-torsion exactly. The proof includes even and composite minimum isogeny degrees.

A reduced elliptic involution acts on the line of hyperelliptic fibres. A rational reduced involution can have lifts that are only defined over a quadratic field. The solver constructs all geometric reduced elliptic involutions exactly and recovers the elliptic factors from rational and nonrational involutions. The outputs describe maps over the algebraic closure, including cases where the elliptic factors are defined over a number field.

## Recovering the finite data

Let \(\sigma\) be a reduced elliptic involution and \(\varphi\) a degree-two invariant function on the hyperelliptic line. Each of its three orbits of Weierstrass points has an image on the quotient line. Write \(L(u)\) for the monic cubic with these values as roots. The two elliptic double covers of the quotient line have branch sets consisting of these three roots and one of the two images \(P_1,P_2\) of the fixed points of \(\sigma\).

Take \(P_i\) as the branch value at the origin of the \(i\)-th elliptic curve. The roots of the model with the origin at infinity then have abscissae \(1/(u-P_i)\) for \(L(u)=0\). Over the algebraic closure one can scale the vertical variable to obtain a monic cubic. The nonzero 2-torsion point labelled by \(u\) on either factor pulls back to the same class in the Jacobian, the sum of the two Weierstrass points in the orbit minus a hyperelliptic fibre. Consequently the kernel of
\[
\Psi_0:E_1\times E_2\longrightarrow J(C)
\]
is the graph of the label bijection \(\gamma:E_1[2]\to E_2[2]\). The pullback of the canonical polarization is twice the product polarization. These statements use the classical degree-two splitting. They also follow from the degree of the sum of the two pullbacks, which is four by their orthogonality under the involution, and the four points in the graph.

For a non-CM elliptic pair, the geometric homomorphism group has rank one if the factors are geometrically isogenous. Its generator \(\rho\), unique up to sign, is cyclic. A noncyclic kernel contains the full \(\ell\)-torsion for some prime \(\ell\) and would make \(\rho\) divisible by \(\ell\). Over any common field of definition, Galois can act on the generator only by sign. Thus its kernel is defined over that field and the quotient by the kernel has \(j\)-invariant \(j(E_2)\). This explains why an exact enumeration of isogenies from a model of \(j(E_1)\) over the common field, with target recognized by \(j(E_2)\), finds the minimum geometric degree \(k\), even when the original elliptic models require a twist. For rational \(j\)-invariants we use a rational model. For nonrational invariants the base is an explicit number field containing both fixed labels.

Here are the two exact calculations that replace the numerical reconstruction.

1. For \(\sigma(x)=(Ax+B)/(Cx-A)\), put \(\delta=A^2+BC\). Solve \(F(Ax+B,Cx-A)=\delta^3 F(x,1)\) for the homogenized binary sextic \(F\), in the disjoint charts \(C=1\) and \(C=0,A=1\). The first chart imposes invertibility with a variable \(z\) and \(z\delta=1\). This eigenvalue condition is exactly the order-two lift condition. Indeed the trace-zero matrix squares to \(\delta I\), and a lift with numerator multiplier \(\eta\) on \(y\) squares to \(\eta^2/\delta^3\). The stabilizer of the six distinct branch points is finite and reduced in characteristic zero. After splitting the points, every permutation has at most one projective transformation and three point conditions determine it scheme-theoretically. The trace and eigenvalue conditions select a subset of this reduced finite scheme. The quotient ring dimension therefore counts the geometric solutions, and rational root extraction in lexicographic order gives every rational solution. A linear form \(A+cB\) separating the remaining solutions gives a primitive element. Reducing \(A\) and \(B\) in that presentation produces exact number-field coefficients for each conjugacy orbit.
2. Compute \(L\) and the fixed-image quadratic by resultants on the invariant line. Use PARI's exact isogeny enumeration when the \(j\)-invariants differ. The program checks the degree of the returned horizontal rational function and substitutes it into the two cubic equations using its derivative to check the vertical coordinate up to a nonzero twist scalar. Equal non-CM \(j\)-invariants give an isomorphism, with horizontal scale determined by the ratio of the short cubic coefficients. Evaluate the transported isogeny on \(1/(u-P_1)\) in the finite algebra \(K[u]/(L(u))\). Convert the image abscissa back to its label using \(P_2+1/x_2\). A point sent to infinity has label \(P_2\), outside the roots of \(L\). A polynomial gcd counts the fixed nonzero labels. The exact computation needs neither a splitting field for \(L\) nor period approximations.

If \(k\) is odd, \(\gamma^{-1}\rho\) permutes the three nonzero points. Let \(s\) be the number of nonzero fixed labels. Identity, transposition and three-cycle give \(s=3,1,0\). If \(k\) is even, the same endomorphism of \(E_1[2]\) has rank one. In this case \(s=0\) when \(\gamma\) sends its kernel line to the image line in \(E_2[2]\), and \(s=1\) otherwise. This follows because a rank-one endomorphism of a two-dimensional space over \(\mathbb F_2\) fixes its image exactly when its image differs from its kernel.

### Involutions over number fields

For a representative of each nonrational Galois orbit, use its primitive-element field in the two resultant computations. If the fixed-image quadratic is irreducible over that field, its two roots lie in the same quadratic extension. Eliminate the two generators with a separating linear form and recover their exact images in a single primitive-element field. The transfer of the branch labels uses this explicit field map. The isogeny check then uses full field coefficient vectors and verifies the exact pullback identity in that field. Conjugate involutions give the same geometric spectrum.

A rational \(j\)-invariant is CM precisely when it belongs to the list of thirteen rational singular moduli. A nonintegral \(j\) is non-CM. For an algebraic integral \(j\), two degree-one good ordinary reductions with distinct geometric endomorphism fields also prove non-CM, since a CM field would inject into both quadratic fields. The implementation searches a finite range for this last certificate and can abstain. [Rational singular moduli](https://aif.centre-mersenne.org/item/10.5802/aif.3520.pdf) · [Ordinary Frobenius fields](https://arxiv.org/abs/1610.09674)

## Degree law for every cyclic minimum isogeny degree

Every elliptic subvariety corresponds to a primitive pair \((a,b)\) modulo simultaneous sign. Its direction in \(E_1\times E_2\) is \(h(P)=(aP,b\rho P)\). The degree of the associated optimal map from \(C\) is
\[
d(a,b)=\frac{2(a^2+kb^2)}{|K_{a,b}|},\qquad K_{a,b}=h^{-1}(\operatorname{graph}\gamma).
\]
Indeed \(h^*\Psi_0^*\Theta\) has degree \(2(a^2+kb^2)\). Dividing by the degree of the isogeny \(E_1\to (\Psi_0 h)(E_1)\) gives the restricted polarization degree, which is also the optimal quotient degree. The kernel lies in \(E_1[2\gcd(a,k)]\), since \(2aP=0\), \(2b\rho P=0\), and the dual isogeny gives \(2bkP=0\). More sharply \(|K_{a,b}|\le4\gcd(a,k)\), because the kernel of \(h\) has order \(\gcd(a,k)\) and the graph has order four. The enumeration uses \(a^2+kb^2\le2kX\) for maps of degree at most \(X\).

Write \(k=2^m k_o\) with \(k_o\) odd. The odd part of \(|K_{a,b}|\) is \(\gcd(a,k_o)\). Indeed an odd-order point lies in the graph condition exactly if \(aP=0=b\rho P\). Primitivity then forces \(\rho P=0\), and the kernel of \(\rho\) is cyclic.

For odd \(k\) the 2-part is \(1+s\) when both entries are odd and one otherwise. For even \(k\) and odd \(a\) it is \(1+s\) when \(b\) is odd and one otherwise. For even \(k\) and even \(a\), put \(d=v_2(a)\) with \(v_2(0)=\infty\) and \(L_2=\min(d,m)\). Its 2-part is
\[
2^{L_2+\mathbf 1_{\{d<m,\ s=0\}\,\cup\,\{d=m,\ s=1\}}}.
\]
To see this, use Smith coordinates in which \(\rho=\operatorname{diag}(1,k)\) on the torsion. Primitivity gives odd \(b\), so the first coordinate has denominator dividing two. Write the second as \(w/2^{L_2+1}\). The first graph equation determines the first coordinate from the parity of \(w\), and the second reads
\[
\bigl(s(a/2^{L_2})+k/2^{L_2}\bigr)w=0\pmod2 .
\]
Here \(s\) is the second target coordinate of \(\gamma(0,1)\), since the first target line is the image of \(\rho\) on 2-torsion. The coefficient is zero in precisely the two cases in the exponent. This computation also shows why the isogeny degree and the order of the overall gluing group are insufficient when \(k\) is even.

The geometric optimal map in this direction is obtained by duality. Choose a base point over the algebraic closure and let \(f_i\) be the centered degree-two quotient maps. The map \([a]f_1+[b]\widehat{\rho}f_2\) has pullback \(\Psi_0h=i\circ p\), where \(p:E_1\to E_1/K_{a,b}\) is the quotient isogeny and \(i\) is an inclusion in the Jacobian. Dualizing gives the centered sum as \(\widehat p\) after the optimal map with target \(E_1/K_{a,b}\). The solver gives the degree and torsion certificate. It does not emit rational functions for every optimal map.

## Counting constant

Write \(N(X)\) for the number of geometric optimal classes through degree \(X\). The [finite-kernel counting theorem](../isotypic-gluing/README.md) applies to \(\Psi_0\circ(1,\rho)\) and the matrix \(2\operatorname{diag}(1,k)\). By the Chinese remainder theorem its residue average factors prime by prime. The leading constant is
\[
N(X)=cX+O_k(\sqrt X\log(2+X)),\quad
c=\frac{B_m(s)}{2\pi\sqrt k}
   \prod_{p^e\parallel k_o}\left(1+\frac{e(p-1)}{p+1}\right),
\qquad
B_m(s)=
\begin{cases}
3+s & m=0,\\
2m+2 & m>0,\ s=0,\\
m+5 & m>0,\ s=1.
\end{cases}
\]
For an odd prime \(p\), conditioned on an admissible pair at \(p\), the probability of \(p^t\mid a\) for \(1\le t\le e\) is \(p^{1-t}/(p+1)\). This gives the odd factor from the successive increments of \(\gcd(a,p^e)\). At 2, average over the three equally likely admissible parity types. For even \(a\) and odd \(b\), conditional \(v_2(a)=d\) has weight \(2^{-d}\) for \(d\ge1\). If \(m>0\) the average of the 2-part on that class is \(2m\) for \(s=0\) and \(m+2\) for \(s=1\). Adding the other two classes gives \(B_m(s)/3\). The geometric factor from the ellipse area and primitive lattice density is \(3/(2\pi\sqrt{k})\).

For odd prime \(k\), this specializes to \((3+s)\sqrt{k}/(\pi(k+1))\). At \(k=1\) a transposition gives \(2/\pi\) and a three-cycle gives \(3/(2\pi)\). Identity at \(k=1\) is excluded for a smooth genus-two curve because the direction \((1,1)\) would have degree one.

## Exact replay and scope

| Isogeny degree | Equations in the replay | Checks |
| ---: | ---: | --- |
| 1 | 2 sextics and 1 quintic reparametrization | Both PR #17 spectra, including one involution whose quotient models use a quadratic field |
| 3 | 6 scrambled sextics | Every gluing bijection, all three cycle types, and an identity example with three geometric reduced involutions and only one rational one |
| 2 and 4 | 4 scrambled sextics | Both rank-one gluing types at each degree |
| Nonisogenous control | 1 sextic | Distinct conjugate quadratic \(j\)-invariants and exact exclusion of an isogeny |

The verifier checks the archived multiplicities through degree 400 for the original eight inputs. For each of 14 isogeny degrees it also enumerates the finite graph independently in Smith coordinates, comparing both the local kernel order and the leading constant with the closed formulas. The degree-two multiplicities agree with twice the independently computed number of geometric reduced elliptic involutions.

Install SymPy and PARI/GP. Use the following commands from this directory.

~~~sh
python -B verify.py
python -B reconstruct.py '[-64,0,84,0,-21,0,1]' --bound 400
~~~

One can pass an alternate GP executable path via the option **--gp**. The historical high-precision route and receipts are archived in [numerical.py](numerical.py), [reconstructions.json](reconstructions.json) and [constants.json](constants.json). With the optional cypari2 dependency, [construct.py](construct.py) regenerates the six labelled \(k=3\) inputs.

Two fresh inputs have **only nonrational reduced elliptic involutions**. The quintic \(y^2=x^5+x^3+2x\) has one orbit of size two and uses a quartic field to label the fixed points. It gives \(k=1\) and the transposition spectrum. The sextic \(y^2=x^6+x^3+2\) has one orbit of size three and uses a degree-six field for the labels and its degree-three isogeny. It gives the identity spectrum. The sextic \(y^2=x^6+x^3+1\) tests integral rational \(j\)-values \(6912,790272\), both certified non-CM.

Curves without elliptic involutions need a different first quotient. CM spectra use the separate [CM method](../cm-spectra/README.md). The output here concerns geometric optimal classes. The [cyclic cubic case](cubic-genus-two.md) also decides the rational spectrum on that family.

**Antecedents.** Classical two-torsion gluing supplies the product presentation, Kani's refined Humbert invariant characterizes geometric elliptic subcovers, and Guerra relates the polarization degrees across an isogeny. The computation here evaluates the finite graph directly from the coefficient input, including the two even-degree cases and their different leading constants. Literature priority for these specific formulas remains open. [Howe, Leprévost and Poonen](https://arxiv.org/abs/math/9809210) · [Kani I](https://doi.org/10.1007/s40316-018-0105-6) · [Kani II](https://doi.org/10.1016/j.jnt.2018.05.011) · [Guerra](https://doi.org/10.13137/2464-8728/27068)
