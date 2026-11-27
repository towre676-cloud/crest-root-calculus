# Equations for every parameter in the full-spectrum family

For $t\ge0$ with $t\not\equiv3\pmod5$, let $L_A(t),L_B(t)$ be the two index-three lattices in the [full-spectrum construction](../opposite-secant-spectra/delay/README.md#an-infinite-family-with-identical-full-spectra). The two curves now have a uniform determinantal equation on a product of two explicitly specified CM elliptic curves. The construction uses roots of class and division polynomials, two Vélu isogenies and a linear system. It gives equations over an explicit finite number field for each parameter. The presentation below is an affine model whose smooth completion has genus two. Expanded hyperelliptic sextics are available separately for $t=0$.

## Why the degree-three construction stops at zero

Write $Q_t=Q_0+36t(x-y)^2$ in parent coordinates. The plane $x=y=u$ has norm $4(3u^2+2uz+2z^2)$. If this were a positive square $(2v)^2$, we would have
$$
 (2z+u)^2+5u^2=2v^2.
$$
Reduction modulo five first forces $5\mid v$ and $5\mid2z+u$. Reduction of the resulting identity divided by five then forces $5\mid u$ and $5\mid z$. Infinite descent excludes every positive square on the plane. Away from the plane the base norm is at least five and $(x-y)^2\ge1$. Every elliptic quotient therefore has degree at least $\lceil\sqrt{36t+5}\rceil$. Exact counts give geometric minimum **nine** for all four curves at $t=1,2$, with two optimal classes on each curve. In particular the uniform construction needs to allow unbounded quotient degrees.

## Finding the two products

Put $N=90(1+18t)$, $r=\sqrt{-N}$ with positive imaginary part and $\mathcal O=\mathbb Z[r]$. For each child find integers $(a,b,n,m,k)$ with
$$
 0\le b<a,\quad a\mid b^2+N,\quad
 \gcd\bigl(a,2b,(b^2+N)/a\bigr)=1,\quad
 n,m,k>0,\quad nm-ak^2=1.
$$
Use $E_0=\mathbb C/\mathcal O$ and $E_I=\mathbb C/I$ for $I=\langle a,b+r\rangle$. The inclusion of $I$ into $\mathcal O$ induces the cyclic degree-$a$ isogeny $h\colon E_I\to E_0$. The degree-one polarizations on the factors are $1$ and $1/a$ in these complex lattice coordinates. Then the Hermitian form
$$
 H=\begin{pmatrix}n&k\\k&m/a\end{pmatrix}
$$
defines a positive principal polarization on $E_0\times E_I$. The lattice check computes the full integral Néron–Severi lattice by imposing integrality of alternating forms on the product basis. It applies the refined Humbert formula to $H$ and requires exact isometry to the specified child and inequivalence to the other child.

This search terminates for the whole stated family. The normalized cofactor form of each child is its Gram adjugate divided by four. For $A$, its first two diagonal entries are $648t+41$ and $648t+81$. For $B$ they are $108t+11$ and $648t+81$. In both cases the displayed integers are coprime because $t\not\equiv3\pmod5$. Indeed the first pair differs by 40 and is odd, while the second pair has difference 15 after multiplying the first entry by six, and $108t+11$ is prime to three. The cofactor form is therefore primitive, and so is its binary restriction to the first two coordinates. This also shows that the binary degree form in the primitive CM realization has content one. Its discriminant is $-360(1+18t)$, obtained from the child Gram determinant $1440(1+18t)$.

In the construction of Kır, the reciprocal form is the primitive normalization of these cofactors. The binary restriction represents infinitely many primes by the prime representation theorem for primitive positive definite binary forms. Thus one can find an odd prime avoiding the child discriminant without a bound conditional on GRH. Using such a prime in the primitive construction gives a product on the two proper ideals of this same order with polarization $(n,m,kh)$. The exact search in `family_equations.py` simply enumerates these possible presentations and verifies the resulting ternary form. The diagonal degrees $n,m$ occur in the child spectrum because each is the positive polarization degree of a coordinate elliptic subvariety. Its complete square-shell enumeration therefore includes the presentation. Conjugating one ideal or changing the sign of a factor absorbs the other choices of sign in the cyclic isogeny and the off-diagonal entry. [Kır, primitive construction and its proof](https://arxiv.org/html/2602.14319v1#S5) · [Prime representation by binary forms, Cox §9](https://www.math.utoronto.ca/~ila/Cox-Primes_of_the_form_x2%2Bny2.pdf)

The first two positive parameters admit the following certificates. $s$ is the number of columns in the defining determinantal matrix.

| $t$ | Curve | $(a,b)$ | $(n,m,k)$ | $s$ |
|---:|---|---|---|---:|
| 1 | $A$ | $(194,6)$ | $(9,345,4)$ | 14,573 |
| 1 | $B$ | $(239,21)$ | $(9,425,4)$ | 17,953 |
| 2 | $A$ | $(242,28)$ | $(9,27,1)$ | 2,690 |
| 2 | $B$ | $(227,23)$ | $(9,101,2)$ | 6,105 |

[Exact product and lattice parameters](family-equations.json)

## Algebraic equations for the factors and the graph

Take $j=j(\mathcal O)$ as the principal singular modulus of discriminant $-4N$, and set
$$
 E_0\colon y_0^2=x_0^3+A_0x_0+B_0,
 \qquad A_0=3j(1728-j),\quad B_0=2j(1728-j)^2.
$$
The modulus is the class-polynomial root labelled by the period $r$. Neither exceptional modulus $0$ nor $1728$ occurs. The two $a$-division points $P_0,R_0$ correspond to $(b+r)/a$ and $1/a$ in $\mathbb C/\mathcal O$. Their $x$-coordinates are the roots of the $a$-division polynomial determined by the normalization
$$
 \frac{x_0(P_0)}{1728-j}=\frac{2g_2(r)}{3g_3(r)}\wp\!\left(\frac{b+r}{a};r\right),\qquad
 \frac{x_0(R_0)}{1728-j}=\frac{2g_2(r)}{3g_3(r)}\wp\!\left(\frac1a;r\right).
$$
The formula follows by rescaling $\wp,\wp'/2$ to the displayed short Weierstrass model. It labels algebraic roots and introduces no transcendental coefficients in the equation. For $t=1,2$ the eight roots are also specified by [rational isolating boxes](family-torsion-roots.json), checked by ball enclosures against all points modulo sign. Either sign of each $y$-coordinate generates the same subgroup.

Apply Vélu to $\langle P_0\rangle$ to obtain $v:E_0\to E_I$. On complex lattices its target is $\mathbb C/(I/a)\simeq\mathbb C/I$. The image of $R_0$ generates $\ker\widehat v$. Apply Vélu a second time to $\langle v(R_0)\rangle$, then scale the resulting short coordinates by $(a^{-2},a^{-3})$ to obtain the dual $h:E_I\to E_0$. With the first and second Vélu maps normalized on differentials, the rescaled composite is multiplication by $a$. All coefficients are now in the field generated by $j$ and the two labelled $a$-division points. [Vélu construction and rational map conventions](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Elliptic_curves.html#ellisogeny)

## The determinantal equation

Write $D_0=\{O_0\}\times E_I$, $D_I=E_0\times\{O_I\}$ and $\Gamma_h=\{(P,Q):P=h(Q)\}$. On the product let
$$
 D=(n+k)D_0+(m+ak)D_I-k\Gamma_h.
$$
The class of $\Gamma_h$ comes from the morphism $(P,Q)\mapsto P-h(Q)$. In the Hermitian convention above it has diagonal $(1,1)$ and off-diagonal entry $-1$. The class of $D_I$ has second diagonal entry $1/a$. Thus $D$ induces exactly $H$. The relation $nm-ak^2=1$ implies $D^2=2$ and the line bundle of $D$ has precisely one independent global section. Its theta divisor is smooth of genus two because the child Humbert form does not represent one. The Jacobian of that divisor, with its principal polarization, is the specified product. [Theta divisors and principal polarizations](https://www.math.stonybrook.edu/~cschnell/pdf/notes/abelian-varieties.pdf) · [Kani's criterion and refined Humbert correspondence](https://doi.org/10.1007/s40316-018-0105-6)

Here is the equation entirely in algebraic operations. For each elliptic factor use the monomials $x^iy^e$ with $i\ge0$, $e\in\{0,1\}$ and $2i+3e$ at most its pole bound. Take the products of the two lists at bounds $n+k$ and $m+ak$, giving $s=(n+k)(m+ak)$ monomials $z_i(P,Q)$. Let $M$ be the matrix of these $s$ columns. For every $0\le j<k$, differentiate each monomial $j$ times in the first factor by the invariant derivation
$$
 \partial=2y_0\,\frac{\partial}{\partial x_0}+(3x_0^2+A_0)\,\frac{\partial}{\partial y_0}.
$$
Then set $P=h(Q)$, reduce powers of $y_I$ using the second Weierstrass equation, clear a common denominator, and take the coefficients of the resulting polynomial pairs in $x_I$. Concatenate these coefficient rows to make $M$. This imposes vanishing of order at least $k$ on $\Gamma_h$, because $\partial$ is transverse to that graph. The kernel is the global section space of $D$, so $\operatorname{rank}M=s-1$.

Away from $D_0,D_I$ and the graph, the two elliptic equations and the condition
$$
 \boxed{\quad \operatorname{rank}\begin{pmatrix} M\\z_1(P,Q)\;\cdots\;z_s(P,Q)\end{pmatrix}<s\quad}
$$
are the affine equations of the required genus-two curve. Equivalently, use all maximal minors, or use the alternating maximal minors of any $s-1$ independent rows of $M$ as the coefficients of a single linear combination of the $z_i$. One can impose the open conditions by adjoining inverses. For example, in affine elliptic coordinates invert the graph-map denominator and the difference $x_0(P)-x_0(h(Q))$. This excludes both signs of the graph and leaves a dense open of the theta divisor. Taking its smooth completion gives precisely the stated curve. The known equality of square representation counts then proves equality of the quotient spectra of these explicit models.

## Replay and format

~~~sh
study=research/elliptic-quotients/full-spectrum-curves
python -B "$study/family_equations.py" --isolate-torsion
python -B "$study/family_equations.py" --parameter 4
~~~

The first command rechecks the principal polarizations, distinct exact Humbert forms, degree-nine minima and all eight isolated root labels for $t=1,2$. It also constructs the section matrix on two exact controls, one with a degree-eleven isogeny and a single vanishing order, and one with a degree-two isogeny and two vanishing orders. Both have one-dimensional kernel. The second command demonstrates the parameter search beyond the four stored models. The displayed determinants specify the large matrices without expanding their coefficients or minors. Expanded sextics for positive $t$ and efficient algorithms for that expansion remain open.
