# Inert prime gluings of a CM elliptic square

At an odd inert prime, one finite-field parameter determines the refined Humbert lattice of a nonscalar graph gluing. A single quadratic character determines its genus. The other parameter can change the geometric genus-two curve without changing the degree of any elliptic direction in the common source square.

## Statement

Let $F$ be an imaginary quadratic field, $O$ its ring of integers, $E=\mathbb C/O$, and $\ell$ an odd inert prime. Identify $E[\ell]$ with $k=O/\ell O=\mathbb F_{\ell^2}$ by dividing elements of $O$ by $\ell$. Conjugation on $k$ is denoted by a bar, and $N$ and $\operatorname{Tr}$ on $k$ are its norm and trace to $\mathbb F_\ell$.

Every $\mathbb F_\ell$-linear map $k\to k$ has a unique expression $\gamma_{a,b}(t)=at+b\bar t$. In the two embedding coordinates after extension to $k$, its matrix has rows $(a,b)$ and $(\bar b,\bar a)$, giving determinant $N(a)-N(b)$. In dimension two the determinant is the multiplier on alternating pairings, so $\gamma_{a,b}$ is an anti-isometry exactly when $N(b)=1+N(a)$. We take $b\ne0$. Write $A_{a,b}=E^2/\operatorname{graph}(\gamma_{a,b})$ with its descended principal polarization, the one pulling back to $\ell I$ on $E^2$.

The refined Humbert lattice of $A_{a,b}$ is
$$
L_a=\{(u,c)\in\mathbb Z\oplus O:\ u+\operatorname{Tr}(ac)=0\bmod\ell\},
\qquad q_0(u,c)=u^2+4N_{F/\mathbb Q}(c).
$$
The trace in the congruence is taken after reduction. The isometry classes and genera at fixed $F,\ell$ are given by
$$
\begin{aligned}
 L_a\simeq L_{a'} &\quad\Longleftrightarrow\quad
 a'\in O^\times a\ \cup\ O^\times\bar a \pmod\ell,\\
 L_a\text{ and }L_{a'}\text{ have the same genus}
 &\quad\Longleftrightarrow\quad
 \left(\frac{1+N(a)}{\ell}\right)=\left(\frac{1+N(a')}{\ell}\right).
\end{aligned}
$$
All these polarized surfaces are Jacobians. Their polarized isomorphisms are exactly those induced by separate unit changes on the two source factors and exchange of the two factors. In particular, the actions on the parameters are
$$
 (a,b)\longmapsto ((v/u)a,\,vu b),\qquad
 (a,b)\longmapsto(-\bar a,b), \qquad u,v\in O^\times.
$$
Here $u,v$ are the two units in the unit-change formula.

### Integral lattice and genus

The analytic lattice of the quotient is
$$
\Lambda_{a,b}=O^2+\ell^{-1}\widetilde{\operatorname{graph}(\gamma_{a,b})},
$$
where the tilde denotes its integral inverse image in $O^2$. Under the product principal polarization, the integral symmetric endomorphisms of $E^2$ are the Hermitian matrices with entries in $O$. On the universal cover their alternating forms are $E_B(x,y)=E_I(x,By)$, with $E_I$ the Riemann form of the product. A Néron–Severi class on the quotient is represented by such a matrix $B$ whose form is integral on $\Lambda_{a,b}$.

The pairings of lifted graph elements with $O^2$ imply that $B$ mod $\ell$ kills the graph. The graph spans $k^2$ as a $k$-space. Otherwise it would be contained in a $k$-line, and its first projection would identify it with the graph of multiplication by a scalar. Thus $B=\ell H$, with $H=\begin{pmatrix}\alpha&c\\\bar c&\delta\end{pmatrix}$ integral Hermitian. Integrality on two lifted graph elements now requires the alternating form for $H$ to vanish on the graph mod $\ell$. For graph elements $(x,\gamma x)$ and $(y,\gamma y)$, the diagonal terms contribute $(\alpha-\delta)$ times the original pairing of $x,y$, because $\det\gamma=-1$. The off-diagonal terms in the Hermitian pairing are
$$
 (ca+\overline{ca})\bar x y+cb\,\bar x\bar y+\overline{cb}\,xy.
$$
The last two terms sum to a scalar fixed by conjugation and contribute zero to the alternating pairing. Therefore the remaining condition is $\alpha-\delta+\operatorname{Tr}(ac)=0$ mod $\ell$.

The canonical class corresponds to $H=I$. Under pullback by the quotient isogeny, intersections of two divisor classes are multiplied by $\ell^2$. The intersection quantities downstairs for the class of $H$ are $D\cdot\theta=\alpha+\delta$ and $D^2=2\det H$. Consequently its refined Humbert value $(D\cdot\theta)^2-2D^2$ is $(\alpha-\delta)^2+4N_{F/\mathbb Q}(c)$. Quotient by $\mathbb Z I$ to obtain the displayed ternary lattice.

For the isometry and genus statements, write $L_0^{\rm amb}=\mathbb Z\oplus O$ with the form $q_0$. With respect to the symmetric bilinear pairing for $q_0$, the normal vector of the defining hyperplane mod $\ell$ is $(1,\bar a/2)$. Its norm is $1+N(a)\ne0$ in $k$. The ambient lattice is unimodular over $\mathbb Z_\ell$. Lifting the normal vector gives an orthogonal splitting into a unimodular plane and the normal line. The local lattice $L_a$ leaves the plane unchanged and multiplies the normal line by $\ell$. The $\ell$-primary discriminant group is cyclic of order $\ell^2$, so it has a unique subgroup of order $\ell$. It follows that the ambient lattice is the unique local integral unimodular overlattice of $L_a$ of index $\ell$. At every other prime, $L_a$ is the ambient lattice.

An isometry $L_a\simeq L_{a'}$ therefore extends to an isometry of the ambient lattice, by uniqueness at $\ell$ and equality at all other primes. Its norm-one vectors are exactly $\pm(1,0)$, so an ambient isometry acts by a sign on the first coordinate. On $O$ it is multiplication by a unit, possibly followed by conjugation. This follows by realizing the norm form on the complex plane and observing that the image of $1$ must be a unit. Acting on the hyperplane condition proves the isometry criterion. For the genus, the unimodular plane has determinant squareclass equal to $\det(q_0)/(1+N(a))$. At odd $\ell$, the rank and determinant squareclass determine each unimodular block up to isometry. The plane and scaled line give the two local blocks, and the other primes contribute no variation. Both signs occur because the norm map of $k^\times$ is onto.

### Polarized isomorphism

The form $q_0$ has ambient norm-one vectors only at $(\pm1,0)$, which are excluded from $L_a$. Thus the principal polarization has no degree-one elliptic subvariety, by the refined Humbert correspondence. A principally polarized complex abelian surface is either a product with its product polarization or a genus-two Jacobian. A product would have degree-one factors, so this is a Jacobian.

The $O$-span of $\Lambda_{a,b}$ is $\ell^{-1}O^2$, since the graph spans $k^2$ over $k$. A polarized isomorphism between two such quotients is complex linear, has entries in $F$ because it sends the standard basis into the target lattice, and preserves the $O$-span. It is therefore an integral unitary matrix for $I$. Each of its columns has squared norm one. Norms of nonzero entries in $O$ are positive integers, so the matrix is monomial with unit entries. Conversely, each such matrix taking the first graph to the second gives a polarized isomorphism. Unit changes replace $\gamma$ by $v\gamma u^{-1}$, and exchanging the factors replaces it by $\gamma^{-1}$, which yield the displayed actions.

### What every degree can see

The proof works at the level of degrees with each elliptic direction labelled by its line in the common source $E^2$. Let $\rho$ be the residue of the slope in $k$, whenever the first coordinate is nonzero in the primitive residue. The restricted kernel has order $\ell$ exactly when $\gamma_{a,b}(t)=\rho t$ for a nonzero $t\in k$. Equivalently,
$$
 N(\rho-a)=N(b)=1+N(a)
 \quad\Longleftrightarrow\quad
 N(\rho)-\operatorname{Tr}(\rho\bar a)=1.
$$
Indeed, $(\rho-a)t=b\bar t$, and the values of $\bar t/t$ are the norm-one subgroup of $k^\times$. When this condition holds the solution space has dimension one over $\mathbb F_\ell$. For all other directions, including the direction with first residue coordinate zero, the restricted kernel has order one. This argument also applies to a primitive direction in a nonprincipal ideal class by the residue transfer in the [CM degree theorem](../cm-spectra/README.md).

It follows that the entire labelled degree function depends only on $a$. Conversely, it determines $a$ with these source labels fixed. Slopes $(1,z)$ for integral lifts $z\in O$ give the full set of residues satisfying the displayed norm equation. If two scalar parts had the same set, subtracting their trace equations would put its $\ell+1$ elements on a line through zero over $\mathbb F_\ell$, which has only $\ell$ elements. Hence the scalar parts agree. The finite-kernel calculation also shows that the leading counting constant is common across the whole nonscalar family, since each scalar part has exactly $\ell+1$ discounted residue slopes out of $\ell^2+1$.

### Separation from scalar gluings

When $b=0$, the anti-isometry condition is $N(a)=-1$. The $O$-stable quotient lattice has basis $(1,\widetilde a)/\ell$, $(0,1)$ over $O$, for any integral lift $\widetilde a$. In this basis the polarization is the integral Hermitian form
$$
 \begin{pmatrix}
 (1+N_{F/\mathbb Q}(\widetilde a))/\ell & \overline{\widetilde a}\\
 \widetilde a & \ell
 \end{pmatrix}
$$
of determinant one. This gives the scalar input to a lattice checker even when its polarization is a product.

For a scalar gluing the restricted kernel has order $\ell^2$ at the single residue slope $\rho=a$, and order one elsewhere. Its total weight in the leading count is $\ell^4+\ell^2$. In the nonscalar case that weight is $\ell^2(\ell+1)+\ell^2-\ell$, with the same volume and primitive-density factors. The scalar weight is strictly larger for every odd $\ell$. Thus no scalar genus-two gluing shares the entire spectrum of a nonscalar one at fixed $F,\ell$.

### The unmarked spectrum determines the genus at fixed level

Put $\varepsilon=\bigl(\frac{1+N(a)}{\ell}\bigr)$. Let $N(X)$ count all geometric optimal classes of degree at most $X$ and $N_{\ell\nmid d}(X)$ count those with degree coprime to $\ell$. Then
$$
 \lim_{X\to\infty}\frac{N_{\ell\nmid d}(X)}{N(X)}
   =\frac{\ell(\ell+\varepsilon)}{\ell^2+2\ell-1}.
$$

To prove this, all undiscounted slopes have degree divisible by $\ell$. Among the $\ell+1$ discounted slopes, a finite residue $\rho$ has degree divisible by $\ell$ exactly when $N(\rho)=-1$. Indeed, the product polarization before scaling by $\ell$ gives norm $N(u)+N(v)$ on a principal primitive direction. For a nonprincipal direction, choose its coordinate ideal $\mathfrak a$ prime to $\ell$, and divide this norm by $N(\mathfrak a)$ as in the [CM degree theorem](../cm-spectra/README.md). The ideal residue transfer changes the norm only by a unit mod $\ell$. Thus divisibility at a finite slope still means $1+N(\rho)=0$.

For $a\ne0$ the two residue equations are equivalent to
$$
 \operatorname{Tr}(\rho\bar a)=-2,\qquad N(\rho\bar a)=-N(a).
$$
The discriminant of a quadratic with this trace and norm is $4(1+N(a))$. If it is a nonsquare, the two conjugate roots in $k$ give the two desired residues. If it is a nonzero square, the roots lie in the base field and have the wrong field trace. For $a=0$, the two norm equations are inconsistent, and $\varepsilon=1$. Hence exactly $1-\varepsilon$ discounted slopes have degree divisible by $\ell$.

In the proof of the CM counting theorem, restricting the admissible residue sum to degrees coprime to $\ell$ keeps the same common density and volume factor. Each residue slope has $\ell^2-1$ primitive pairs mod $\ell$, and the weight in dimension two is $\ell^2$ for a discounted slope and one otherwise. The restricted sum has weight $\ell^2(\ell+\varepsilon)$, and the unrestricted sum has weight $\ell^2(\ell+1)+(\ell^2-\ell)$. Their ratio proves the claim. Thus any collision of unmarked spectra for a fixed $F,\ell$ must be between forms in the same genus.

### Number of curves in a Humbert class

When $O^\times=\{\pm1\}$ and $a$ is nonreal, its norm equation gives $\ell+1$ possible $b$ and none of the polarized isomorphisms above identify distinct $b$ at that fixed $a$. For real $a$ there are $(\ell+1)/2$ classes from identifying $b$ with $-b$. By Torelli these are counts of nonisomorphic curves. In this unit case the number of Humbert isometry classes at level $\ell$ is
$$
 \frac{\ell^2+\ell-1-\left(\frac{-1}{\ell}\right)}{4}.
$$
Indeed, the admissible scalar parts have count $\ell^2-\ell-1$. The actions $a\mapsto-a$, $a\mapsto\bar a$ and $a\mapsto-\bar a$ have, respectively, $1$, $\ell-1-\bigl(\frac{-1}{\ell}\bigr)$ and $\ell$ fixed admissible elements. In the last case the trace-zero line has no element of norm $-1$, because the square of a nonzero trace-zero element is a nonsquare in the base field. Burnside's formula gives the count.

## Exact checks and handoff

~~~sh
python -B research/elliptic-quotients/inert-gluing/verify.py
python -B research/elliptic-quotients/inert-gluing/verify.py --order -7 --prime 5 --forms
~~~

The replay evaluates the alternating forms on the actual finite graph, compares the set of restricted kernels with the norm equation, and prints the Gram matrices in the basis $(\ell,0),(-\operatorname{Tr}(a),1),(-\operatorname{Tr}(a\omega),\omega)$ when requested, with primitive counts indexed from degree one. It checks ten pairs of order and prime, covering discriminants $-3,-4,-7$. Its $j=1728$ checks use odd inert levels, so they supplement the ramified level-two replay.

For discriminant $-7$ and $\ell=5$, take $\omega=(1+\sqrt{-7})/2$. The 114 nonscalar graph presentations give **36 geometric curves, seven Humbert isometry classes and two genera**. The parameter $a=\omega$ has six curves with identical Humbert form and identical degrees even with their source directions labelled. The parameters $a=0$ and $a=2+3\omega$ give inequivalent forms in the same genus. They have **two and zero** degree-two quotients respectively. Indeed $N(2+3\omega)=28$, so $1+N(a)$ is a square mod $5$ for both. A vector of value $4$ must have either $u=\pm2,c=0$ or $u=0,c=\pm1$. The first pair fails the congruence for both forms. The second pair succeeds just when $\operatorname{Tr}(a)=0\bmod5$.

At $D_F=-4,\ell=19$, the inequivalent forms with $a=3i$ and $a=7i$ have the same genus and the same primitive square counts through degree $25$. They first differ at degree $26$, with **four versus zero** classes. Thus a cutoff at the gluing prime can miss the first difference. In the ten full controls, every pair of distinct same-genus forms separates at some degree at most $2\ell$.

### Uniform level for the square comparison

In the convention $q(x)=x^{\mathsf T}Gx$, the theta level is the least positive integer $M$ for which $M(2G)^{-1}$ is integral with even diagonal. The ambient lattice has level $4|D_F|$. Indeed, for an integral basis $1,\omega$ with trace $T$ and norm $U$, the binary part of $(2G)^{-1}$ is
$$
 \frac{1}{4|D_F|}\begin{pmatrix}2U&-T\\-T&2\end{pmatrix}.
$$
At primes other than $\ell$ the lattice is unchanged. At $\ell$ the unimodular normal line gains a factor $\ell^2$ in its Gram matrix, so the exact level for every $L_a$ is
$$
 M=4|D_F|\ell^2.
$$
The checker confirms this directly from the inverse Gram matrix. Two forms in the same genus have a common theta character and their theta difference is cuspidal. Its $t=1$ Shimura lift is a holomorphic weight-two form at level $M/2$, and may be compared at the common level $M$. The coefficient of index $m$ in the lift is a triangular divisor sum of the square-index coefficients through $m^2$, with diagonal coefficient one. A conservative weight-two Sturm index is
$$
 B=\left\lfloor \frac{M}{6}\prod_{p\mid M}(1+1/p)\right\rfloor .
$$
The lift of the cuspidal difference has zero constant coefficient in this expansion. Equality of ordinary square counts through $B$ therefore proves equality in all degrees. Primitive counts are equivalent to ordinary square counts by taking out a common divisor of the vector. [Shimura lift and level](https://doi.org/10.1017/S002776300001895X)

Kani's work gives the general graph-to-Humbert calculation and the correspondence with primitive square representations. Kani and Kır count curves with a given refined Humbert invariant across the CM product locus. The scope here is the explicit classification at one fixed inert level, including what the degree function of the marked source directions retains. [Kani I](https://doi.org/10.1007/s40316-018-0105-6) · [Kani II](https://doi.org/10.1016/j.jnt.2018.05.011) · [Kani and Kır](https://www.queensu.ca/academia/kani/sites/acadkaniwww/files/uploaded_files/Preprints/cardHq16.pdf)
