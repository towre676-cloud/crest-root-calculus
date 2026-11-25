# Optimal elliptic quotients with repeated factors

A Frobenius polynomial can identify a repeated elliptic factor without identifying any one of its elliptic directions. Those directions form a lattice, and the polarization assigns degrees to them. The finite kernel of an isogeny changes the degree by an explicit, periodic factor. We determine this factor in arbitrary dimension for a non-CM isotypic abelian variety with a supplied isogeny. We then use it to give every optimal degree, every multiplicity, the field of definition, and maps for an infinite family of genus-two curves. The repeated-factor example with squared roots $1,4,16$ is its member $t=2$.

## The degree and its finite correction

Work over an algebraically closed field of characteristic zero. Let $E$ be an elliptic curve with $\operatorname{End}(E)=\mathbb Z$, let $(A,L)$ be a polarized abelian variety, and let
$$
\Psi:E^r\longrightarrow A
$$
be an isogeny. Write $K=\ker\Psi$ and choose $n\ge 1$ killing $K$. Identify $E$ with its dual using the principal polarization. Under the product identification, the polarization of $\Psi^*L$ is an integral symmetric positive definite matrix $M$.

**Theorem.** The elliptic subvarieties of $A$ are in bijection with primitive integer vectors $v\in\mathbb Z^r$ modulo sign. If $h_v:E\to E^r$ is $P\mapsto(v_1P,\ldots,v_rP)$, put
$$
K_v=\{P\in E[n]:h_v(P)\in K\}.
$$
Then the degree of $L$ restricted to the associated elliptic subvariety is
$$
d(v)=\frac{v^\mathsf{T} M v}{|K_v|}.
$$
For a Jacobian with its canonical principal polarization this is the degree of the associated optimal elliptic quotient of the curve, up to translation and target isomorphism.

To prove this, the identity $\operatorname{End}(E)=\mathbb Z$ identifies elliptic subvarieties of $E^r$ with one-dimensional $\mathbb Q$-subspaces of $\mathbb Q^r$. The primitive generator gives a closed immersion, since a Bézout linear combination of its coordinates is one. An isogeny gives a bijection on connected abelian subvarieties, with inverse the identity component of the inverse image. The kernel of $\Psi h_v$ lies in $E[n]$ by the same Bézout identity. Restriction pulls $L$ back to a degree-$v^\mathsf{T}Mv$ line bundle on $E$. Dividing by the degree $|K_v|$ of $E\to(\Psi h_v)(E)$ gives the formula. For a Jacobian, dualizing the elliptic inclusion using the principal polarizations gives a quotient with connected kernel. Composing with an Abel–Jacobi map gives the optimal elliptic map, whose degree is the degree of the restricted polarizing line bundle. See the [polarized-reduction argument](../polarized-reduction/README.md).

In particular, the correction $|K_v|$ depends only on $v\bmod n$ and the embedded finite subgroup $K$. Two finite subgroups of the same order need not give the same corrections. The calculation retains the subgroup inside $E[n]^r$.

### A counting constant

Let $r\ge2$, write $B_M=\{x\in\mathbb R^r:x^\mathsf{T}Mx\le1\}$, and let $N_A(X)$ count the elliptic subvarieties of degree at most $X$. Sum over residues $a\in(\mathbb Z/n\mathbb Z)^r$ such that $\gcd(n,a_1,\ldots,a_r)=1$. Then
$$
N_A(X)=c_A X^{r/2}
  +O_{M,n}\!\left(X^{(r-1)/2}L_r(X)\right),
\quad
c_A=\frac{\operatorname{vol}(B_M)}
{2\,n^r\,\zeta(r)\prod_{\ell\mid n}(1-\ell^{-r})}
\sum_a|K_a|^{r/2},
$$
where $\ell$ runs over prime divisors, $L_2(X)=\log(2+X)$, and $L_r(X)=1$ for $r\ge3$.

In a fixed admissible residue class, primitivity is detected by Möbius inversion on common divisors coprime to $n$. The region is $\sqrt{X|K_a|}\,B_M$. After division by $m$, the congruence is still one class modulo $n$. The lattice estimate for this ellipsoid has main term its volume divided by $n^r$ and error $O_{M,n}((\sqrt X/m)^{r-1}+1)$. Sum over $m=O_{M,n}(\sqrt X)$, excluding the zero vector, and extend the main sum to infinity. This gives the Euler product and the stated error. Dividing by two identifies signs. The expression computes the leading coefficient from the matrix and finite subgroup without enumerating curves in $A$.

## A complete spectrum over $\overline{\mathbb Q}$ and over $\mathbb Q$

For $t\in\mathbb Q,\ t>1$, let
$$
C_t:y^2=(x^2-1)(x^2-t^2)(x^2-t^4),\qquad
E_t:Y^2=(u-1)(u-t^2)(u-t^4).
$$
This curve is smooth and has two rational points at infinity. Pick one as base point. The elliptic curve has Legendre parameter $1+t^2$ and no complex multiplication. Indeed, writing $t=a/b$ in lowest terms, its $j$-invariant is
$$
j=256\,\frac{(a^4+a^2b^2+b^4)^3}{a^4b^4(a^2+b^2)^2}.
$$
No odd prime dividing $ab(a^2+b^2)$ divides the numerator. Such a prime exists since $a>b>0$. If both $a,b$ are odd, $a>1$ has an odd prime factor. If they have opposite parity, $a^2+b^2$ is odd and greater than one. Thus $j$ is not integral, whereas CM $j$-invariants are algebraic integers.

Over $k=\mathbb Q(i)$, define degree-two maps to $E_t$ by
$$
f_+=(x^2,y),\qquad f_-=(t^4/x^2,\ it^3 y/x^3).
$$
A map appearing in a group sum below is centered by subtracting its value at the chosen source point. Involution $x\mapsto-x$ acts with opposite signs on the two maps on Jacobians, so their cross-adjoint compositions are zero. Each pullback is an elliptic embedding because a degree-two cover of a smooth genus-two curve cannot factor through a nontrivial isogeny. The sum of the pullbacks gives a polarized $(2,2)$-isogeny $\Psi:E_t^2\to J(C_t)$ with $M=2I$.

For $r\in\{1,t^2,t^4\}$, the difference of the sum of the two Weierstrass points over $x^2=r$ and a hyperelliptic fibre defines a class $\Delta_r$. The pullback of the nonzero two-torsion point at $u=r$ by $f_+$ is $\Delta_r$. The analogous pullback by $f_-$ is $\Delta_{t^4/r}$. The two choices of hyperelliptic fibre used in these computations are linearly equivalent. It follows that the kernel of $\Psi$ is the graph of the transposition on $E_t[2]$ interchanging $1$ and $t^4$ and fixing $t^2$. The graph has order four, the degree of $\Psi$. For a primitive direction $(a,b)$, its kernel intersection has order two when both coordinates are odd, and order one otherwise. Hence
$$
d(a,b)=\begin{cases}
 a^2+b^2,&a,b\ \text{odd},\\
 2(a^2+b^2),&a,b\ \text{of opposite parity}.
\end{cases}
$$

Consequently, the degree $2m$ occurs exactly when every prime factor of $m$ is $1\bmod4$, allowing $m=1$. Its geometric multiplicity is $2^{\omega(m)+2}$, where $\omega$ counts distinct prime factors. Here optimality means the map does not factor through a further elliptic isogeny. The opposite-parity vectors are precisely the primitive representations of $m$ as a sum of two squares, modulo sign. The odd-odd vectors are in bijection with the same set by $(a,b)\mapsto((a+b)/2,(a-b)/2)$. For odd $m$, unique factorization in the Gaussian integers shows that primitivity forces all prime factors to split and all of each prime power to go to one conjugate. The units and the overall sign give $2^{\omega(m)+1}$ classes in each parity class. The total number of geometric classes through degree $X$ is $2X/\pi+O(\sqrt X\log(2+X))$ by the counting theorem.

Complex conjugation on the two elliptic directions sends $(a,b)$ to $(a,-b)$, because $f_-$ is conjugate to its negative on Jacobians. A rational optimal map must have a stable line, so only the two axes descend. They do descend. The second can be written over $\mathbb Q$ as $(x^2,xy)$ to the smooth pointed quartic $v^2=u(u-1)(u-t^2)(u-t^4)$, translated if necessary. Thus there are exactly two rational optimal classes, both of degree two. Every geometric optimal class has a representative over $k$ by the formulas below.

### Constructing every map

If $a,b$ have opposite parity, use the centered map $[a]f_++[b]f_-$ to $E_t$. The pullback on Jacobians is $\Psi h_{(a,b)}$, and has trivial kernel.

For the other parity class put $h=x-t^2/x$. The common target is the smooth genus-one model
$$
E_t^\star:\ v^2=u(u-2it)\bigl(u^2-(t^2-1)^2\bigr).
$$
Use the point $(0,0)$ as origin, and take the centered versions of the maps
$$
g_+=(h,\ y(x-it)/x^2),\qquad
g_-=(-h,\ -y(x+it)/x^2).
$$
Their pulled-back differentials $du/v$ are $(x+it)\,dx/y$ and $(x-it)\,dx/y$, respectively. The corresponding differentials of $f_+,f_-$ are $2x\,dx/y$ and $2it\,dx/y$. Composing the two $g$ maps with one common isogeny $E_t^\star\to E_t$ multiplies their pulled-back differential forms by the same scalar. Thus their quotient directions are $(1,1)$ and $(1,-1)$, and the pullback for
$$
[(a+b)/2]g_+ + [(a-b)/2]g_-
$$
has direction $(a,b)$. This sum is optimal. To see this last claim without assuming that a primitive slope guarantees connectedness, let $H$ be the fixed subgroup of order two of the gluing transposition. Both diagonal embeddings have source $E_t/H$. By the same kernel formula, this is also the source of the elliptic inclusion for every odd-odd direction. The degree-two maps $g_\pm$ are optimal, so they identify $E_t^\star$ with $E_t/H$ in each direction. An automorphism of this non-CM elliptic curve fixing its origin is $\pm1$. The displayed pulled-back differentials fix the relative signs. Equivalently, the two adjoint inclusions, precomposed with $E_t\to E_t/H$, are the two signed directions with the same common sign. The combination in the display therefore factors the exact homomorphism $\Psi h_{(a,b)}$ through $E_t/H$, and the resulting inclusion has trivial kernel. These are constructive maps using only four degree-two formulas and elliptic addition.

## What the torsion operator can recover

For a principally polarized $J$, suppose an optimal degree-$d$ quotient $q:J\to E_1$ has adjoint inclusion $i$. Set $E_2=\ker q$ and $U=iq=d\,e$, where $e$ is the symmetric projector. Over an algebraically closed field of characteristic prime to $d$,
$$
\operatorname{im}(U\!\mid_{J[d]})
   =i(E_1[d])=i(E_1)\cap E_2.
$$
The first equality uses surjectivity of $q:J[d]\to E_1[d]$: lift a $d$-torsion point and subtract a $d$-division point in the connected kernel. The second equality follows from $qi=[d]$. In dimension two, identify the two elliptic factors by their inclusions. The recovered subgroup determines the graph in $E_1[d]\times E_2[d]$ of the polarization anti-isometry, using $(P,-P)$ for each common point in $J$.

At a good prime where the rational projector specializes and $d$ is prime to the residue characteristic, the same formula reconstructs the gluing on that special fibre. In the multiplicity-one setting the projector is available from the Frobenius polynomial together with its action on torsion. Evaluating the integral operator $d\,e$ may require a deeper torsion level than $d$ if its expression in Frobenius has a denominator. One can use torsion lifts and check that the resulting action is independent of the lift. This statement recovers the finite subgroup in compatible torsion coordinates. Lifting it to explicit characteristic-zero quotient maps from an arbitrary source equation requires additional reconstruction. In the family above the four explicit degree-two maps supply it directly.

## The same polarized isogeny data with different spectra

The effect of the embedded finite kernel can be isolated while keeping $E_t^2$, $M=2I$, and the kernel order fixed. Write $R={1,t^2,t^4}$, $W=t^4+t^2+1$, and put
$$
B=\frac{t^4(t^2+2)}{W},\qquad
D=\frac{t^2(2t^2+1)}{W},\qquad
S=-\frac{t^4(t^4-1)^2}{W^2},\qquad
\eta^2=\prod_{r\in R}(D-r).
$$
The fractional linear map $T(u)=D+S/(u-B)$ cycles $1\mapsto t^2\mapsto t^4\mapsto1$. Define
$$
C_t^{\rm cyc}:y^2=\prod_{r\in R}(x^2+B-r).
$$
This sextic is separable since $B\notin R$ and the three roots in $x^2$ are distinct. The two maps to $E_t$ over $\mathbb Q(\eta)$ are
$$
(u,Y)=(x^2+B,y),\qquad (u,Y)=(D+S/x^2,\eta y/x^3).
$$
For the second map, the zeros of each factor $D+S/x^2-r$ occur at the squared branch value $T^{-1}(r)-B$. The product identity follows by matching the coefficient of $x^6$ after multiplying by $x^6$. Both maps have degree two. The same orthogonality and Weierstrass pullback arguments give $M=2I$ and identify the kernel as the graph of the three-cycle (or its inverse when the factors are exchanged). These formulas give a direct model of the classical two-torsion gluing construction.

The three-cycle fixes only zero on $E_t[2]$. Thus $|K_v|=1$ for every primitive direction of the companion and its optimal degrees are $2(a^2+b^2)$ with $\gcd(a,b)=1$. Both maps have a value at either point at infinity, so their centered sums construct all optimal quotients. With every prime factor of $m$ congruent to $1$ modulo $4$, the multiplicity is $2^{\omega(m)+1}$ at each of $2m$ and $4m$. These are the only degrees, including the cases $m=1$. For an odd-odd primitive direction, $(a+b)/2,(a-b)/2$ give a primitive opposite-parity representation of $m=(a^2+b^2)/2$, which proves the multiplicity at $4m$ from that at $2m$.

| Quantity | $C_t$ | $C_t^{\rm cyc}$ |
| --- | ---: | ---: |
| Geometric degree-two classes | 4 | 2 |
| Geometric degree-four classes | 0 | 2 |
| Leading constant for the number through degree $X$ | $2/\pi$ | $3/(2\pi)$ |

Over $k=\mathbb Q(i,\eta)$ the two Jacobians are isogenous via the same self-product and both pull back their canonical polarization to $2I$. Their Frobenius polynomials coincide at every place of common good reduction over $k$. At $t=2$, one may take $k=\mathbb Q(i,\sqrt{35})$ and an integral companion model $Y^2=(X^2+175)(X^2+28)(X^2-560)$. The degree multiplicities and the different asymptotic constants show how much information can change while the polarized rational isogeny data and the kernel order agree.

## Replay and antecedents

~~~sh
python -B research/elliptic-quotients/isotypic-gluing/verify.py
python -B research/elliptic-quotients/isotypic-gluing/spectrum.py --parameter 2 --bound 130
python -B research/elliptic-quotients/isotypic-gluing/spectrum.py --companion --bound 130
python -B research/elliptic-quotients/isotypic-gluing/spectrum.py --direction 3 1
~~~

The verifier checks the original four map identities and their pulled-back differentials symbolically, along with both companion map identities. It computes the two finite kernel corrections, compares primitive directions with the factorization-based multiplicities at every degree through 10,000 and compares point counts at seven split primes for $t=2$.

Kani's refined Humbert invariant already relates optimal genus-two covers to primitive representations and is computable from a Jacobian presentation. Guerra gives degree parametrizations and counting estimates for product abelian varieties. His §7.1 gives the same projection formula for degrees across a polarized isogeny and estimates the resulting counts using bounds on the restricted isogeny degree. Howe, Leprévost and Poonen give general explicit formulas for two-torsion gluing. The refinement here evaluates that restricted degree as a function of the residue of each slope, computes the resulting leading constant and uses it to distinguish two explicit families across their entire spectra. Priority for this refinement remains open. [Kani, Part I](https://doi.org/10.1007/s40316-018-0105-6) · [Kani, Part II](https://doi.org/10.1016/j.jnt.2018.05.011) · [Guerra](https://doi.org/10.13137/2464-8728/27068) · [Howe, Leprévost and Poonen](https://arxiv.org/abs/math/9809210)
