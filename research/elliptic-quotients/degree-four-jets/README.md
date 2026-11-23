# Degree-four elliptic maps from 57 local coefficients

For any smooth plane quartic in characteristic zero, degree-four elliptic maps can be recovered by a four-parameter rank condition. At an unramified source point the unknowns specify a differential and a Weierstrass equation. A matrix with 57 rows and 20 columns has rank 17 or 18 exactly when those parameters give a degree-four map. Its kernel returns the rational functions. Three noncollinear source points suffice to cover all classes over the algebraic closure.

On the previously classified example at $a=3$, the program reconstructs cubics for **all four primitive geometric classes**, including the conjugate pair defined over $\mathbb Q(\sqrt{-2})$. It starts from the differential and target parameters for each class. It gives rank 18 for the primitive classes and rank 17 for a degree-four map that factors through an involution. It also recognizes a double cover and rejects a change of the target constant. The implementation evaluates exact rational charts and parameter tuples. Exhaustively finding the parameter tuples of a new quartic requires solving a polynomial system. The [projective method](../degree-four-incidence/README.md) gives a separate equation-only formulation. [Example input](example.json) · [Reconstruction output](receipt.json) · [Independent class count](../a3-spectrum/README.md)

## The normalized local equation

Let $K=K_C=\mathcal O_C(1)$. Work at the origin of an affine chart $F(t,s)=0$ with $F_s(0,0)\ne0$, so $s$ is a power series in $t$. Given an elliptic map unramified here, translate its value to the origin on its target. A short Weierstrass model admits a change $(X,Y)\mapsto(c^2X,c^3Y)$. This gives a unique choice for which the pulled-back invariant differential has constant coefficient one. All such source differentials and targets have the form

$$
\eta=h(t)dt=\frac{F_s(0,0)+ut+vs}{F_s(t,s)}dt,\qquad
Y^2=X^3+AX+B,\qquad 4A^3+27B^2\ne0.
$$

The span of $dt/F_s,t\,dt/F_s,s\,dt/F_s$ is the canonical linear system by plane adjunction. Write $X=t^{-2}R(t)$ with $R(0)=1$. Substituting $Y=dX/(2\eta)$ gives

$$
(tR'-2R)^2=4h^2(R^3+A t^4R+B t^6).
$$

The coefficient of $R_n$ in degree $n\ge1$ of the left side minus the right side is $-4(n+1)$. Each Taylor coefficient is therefore determined recursively, with denominators independent of $u,v,A,B$. This gives polynomial entries in four parameters over the field of the chart.

## The rank criterion and degree bound

Use the ten monomials of degree at most three in $(t,s)$ as the dehomogenized basis of $H^0(3K)$. For two cubics $N,D$, take the coefficients through $t^{56}$ of

$$t^2N-RD.$$

This is the stated linear map on a space of dimension 20. Suppose its kernel contains a nonzero pair. A nonzero cubic section has at most 12 zeros counted with multiplicity. This proves $D\ne0$ and $2\le d_0=\operatorname{ord}_P D\le12$. Locally write $N=DX+\epsilon$, where $X$ is the formal expansion and $\operatorname{ord}_P\epsilon\ge55$. The following is a global section of $(3K)^{\otimes4}\otimes K^{\otimes2}$ of degree 56, using $3K$ as additive notation for the third power of $K$:

$$
(D\,dN-N\,dD)^2-4D(N^3+A N D^2+B D^3)\eta^2.
$$

The Wronskian in this expression is independent of the local frame with the indicated tensor weights. Replacing $N$ by $DX$ makes the expression vanish formally. The leading Wronskian has order $2d_0-3$. Its error has order at least $d_0+54$. On the polynomial side the error has order at least $d_0+2(d_0-2)+55$. Thus both errors have order at least $3d_0+51\ge57$. A section of degree 56 with such a zero vanishes identically. It follows that $x=N/D$ and $y=dx/(2\eta)$ give a nonconstant morphism to the smooth target. Its $x$-coordinate has a pole of order two at the origin of the chart.

Every kernel pair represents this same $x$. The cross product of any two pairs is a section of $6K$ of degree 24, and the local error estimate gives vanishing of order greater than 24. Conversely, for the resulting map $f$, the normalized formal expansion agrees with $x$ by uniqueness in the recurrence. Therefore the kernel is precisely the space of possible denominators $H^0(3K-(x)_\infty)$.

If $f$ has degree $d$, the pole divisor of $x$ has degree $2d$. Riemann–Roch gives nullity 6 at $d=2$, nullity 4 at $d=3$, and nullity 2 or 3 at $d=4$. For $d\ge5$ the denominator bundle has degree at most two, with at most one section because a smooth plane quartic is nonhyperelliptic. Degree one is impossible. This establishes both directions of the rank criterion. [Riemann–Roch, Stacks Project](https://stacks.math.columbia.edu/tag/0B5B)

## Exhaustion and primitivity

Two distinct points on a transverse line and a third point outside that line can be chosen algebraically on any smooth quartic. At least one lies outside the zero line of each nonzero canonical differential. On its chart the map is unramified and the preceding normalization applies. To specify rank 17 or 18 algebraically, impose vanishing of the size-19 minors and require at least one nonzero size-17 minor. The target discriminant must also be nonzero. Exact elimination over the algebraic numbers thus defines an all-input procedure. There is a smaller equivalent matrix. The ten columns for $N$ have rank ten independently of the parameters, because a cubic has zero order at most 12. Eliminating them leaves 47 rows and ten columns for $D$, with allowed ranks seven and eight.

A degree-four map is imprimitive exactly when its differential vanishes along one of the involution axes of the canonical quartic. Indeed an intermediate curve would give a degree-two quotient by an involution. Each such involution fixes a line, meets the quartic and has elliptic quotient. Its invariant differential vanishes at the four intersections with that axis. Conversely, if the differential belongs to this invariant line, the difference between the map and its composition with the involution has zero differential and is constant in characteristic zero. Evaluating at a fixed point shows that constant is the identity. The map descends to the double quotient. The [involution atlas](../involution-atlas/README.md) computes these axes, including axes defined over extensions.

Primitive maps with the same projective differential line induce quotients of the Jacobian with the same connected kernel. Over the complex numbers that kernel is the connected abelian subvariety with the specified tangent hyperplane. The primitive quotient is unique up to target isomorphism and translation, so comparing lines also removes duplicates from charts. Finiteness follows independently from the Rosati form on the endomorphism lattice. For an optimal map of degree four the integral symmetric endomorphism $f^*f_*$ has eigenvalues $4,4,0,0,0,0$ on first homology, and therefore fixed squared Rosati trace norm 32. There are finitely many such lattice points. Imprimitive degree-four maps arise from the finitely many involutions and order-two isogeny kernels.

## The four reconstructed classes

In the [notation of the example](../a3-spectrum/README.md), the four classes are $P+Q,P-Q,P+S,P-S$. The respective differential lines are $U+V,U-V,U-\sqrt{-2}W,U+\sqrt{-2}W$. For the third generator take $X=2(r-2)/(r-1)$ and $Y=\sqrt{-2}\,y_R/(r-1)^2$ where $y_R=r(U/W)(V/W)$. The quotient of its invariant differential by that of $P$ is $-\sqrt{-2}W/U$. This gives the two conjugate directions.

For the rational classes take $t=W/V$ and $s=U/V-1$ or $s=U/V+1$, using the point where the respective differential is nonzero. For the other two use $t=\sqrt{-2}\,W/V$, $s=U/V-1$. The transformed quartic still has rational coefficients since it is even in $W$. The last chart allows rational arithmetic in the verifier while exhibiting the field of the transformation back to the original plane. Translation takes each chosen point to the elliptic origin. The verifier computes the five involution axes from the original equation and checks the four primitive differential lines against all of them. The program recovers $x=N/D$ from a kernel vector and independently reduces the global Weierstrass identity modulo the original chart. The formula $y=dx/(2\eta)$ recovers its other coordinate.

From this folder with SymPy installed, run

```sh
python reconstruct.py example.json
python verify.py
```

In an input file, `chart` lists $[i,j,c]$ for the term $c\,t^i s^j$. The chart must have smooth projective closure, pass through the origin and satisfy $F_s(0,0)\ne0$. The entries in `parameters` are $u,v,A,B$, given as integers or exact rational strings. A negative test here concerns that normalized tuple. Finding every tuple on a fresh positive quartic by unrestricted symbolic elimination remains a computational challenge.
