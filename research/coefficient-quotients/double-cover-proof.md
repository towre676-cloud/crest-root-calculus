# Degree-two quotients from one hyperelliptic equation

Given an expanded polynomial over $\mathbb Q$, remove square factors and keep its leading unit. Write the resulting curve as $C:y^2=F(x)$ with $F$ squarefree and $g(C)\ge2$. This search returns every rational degree-two morphism from $C$ to a positive-genus curve, up to rational isomorphism of the target. It allows even- and odd-degree equations, and any fibre over infinity.

A degree-two extension of function fields in characteristic zero has a unique nonidentity deck involution. Conversely, any rational involution defines a rational degree-two quotient by its fixed field. The hyperelliptic involution gives the genus-zero quotient. Every other one induces a nonidentity involution of the projective $x$-line because the hyperelliptic involution is unique and central. Extra involutions and their relation to elliptic subcovers are classical. The implementation here computes the involutions over the input field and writes their exact quotient formulas. [Gutiérrez and Shaska](https://doi.org/10.1112/S1461157000000917) · [Shaska, elliptic subcovers](https://www.math.uniri.hr/NATO-ASI/abstracts/shaska.pdf)

## Two projective charts

An affine nonidentity involution has the form $x\mapsto 2a-x$. The ratio $F(2a-x)/F(x)$ is the constant $(-1)^{\deg F}$ if it preserves the branch divisor. An order-two lift needs this ratio to be one. Thus the affine chart only needs even degree $n$. Comparing the coefficient of degree $n-1$ in $F(a+z)$ determines $a=-F_{n-1}/(nF_n)$. Test whether that shifted polynomial is even. If so the two lifts act on $y$ by plus and minus one.

A non-affine projective involution has a unique form
$$
 x\longmapsto M(x)=a+\frac{k}{x-a},\qquad k\ne0.
$$
Indeed its projective matrix has trace zero and can be scaled to have lower-left entry one. Put $h=g(C)+1$ and write $F(a+z)=\sum_{j=0}^{2h}p_j(a)z^j$, with $p_{2h}=0$ in odd degree. Invariance of the homogeneous branch form says that the ratio between $z^{2h}F(a+k/z)$ and $F(a+z)$ is a nonzero constant. A lift acts on $y$ by $c/(x-a)^h$ and the constant is $c^2$. Its square acts on $y$ by $c^2/k^h$. An order-two lift therefore requires
$$
 z^{2h}F(a+k/z)=k^h F(a+z),\qquad c^2=k^h.
$$

For $j<h$, matching coefficients and cancelling known powers of $k$ gives
$$
 p_{2h-j}(a)k^{h-j}=p_j(a).
$$
The middle coefficient is automatic, and coefficient pairs above the middle follow from these equations. At $k=0$ the first two equations would give $p_0(a)=p_1(a)=0$, contrary to squarefreeness. The two-variable system therefore needs no extra variable to invert $k$. It is finite since each point gives a distinct projective permutation of the branch set, and a projective transformation is determined by its action on three distinct points. Exact rational-point reconstruction supplies the complete list of rational $(a,k)$. For each pair, the rational-square condition gives precisely the two signed lifts or an empty list. This last test distinguishes a rational symmetry of the branch set from a rational involution of the curve.

## Construct the quotient and its map

Let $\sigma(x)=M(x)$ and $\sigma(y)=c(x)y$. Take $u=x+M(x)$ unless this is constant, in which case use $u=xM(x)$. The second choice is nonconstant in the exceptional case. Both are fixed and have degree at most two as rational functions of $x$. Invariance and nonconstancy force degree two, so $u$ generates the fixed subfield of $\mathbb Q(x)$.

Take $v_0=y+\sigma(y)$, or use $xy+\sigma(xy)$ when the first is zero. The second is nonzero in that case since $M$ is nonidentity. Both choices are fixed, and $v_0$ is anti-invariant under the hyperelliptic involution. Its square descends to some $N(u)/D(u)$, obtained by reducing in the quadratic algebra defined by the equation for $u(x)$. The field $\mathbb Q(u,v_0)$ has degree two over $\mathbb Q(u)$ while $\mathbb Q(x,y)$ has degree four over $\mathbb Q(u)$. Thus $\mathbb Q(u,v_0)$ is the full fixed field of the deck involution.

Write $ND=A_0^2G$ with $G$ squarefree and put $v=v_0D(u)/A_0(u)$. This gives the polynomial target $v^2=G(u)$. The program verifies the exact identity in the original source equation, including its initial square factors. Distinct deck involutions give distinct maps up to target isomorphism, since the deck involution is uniquely determined by the map. When the output has genus one, its equation specifies a genus-one curve over $\mathbb Q$. Existence of a rational point to use as an elliptic origin is a separate test.

## A fibre with two points

For $F=(x^2-1)(x^2-4)(x^2-9)$, the map $(x,y)\mapsto(-x^2,y)$ gives a degree-two elliptic quotient. Both points at infinity on the source map to the unique point at infinity on the cubic target. The search also returns the quotient from $(x,y)\mapsto(-x,-y)$ acting on the function field. After a projective coordinate change that moves a branch point to infinity, the same deck transformations are recovered in the non-affine chart. [Runnable input](double-cover-example.json) · [Independent branch-point replay](verify_double_covers.py)
