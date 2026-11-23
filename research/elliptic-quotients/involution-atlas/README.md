# Every degree-two elliptic quotient of a smooth plane quartic

A smooth rational plane quartic determines a finite list of involutions over the algebraic closure. The program computes their exact fields, centers, axes and quotient equations directly from the coefficients. Every degree-two map to a genus-one curve occurs, up to isomorphism of its target. A rational involution yields a genus-one quotient over the rationals. A rational point on the target supplies an elliptic origin. The elimination has no center-height cutoff.

On the twelve previously selected Sutherland rows, ten have zero geometric degree-two quotients. The other two have exactly one, in each case defined over the rationals and with a rational point on its target. These numbers concern the twelve inputs selected before the original atlas experiments. [Original sample](../../../applications/atlas/example.json) · [Sutherland's database](https://math.mit.edu/~drew/genus3curves.html)

| Quartic | Geometric degree-two classes | Rational involutions | Quotient information |
|---|---:|---:|---|
| Row 36518 | 1 | 1 | $j=0$ |
| Row 23640 | 1 | 1 | $j=13997521/225$ |
| Other ten blind rows | 0 each | 0 each | All geometric double covers to genus one excluded |
| Fermat control | 15 | 9 | Twelve hyperflex centers removed |
| Klein control | 21 | 0 | Galois orbit sizes $3,6,6,6$ |
| Primitive degree-four example, $a=3$ | 5 | 3 | Three rational and two conjugate over $\mathbb Q(\sqrt{-2})$ |

The two classical controls recover known geometric totals. The degree-three atlas answers a different question on these same inputs. In particular the Klein control has a rational degree-three quotient. [Bielliptic plane quartics, Badr and Bars](https://mat.uab.cat/~francesc/smoothbiellip.pdf) · [Outputs for every row](receipt.json)

## Characterization by the first polar

Work in characteristic zero. Write $F$ for the homogeneous quartic, $P$ for a nonzero column of three possible center coordinates, and $L=\nabla F(P)$ for its polar covector. The canonical model is the given plane quartic, so each automorphism acts by a projectivity. Each nonidentity involution has an isolated center and a fixed line. In coordinates $[u:v:t]$ along the axis and center, its action is $t\mapsto-t$. For the normalized sign matrix the preserving scalar on the quartic must be $+1$, since an odd quartic would contain the fixed line. The center lies off the curve because an even quartic through the center would be singular there. The quartic is therefore

$$
f t^4+A(u,v)t^2+B(u,v),\qquad f=F(P)\ne0 .
$$

Its restriction $B$ to the axis has four simple roots: a multiple root would be a singular point of the plane curve there. Riemann–Hurwitz gives genus one for the quotient. Conversely, every separable degree-two map is the quotient by its unique deck involution. The classical geometry here goes back well before this computation. [Badr and Bars](https://mat.uab.cat/~francesc/smoothbiellip.pdf)

Taylor expansion of $F(X+tP)$ gives the cubic-in-$t$ coefficient $\nabla F(P)\cdot X$. Euler’s identity gives $L\cdot P=4F(P)$. Thus the axis must be $L\cdot X=0$. On this line, evenness is equivalent to the vanishing of the linear-in-$t$ coefficient $D_P F(X)=P\cdot\nabla F(X)$. This reduces discovery to two projective unknowns and the polar criterion

$$
F(P)\ne0,\qquad (D_P F)(L\times H)=0\quad\text{identically in }H .
$$

The cross-product matrix at nonzero $L$ has image the full axis, even for an isotropic $L$. The last expression is a homogeneous cubic in $H$. Ten evaluations at the nonnegative integer triples summing to three determine it, by interpolation on the affine plane of sum three and homogeneity. Each is a homogeneous degree-ten polynomial in $P$. The recovered involution matrix is $2F(P)I-P L^t$, with gradients written as columns. Each output is checked by direct quartic substitution in its number field.

## Exact solution and removal of extraneous points

Use three disjoint center charts $[a:b:1]$, $[a:1:0]$ and $[1:0:0]$. On the source curve itself, the additional polar solutions are exactly hyperflexes. In fact, restrict $F(P+sV)$ to the tangent, with coefficients $c_i$ in powers of $s$. The constant and linear coefficients vanish. There $D_PF(P+sV)=2c_2s^2+c_3s^3$, so the polar condition says that the tangent cuts $4P$. Hyperflexes are finite in characteristic zero (the Hessian cuts the flex divisor). The remaining centers are finite by finiteness of the automorphism group of a curve of genus at least two.

On the two-variable chart, a Gröbner basis and univariate relations for $a,b$ give an exact finite algebra. Adding the squarefree parts of both relations preserves every geometric point and yields a reduced algebra, since the resulting ring is a quotient of the tensor product of two finite separable algebras. Find a separating coordinate $T=a+kb$. For a finite set with at most $N$ geometric points, at most $N(N-1)/2$ values of $k$ can identify any pair. Lexicographic elimination then gives a squarefree polynomial in $T$ and a polynomial for $b$ in $T$. Each irreducible factor defines one Galois orbit; the separating coordinate generates its residue field. The one-variable chart uses a univariate gcd. Discard the orbits satisfying $F(P)=0$. The count is exact over the algebraic closure, including for inputs with nonrational involutions. No uniform practical runtime is asserted.

For each surviving center, choose two basis columns on its axis and use $P$ as third column. In this frame the quotient is $w^2=A^2-4fB$ with $w=2ft^2+A$. The program verifies the homogeneous identity and squarefreeness in the orbit field. Binary quartic coefficients in the receipt are in increasing powers of $u$ on $v=1$.

## Independent primitivity check on the degree-four cover

For the [$a=3$ construction](../primitive-degree-four/README.md), the pullback of the target differential has canonical line $U+V=0$. For an involution, the invariant differential line is its fixed axis in the canonical plane. Indeed $M=2F(P)I-PL^t$ has eigenvalues $d,d,-d$, $d=2F(P)$, and $F(MX)=d^4F(X)$. In the plane-residue presentation the action on linear forms is $(\det M/d^4)M^t=-M^t/d$. Its invariant subspace consists of multiples of $L$. The five axes in the receipt have forms proportional to $U$, $V$, $W$, $18TV-36W$ with $T^2=-2$ (the last expression represents two axes). None is $U+V$. Thus no degree-two deck transformation can fix this degree-four map. This gives an independent exclusion of every quadratic intermediate field.

## Replay

With Python 3.11 or newer and SymPy installed, run these commands from this directory.

```sh
python classify.py ../row-23640/quartic.json
python verify.py
```

The first command accepts any smooth rational ternary quartic in the sparse monomial format used by the other discovery programs. The second recomputes the full twelve-row atlas and three controls and compares every orbit and quotient equation with the stored receipt.
