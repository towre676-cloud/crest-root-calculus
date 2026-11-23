# Repeated cells preserve the normalized spectral curve

The exact periodic-chain calculation has a direct connection with algebraic
curves. A cell repeated many times produces a polynomial of increasing
degree, but its associated Bloch spectral curve acquires square factors
rather than new geometric genus. Removing those factors recovers the same
smooth curve. This is a classical Chebyshev and Floquet identity, supplied
here with a proof and executable certificates; it is not claimed as a new
map-discovery theorem.

## Definitions and the identity

Use the cell data and matrices from [proof.md](proof.md). The cell has length
$m\geq1$, integer diagonal coefficients $a_j$, and nonzero integer bonds
$b_j$. Define

$$
B=\prod_{j=0}^{m-1}b_j\ne0,\qquad
D(x)=\operatorname{tr}\prod_{j=m-1}^{0}
\begin{pmatrix}x-a_j&-b_{j-1}\\ b_j&0\end{pmatrix}.
$$

Here $b_{-1}=b_{m-1}$ and the product order is displayed explicitly.
$D$ is monic of degree $m$. Define polynomials in two formal variables
$d,b$ by

$$
R_0=2,\qquad R_1=d,\qquad
R_n=dR_{n-1}-b^2R_{n-2},
$$

and

$$
S_{-1}=0,\qquad S_0=1,\qquad
S_k=dS_{k-1}-b^2S_{k-2}\quad(k\geq1).
$$

Thus $S_1=d$, $S_2=d^2-b^2$, and
$S_3=d^3-2b^2d$. All expressions belong to $\mathbb Z[d,b]$;
the implementation needs no division or floating-point approximation.

**Theorem.** For every positive integer $N$,

$$
\boxed{R_N(d,b)^2-4b^{2N}
       =S_{N-1}(d,b)^2(d^2-4b^2).} \tag{1}
$$

**Proof.** Work first in the splitting field of
$t^2-dt+b^2$ over $\mathbb Q(d,b)$, with roots $\lambda,\mu$.
The recurrences and their initial conditions give

$$
R_N=\lambda^N+\mu^N,\qquad
S_{N-1}=\frac{\lambda^N-\mu^N}{\lambda-\mu}.
$$

Since $\lambda\mu=b^2$ and
$(\lambda-\mu)^2=d^2-4b^2$, subtracting
$4\lambda^N\mu^N$ from $R_N^2$ proves (1). Both sides are
polynomials in $d,b$, so the equality holds in $\mathbb Z[d,b]$
and survives every specialization, including a double-root
specialization where the intermediate quotient would be undefined.
$\square$

Equivalently,

$$
R_N(d,b)=2b^NT_N\!\left(\frac d{2b}\right),\qquad
S_{N-1}(d,b)=b^{N-1}U_{N-1}\!\left(\frac d{2b}\right).
$$

The standard definitions of $T_N$ and $U_{N-1}$ by powers of a parameter
and its reciprocal appear in the
[NIST DLMF, equations 18.5.1–18.5.2](https://dlmf.nist.gov/18.5).
Formula (1) is the corresponding polynomial identity with denominators
cleared. The proof above establishes it for all $N$; the accompanying
finite symbolic checks verify the implementation, not the universal
quantifier by enumeration.

## Normalization, rather than an increase of genus

For the original cell, put

$$
F_1(x)=D(x)^2-4B^2.
$$

For the same infinite chain described by a supercell of $N$ copies, put

$$
F_N(x)=R_N(D(x),B)^2-4B^{2N},\qquad
Q_N(x)=S_{N-1}(D(x),B).
$$

The larger plane equation is
$\mathcal C_N:Y^2=F_N(x)$, while the original equation is
$\mathcal C_1:y^2=F_1(x)$. Equation (1) gives

$$
F_N=Q_N^2F_1. \tag{2}
$$

Over characteristic zero, $F_1$ is not a square even over the algebraic
closure: the coprime polynomials $D-2B$ and $D+2B$ could have a square
product only if each were a square up to a nonzero scalar. Their square
roots would then satisfy a nonzero constant difference of squares,
forcing both roots to be constant, contrary to $\deg D=m\geq1$.
The curves therefore have genuine quadratic function fields.

On the dense open set where $Q_N\ne0$, the maps

$$
(x,Y)\longmapsto(x,Y/Q_N(x)),\qquad
(x,y)\longmapsto(x,Q_N(x)y)
$$

are inverse. Consequently their smooth projective normalizations are
isomorphic over $\mathbb Q$. This is a **degree-one birational
identification**. It is not an $N$-fold map between different smooth
curves, and it is not an isogeny between Jacobians of different genera.
The displayed inverse rational formula is undefined at zeros of
$Q_N$; those points are handled by normalization, not by discarding
spectral information.

There may already be closed gaps or a nonminimal cell description.
Write

$$
F_1=H^2F_{\mathrm{sf}},
$$

where $H$ is monic and $F_{\mathrm{sf}}$ is monic and squarefree.
The word squarefree here means that an irreducible factor with even
multiplicity is removed completely, and one with odd multiplicity is
retained once. Taking the radical, which retains every distinct factor,
would give the wrong curve in general. Then

$$
F_N=(Q_NH)^2F_{\mathrm{sf}},\qquad
\widetilde{\mathcal C}_N:\ v^2=F_{\mathrm{sf}}(x),\qquad
v=\frac{Y}{Q_NH}. \tag{3}
$$

The affine curve in (3) is smooth. Its map to $Y^2=F_N$ is finite and
birational: $v$ is integral, satisfying the monic equation
$v^2-F_{\mathrm{sf}}=0$, and the two function fields agree.
It is the affine normalization; the smooth projective model also
resolves the behavior at infinity.

If $\deg F_{\mathrm{sf}}=2g+2$, Riemann–Hurwitz for its double cover
of the $x$-line gives genus $g$. In particular, if $F_1$ is squarefree,
the genus is **$m-1$, independently of $N$**. The raw degree $2mN$
of $F_N$ cannot be inserted into the squarefree genus formula:
when $N>1$, the factor $Q_N^2$ records singularities. Gap closing can
reduce the genus below $m-1$ even before the cell is repeated.

For real Jacobi coefficients, periodic Floquet theory relates the
bands to $|D(x)/B|\leq2$. Coincident band edges produce the square
factors associated with closed gaps. This setting and its finite-gap
spectral curves are classical; see Gerald Teschl,
[*Jacobi Operators and Completely Integrable Nonlinear Lattices*](https://www.mat.univie.ac.at/~gerald/ftp/book-jac/jacop.pdf),
Chapter 7 and the discussion of finite-gap operators.

## Exact examples

For the three-site spring cell used in the main panel, take
stiffnesses $(1,4,2)$ and pinning constants $(0,2,1)$, with unit
masses. Then

$$
a=(3,7,7),\qquad b=(-1,-4,-2),\qquad B=-8,
$$

and direct matrix multiplication gives

$$
D=x^3-17x^2+70x-64.
$$

Its spectral polynomial is

$$
F_1=x^6-34x^5+429x^4-2508x^3
       +7076x^2-8960x+3840.
$$

The exact gcd $\gcd(F_1,F_1')=1$, so the smooth curve has genus two.
Five copies give a degree-30 radicand but the same genus-two
normalization. More explicitly,

$$
R_5=D^5-5B^2D^3+5B^4D,
\qquad Q_5=D^4-3B^2D^2+B^4,
$$

and $F_5=Q_5^2F_1$ holds coefficient by coefficient over
$\mathbb Z$.

An important negative control is a uniform chain represented by
the nonminimal three-site cell $a=(2,2,2)$, $b=(-1,-1,-1)$.
Its discriminant is $D=x^3-6x^2+9x-2$ and

$$
D^2-4=x(x-4)(x-1)^2(x-3)^2.
$$

The normalized curve is $v^2=x(x-4)$, of genus zero. Counting the
six-degree presentation as genus two would incorrectly retain two
closed gaps. Repeating this cell still gives the same genus-zero
normalization.

The companion checker supplies five exact examples:

| Cell | Repetitions | Raw radicand degree | Squarefree radicand degree | Normalization genus |
|---|---:|---:|---:|---:|
| Two-site spring cell | 4 | 16 | 4 | 1 |
| Three-site spring cell | 5 | 30 | 6 | 2 |
| Four-site spring cell | 3 | 24 | 8 | 3 |
| Three-site signed-bond control | 3 | 18 | 6 | 2 |
| Uniform nonminimal three-site cell | 4 | 24 | 2 | 0 |

All coefficients, square multipliers, detected periods, and exact
Sturm counts of the real branch points are printed by:

```sh
python research/periodic-spectra/spectral_curve.py
```

It checks (1) symbolically for $N=1,\ldots,12$, independently forms
each cell discriminant by symbolic matrix multiplication, compares
that result with the solver, verifies both square-factor
decompositions, and checks the listed genera. These calculations
are part of verification and are excluded from timed solver work.

## Connection to the research program

The finite ring benchmark returns
$P_+(x)=R_N(D,B)-2B^N$, including every coefficient and eigenvalue
multiplicity. The curve radicand is instead

$$
F_N=P_+P_-,\qquad P_-(x)=R_N(D,B)+2B^N,
$$

where $P_-$ is the antiperiodic boundary polynomial. The curve
certificate therefore complements the finite spectrum; it does not
substitute a smaller output for the benchmark's required polynomial.

After normalization, the squarefree equations above are explicit
inputs suitable for the project's
[hyperelliptic correspondence tools](https://github.com/towre676-cloud/crest-root-calculus/blob/main/research/bounded-correspondences/README.md).
The example gives an exact bridge from a periodic quadratic energy
model to the existing curve language. It does not demonstrate a new
nontrivial map between inequivalent materials, an automatic discovery
of additional symmetry, or a new source of spectral speedup. Those
would require separate mathematical and computational evidence.
