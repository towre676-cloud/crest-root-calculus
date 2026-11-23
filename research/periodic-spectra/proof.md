# Exact spectra of periodic Jacobi rings

## The task and its scope

The input consists of a cell of length $m\geq 1$, integer diagonal entries
$a_0,\ldots,a_{m-1}$, nonzero signed integer bonds
$b_0,\ldots,b_{m-1}$, and a positive integer number $N$ of repetitions. Put
$L=mN$ and require $L\geq 3$. Extend both coefficient lists periodically to
length $L$. The matrix $A\in\operatorname{Sym}_L(\mathbb Z)$ has

$$
A_{ii}=a_i,\qquad A_{i,i+1}=A_{i+1,i}=b_i
\quad(i\bmod L),
$$

with all other entries zero. Every edge is added once. The output is the
**entire monic integer polynomial** $\det(xI_L-A)$, with its $L+1$
coefficients, including zeros and algebraic multiplicities. A factorization,
a list of approximate eigenvalues, or a few traces would be a different
output contract.

The reduction below is classical transfer-matrix/Floquet theory and the
Chebyshev trace identity. The contribution tested here is a reproducible
implementation, detection of an actual repeated cell, exact certificates,
and comparison against both generic and structure-aware PARI/GP. The
identity itself is not asserted to be new.

## The cell certificate

Define the polynomial matrices

$$
C_j(x)=\begin{pmatrix}x-a_j&-b_{j-1}\\ b_j&0\end{pmatrix},
\qquad C(x)=C_{m-1}(x)\cdots C_0(x),
$$

where $b_{-1}=b_{m-1}$. Put

$$
B=\prod_{j=0}^{m-1}b_j,\qquad D(x)=\operatorname{tr}C(x).
$$

Then $D\in\mathbb Z[x]$ is monic of degree $m$ and
$\det C=B^2$. The assertion about the determinant follows by multiplication:
$\det C_j=b_jb_{j-1}$. The leading term of the trace comes from selecting
$x$ in the top-left entry of every matrix; all other contributions have
smaller degree.

Let

$$
R_0(x)=2,\qquad R_1(x)=D(x),\qquad
R_n(x)=D(x)R_{n-1}(x)-B^2R_{n-2}(x).
$$

**Theorem.** For every input above,

$$
\boxed{\det(xI_L-A)=R_N(x)-2B^N.}
$$

The product $B$ is signed. Replacing $B^N$ by $|B|^N$ is incorrect in
general. The recurrence uses $B^2$ and is valid for either sign of $B$.

### Proof with multiplicities

Write $P_{r:s}(x)$ for the determinant of the open tridiagonal path on
vertices $r,r+1,\ldots,s$. An empty path has determinant one. Its recurrence
is

$$
P_{r:s}=(x-a_s)P_{r:s-1}-b_{s-1}^2P_{r:s-2}.
$$

Expanding a cyclic determinant by permutations gives, for $L\geq 3$,

$$
\det(xI_L-A)
=P_{0:L-1}-b_{L-1}^2P_{1:L-2}-2\prod_{j=0}^{L-1}b_j. \tag{1}
$$

Indeed, terms omitting the wrap edge form the first path determinant.
Terms using that edge as a transposition contribute the second term.
The two oriented full cycles each contribute minus the product of the
bonds. These exhaust the nonzero permutation terms on a cycle; the
same formula holds for the triangle $L=3$.

Induction in the path recurrence, or multiplication of the displayed
$2\times2$ matrices, gives

$$
\operatorname{tr}(C_{L-1}\cdots C_0)
=P_{0:L-1}-b_{L-1}^2P_{1:L-2}. \tag{2}
$$

Periodicity makes the product on the left equal to $C(x)^N$, and the
product of all bonds is $B^N$. Cayley–Hamilton for $C$ yields

$$
C^2-DC+B^2I_2=0.
$$

Taking traces after multiplication by successive powers of $C$ proves
$\operatorname{tr}(C^n)=R_n$. Equations (1) and (2) establish the theorem
as an identity in $\mathbb Z[x]$. The argument proves all multiplicities;
it does not merely compare sets of eigenvalues. $\square$

This proof also supplies a separate verifier: evaluate the two open-path
continuants in (1) on the **full expanded ring**, without using the
declared or detected cell. The verifier in this directory uses that
method and compares every integer coefficient. Small cases are also
checked by direct symbolic matrix determinants.

## The root-of-unity and hyperelliptic descriptions

The recurrence is equivalently

$$
R_N=2B^N T_N\!\left(\frac{D}{2B}\right),
$$

where $T_N$ is the Chebyshev polynomial of the first kind. The rational
notation is only explanatory: the integer recurrence avoids divisions.
For a Bloch multiplier $z$, the cell's eigenvalue equation is

$$
D(x)=B(z+z^{-1}). \tag{3}
$$

Selecting $z^N=1$ reconstructs the spectrum of the finite ring. In
particular,

$$
\det(xI_L-A)
=\prod_{\zeta^N=1}\bigl(D(x)-B(\zeta+\zeta^{-1})\bigr). \tag{4}
$$

One proof of (4) uses the monodromy $C/B$, which has determinant one:
$\det(C/B-zI)=0$ is (3), and the cyclic boundary condition makes $z$
an $N$th root of unity. Alternatively, (4) follows directly from the
Chebyshev polynomial identity, so it remains valid when eigenvalues
coincide.

For $d\geq3$, let $\Psi_d(u)\in\mathbb Z[u]$ be the minimal polynomial
of $\zeta_d+\zeta_d^{-1}$ for a primitive $d$th root of unity. It has
degree $\varphi(d)/2$. Pairing inverse roots in (4) gives the integral
factorization

$$
\det(xI_L-A)
=(D-2B)(D+2B)^{\mathbf1_{2\mid N}}
\prod_{\substack{d\mid N\\d\geq3}}
\left[B^{\varphi(d)/2}\Psi_d(D/B)\right]^2. \tag{5}
$$

This identifies an exact connection with Crest–Root's existing
character decompositions: local character sectors assemble into a
global polynomial, and exact identities certify the assembly. That
connection does not make a general lattice problem a Jacobi ring.

Eliminating $z$ from (3) produces the spectral curve

$$
Y^2=D(x)^2-4B^2,\qquad Y=B(z-z^{-1}). \tag{6}
$$

If the polynomial on the right is squarefree, its smooth projective
model has genus $m-1$. If gaps close and square factors occur, one
must first remove square factors and normalize; the genus can be
smaller. Existing hyperelliptic-map machinery can therefore be
considered on appropriate smooth spectral curves. This experiment
does **not** prove that such maps improve spectral computation or
classify all periodic materials.

The [repeated-cell normalization identity](spectral-curve.md) proves
precisely how repeating a cell introduces square factors into (6),
while preserving the normalized spectral curve. Its independent
symbolic checks run with this directory's default verifier.

## Why the validation checks matter

The matrix recognizer must check a square matrix of integer entries,
symmetry, zeros away from the cycle, and nonzero cycle bonds. Its
reported period $p\mid L$ must satisfy both
$a_{i+p}=a_i$ and $b_{i+p}=b_i$ for every $i$. A repeated diagonal
alone is insufficient. A cell of length one or two is valid when
the expanded ring has at least three vertices. A ring of length one
or two is rejected: ordinary adjacency-matrix conventions for
loops and parallel edges would otherwise require a separate
definition.

A changed diagonal or bond can destroy a proposed period. The
detector must then choose a larger true period, possibly $L$ itself.
The theorem still applies with $N=1$; the computation remains an
exact cyclic-continuant calculation, but there is no repetition
benefit. A defect cannot be silently erased, averaged, or treated as
roundoff. Likewise, a chord or a nonsymmetric entry lies outside
this matrix class and must be rejected rather than projected to a
ring.

The recurrence uses exact integers throughout. This avoids
floating-point sign decisions and roundoff during the calculation.
It does not remove uncertainty in a physical model's measured
coefficients. Rational inputs could be handled after recording an
explicit common-denominator scaling, but this experiment's input
contract is integral.

## A concrete mechanical interpretation

For identical unit masses, choose positive integer spring constants
$k_i$ on the cycle and nonnegative integer pinning constants $p_i$.
Set

$$
b_i=-k_i,\qquad a_i=k_{i-1}+k_i+p_i.
$$

Then $A$ is the stiffness matrix of the harmonic energy

$$
E(u)=\frac12\sum_i k_i(u_{i+1}-u_i)^2
 +\frac12\sum_i p_i u_i^2,
\qquad 2E(u)=u^TAu.
$$

Every eigenvalue is nonnegative. When all $p_i=0$, connectedness of
the cycle makes the zero eigenspace exactly the constant vectors;
the characteristic polynomial has a simple factor $x$. If at least
one $p_i>0$, the quadratic form is positive definite: a vector of
zero energy would have to be constant and vanish at the pinned
site. This provides useful exact controls independently of the
spectral algorithm. Negative pinning values instead model unstable
quadratic directions and are allowed by the algebraic input class.

Thus the experiment is relevant to a periodic harmonic chain, a
one-dimensional periodic quadratic energy landscape, and special
symmetric reaction–diffusion linearizations. It does not supply
first-principles force constants, nonidentical mass dynamics,
three-dimensional crystal symmetries, molecular enumeration, or
pharmaceutical predictions. Those require additional models and
baselines.

## Comparison and attribution

PARI/GP's default `charpoly(A)` is a genuine optimized baseline;
for integer matrices its documented default selects modular
Hessenberg reconstruction. A structure-aware GP comparator must
also receive the same matrix class and implement the same
transfer/continuant reduction. A gain over the former establishes
a benefit of using the structure. A gain over the latter would
need to survive implementation and repeated-timing checks before
being described as an implementation advantage. Neither alone
would establish a new asymptotic theorem or superiority over
PARI/GP in general.

Primary references for the classical background and comparator:

* Gerald Teschl, *Jacobi Operators and Completely Integrable Nonlinear
  Lattices*, AMS Mathematical Surveys and Monographs 72 (2000),
  Chapter 7, periodic Jacobi operators and Floquet theory.
  [Author's book](https://www.mat.univie.ac.at/~gerald/ftp/book-jac/jacop.pdf).
* A. Aghaei, K. Dayal and R. S. Elliott, “Symmetry-adapted phonon
  analysis of nanotubes,” *Journal of the Mechanics and Physics of
  Solids* 61 (2013), 557–578.
  [DOI](https://doi.org/10.1016/j.jmps.2012.09.008).
* PARI/GP documentation, `charpoly`, linear algebra chapter.
  [Official documentation](https://pari.math.u-bordeaux.fr/dochtml/html/Linear_algebra__Vectors__matrices_and_sets.html#charpoly).

The proof above is self-contained and fixes the exact conventions
used by the implementation, including signed bonds and triangles.
