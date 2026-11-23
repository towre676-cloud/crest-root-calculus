# Point-free Riemann–Roch spaces and symmetric powers

Given projective equations for a smooth geometrically integral curve over $\mathbb Q$ and finite homogeneous ideals for a rational divisor $D=D_+-D_-$, the [embedded-curve algorithm](embedded_curve.py) computes a basis of the complete space
$$L(D)=\{f\in\mathbb Q(C):\operatorname{div}(f)+D\ge0\}\cup\{0\}.$$
The same input equations give a finite affine atlas of every specified symmetric power. There is a faster [plane-curve implementation](riemann_roch.py) and an [exhaustive search by chord height](symmetric_power.py) for the rational symmetric square of a smooth plane curve. Divisors can have repeated points, closed points of higher degree, and common support in their positive and negative presentations. Every computation runs over the rationals. No rational base point or target library is required.

Two applications use these spaces to reconstruct maps. A bitangent and its first normal derivative identify all rationally defined, geometrically cyclic degree-three genus-one quotients of any smooth rational plane quartic. A family with nonsplit deck group has rational elliptic-map degrees exactly $2n^2,6n^2,12n^2$ and geometric degrees exactly the multiples of two or three. [Cubic recognition and reconstruction](cyclic-quotients.md) · [Full family spectrum and proof](family-spectrum.md)

## Using the programs

From this folder, with Python 3.11 or later with SymPy and NumPy installed, run

```sh
python riemann_roch.py fermat.json --duality
python embedded_curve.py space-quartic.json
python symmetric_power.py fermat.json --chord-height 2
python symmetric_power.py fermat.json --atlas-degree 3 --frame-index 2 --slope 1
python embedded_curve.py space-quartic.json --atlas-degree 2
python embedded_curve.py space-quartic.json --search-degree 2 --height 1
python cyclic_quotients.py nonsplit.json
python fermat_divisors.py 6
python family_spectrum.py 12
python verify.py
```

A plane equation uses the repository's `terms` rows $[i,j,k,c]$ for $cX^iY^jZ^k$. An embedded equation uses `coordinates` for the number of homogeneous coordinates and `equations` for the list of forms. Each term has one nonnegative exponent per coordinate and a rational coefficient. Rational coefficients may be strings such as `"-2/3"`. The source for the general embedding calculation is a smooth geometrically integral projective curve. It checks dimension, arithmetic genus, the projective Jacobian criterion and geometric connectedness. A separate source-integrality certificate is unnecessary.

Each entry of `projective_positive` or `projective_negative` has `generators`, a list of homogeneous forms defining a finite scheme on the curve, and optional positive integer `multiplicity`. Adding entries adds divisors. Omitting a side represents the zero divisor. Plane inputs can alternatively specify affine ideals in `positive` and `negative`, whose term rows are $[i,j,c]$ in $u,v$. The optional invertible rational `frame` gives the coordinate forms $(U,V,W)$ as its rows, and $u=U/W$, $v=V/W$. For projective data, the plane program chooses a chart avoiding the entire supplied support and rejects a specified chart that loses any of it. Results give a common homogeneous denominator and independent homogeneous numerators. Every linear combination over $\mathbb Q$ is a member of $L(D)$. For effective divisors $D_1,D_2$ of equal degree, $L(D_1-D_2)$ supplies a rational equivalence function precisely when the two divisors are linearly equivalent.

## Complete spaces in an arbitrary projective embedding

Write $C\subset\mathbb P^r$, let $I$ be its homogeneous ideal and let $I_+,I_-$ denote the finite divisor ideal sheaves. Choose a homogeneous section $H$ of degree $n$ vanishing on $D_+$ with the prescribed multiplicities. A product of nonzero forms, one from each effective divisor ideal and raised to its multiplicity, suffices. Multiplying by a coordinate nonzero on the curve allows $n$ to increase. At any point, the smooth local ring is a discrete valuation ring. For a numerator $N$ of degree $n$, this gives the exact local criterion
$$ N/H\in L(D)\quad\Longleftrightarrow\quad N I_+\subseteq H I_- .$$
Indeed it says $v(N)+D_+(P)\ge v(H)+D_-(P)$ at each point. Conversely, if $f\in L(D)$, then $Hf$ is a global section of $\mathcal O_C(n)$. Thus completeness requires all sections of this twist to be represented in the numerator calculation.

The program obtains the Hilbert series of the given homogeneous ideal from its initial monomials. Inclusion–exclusion on their least common multiples gives $Q(s)/(1-s)^2$. The Hilbert polynomial is $\delta n+1-g$, with $\delta=Q(1)$ and $g=1-\delta+Q'(1)$. Finite projective ideals are handled by the same calculation with a simple pole. Its numerator at one gives their length. On the curve, choose $n\delta>2g-2$, so Riemann–Roch and duality give $h^0(\mathcal O_C(n))=n\delta+1-g$. Reduce the homogeneous standard monomials of degree $n$ on every standard affine chart and take the rank of their joint coefficient matrix. This rank measures the actual space of ambient restrictions. Increase $n$ until the rank equals the known dimension. Eventual surjectivity follows from Serre vanishing applied to the ideal sheaf. In particular this algorithm also works in embeddings that omit sections in small degree. [Cohomology of projective schemes](https://stacks.math.columbia.edu/tag/01YR) · [Cohomology on projective space](https://stacks.math.columbia.edu/tag/089X)

The source check certifies geometric connectedness directly for complete intersections of positive dimension. For other presentations let $B$ be the total degree of the least common multiple of the leading Gröbner monomials. The Taylor resolution gives a regularity bound $B$ for the initial quotient, hence also for the original graded quotient by Gröbner degeneration. Thus every global section of $\mathcal O_C(B)$ lifts to a form. Choose a rational linear form $L$ whose hyperplane slice has length exactly $\delta$. The Hilbert polynomial exact sequence for multiplication by $L$ shows its sheaf kernel has zero Hilbert polynomial, hence vanishes. A generic linear form is suitable, and testing the moment-curve forms $L_i$ eventually avoids every associated point. The kernel calculation with denominator $L^B$ and divisor zero therefore computes $H^0(C,\mathcal O_C)$ before using connectedness. Its dimension must be one, which persists over the algebraic closure and excludes a disjoint or isolated component. Every closed point of the connected one-dimensional scheme has local dimension one. The verified tangent dimension at most one then makes every such local ring regular. Smoothness and geometric connectedness imply geometric integrality. [Bayer, Peeva and Sturmfels on monomial resolutions](https://arxiv.org/abs/alg-geom/9610012)

In each standard affine chart the algorithm forms the positive and negative divisor ideals by multiplication. For each numerator-basis monomial and each generator of $I_+$ it reduces their product modulo $(I_C,HI_-)$. Coefficient comparison on all charts gives one rational matrix. Its kernel is $L(D)$ by the local criterion. Checking all charts includes every pole introduced by increasing the degree of $H$.

The [space quartic](space-quartic.json) is the image $[s^4:s^3t:st^3:t^4]$. Its four linear coordinates omit $s^2t^2$. The numerator-rank test detects the missing section, raises the degree to two, and recovers $1,t/s,(t/s)^2,(t/s)^3,(t/s)^4$ for the hyperplane divisor cut out by $s^4=0$. The independent parameter substitution checks all five. A second test uses a smooth genus-one intersection of quadrics in $\mathbb P^3$ with no real point.

## Plane optimization

For a smooth plane curve of degree $d$, every section of $\mathcal O_C(n)$ for $n\ge0$ lifts to a form by the exact sequence for the plane hypersurface and $H^1(\mathbb P^2,\mathcal O(n-d))=0$. Choose affine coordinates $t=u+sv,v$ so that $[0:1:0]$ in these coordinates misses $C$. At most $d$ rational integer choices of $s$ are forbidden. The affine equation has degree $d$ in $v$ with constant leading coefficient. The monomials $t^iv^j$ with $i+j\le n$, $0\le j<d$ give the independent ambient sections in that degree.

For $D_+$ in this affine chart, the nonzero monic $t$-eliminant of its ideal bounds a terminating search for a smaller denominator. Its homogenization avoids all boundary points. At each degree take the kernel of evaluation modulo the divisor ideal. The leading forms of the candidate denominators must collectively avoid a common zero on the curve at infinity. A univariate gcd tests this. If $h_0,\ldots,h_{k-1}$ form the kernel basis, each of the at most $d$ boundary points forbids at most $k-1$ choices of the integer $c$ in $\sum c^j h_j$. A suitable denominator therefore occurs among the first $d(k-1)+1$ choices. All zeros of that denominator and all prescribed divisor points are in one affine chart, where the same ideal inclusion computes the space. The optional duality check computes $L(K-D)$ independently using $K=(d-3)H_{\rm line}$.

## Symmetric powers and rational divisor searches

The symmetric power $C^{(m)}$ of a smooth curve parametrizes length-$m$ effective divisors including nonreduced ones. A divisor over an infinite field lies in at least one chart $L_i\ne0$, where $L_i=x_0+i x_1+\cdots+i^r x_r$, $0\le i\le rm$. Each of the at most $m$ geometric support points rules out at most $r$ choices. Put $z_j=x_j/L_i$. In the chart choose $t=z_1+s z_2+\cdots+s^{r-1}z_r$ so that different support points have different values and $dt$ is nonzero along the tangent at every multiple point. Each pair or tangent condition rules out at most $r-1$ choices. With $k$ distinct geometric support points, at most $\binom{k}{2}+(m-k)\le\binom{m}{2}$ such conditions occur. Thus $0\le s\le(r-1)\binom{m}{2}$ suffices.

On this chart $t$ generates the divisor algebra. Locally at a point of multiplicity $j$, $t-t(P)$ is a uniformizer modulo its $j$th power. Chinese remaindering combines these local descriptions. There is a unique monic $h(t)$ of degree $m$ and unique polynomials $g_2(t),\ldots,g_r(t)$ of degree below $m$ for $z_2,\ldots,z_r$. Substitute $z_1=t-\sum_{j=2}^rs^{j-1}g_j(t)$ and the other coordinates in every source equation. All $m$ remainder coefficients modulo $h$ must vanish; conversely their vanishing constructs the divisor algebra $\mathbb Q[t]/h$ with that embedding. These polynomial equations are the chart itself: on a family the power basis $1,t,\ldots,t^{m-1}$ is an open basis condition, and monic division is valid over the parameter ring. The [general program](embedded_curve.py) prints the equations of any chart. The [plane program](symmetric_power.py) also encodes a supplied divisor in a separating chart.

For every degree and embedded source, `--search-degree m --height B` enumerates rational divisors whose first admissible chart has coefficient height at most $B$. The coefficients of $h$ below its leading term and of the interpolants are written with a common positive denominator in primitive integral form. Height is the largest absolute numerator or denominator. The search tests the chart equations by homogeneous integer arithmetic and uses invertibility and a power-basis determinant to choose a unique chart for every divisor. At most $(rm+1)((r-1)\binom{m}{2}+1)B(2B+1)^{rm}$ tuples are considered. Increasing $B$ eventually reaches each rational divisor.

For plane curves, any rational length-two divisor spans a unique rational line, including a doubled point and a pair of conjugate points. For each primitive integer line of height at most $B$, factor its binary section over $\mathbb Q$ and enumerate all quadratic subdivisors with allowed multiplicities. Every result is returned with its projective generators. This is exhaustive through chord height $B$. The atlas and increasing-height search leave the general arithmetic problem of deciding all rational points of a symmetric power open.

## All rational symmetric powers in a no-real-point example

For $C:X^4+Y^4+Z^4=0$, the height-one chord search recovers four degree-two contact divisors $D_0,\dots,D_3$ on $X\pm Y\pm Z=0$. In the program they are ordered by the integer line triple, beginning with $(1,-1,-1)$. Each satisfies $2D_i=\operatorname{div}_C(L_i)$ as a hyperplane section. The [Picard-group theorem of Yasuhiro Ishitsuka, Tetsushi Ito and Tatsuya Ohshita](https://arxiv.org/abs/2107.05902), Theorem 1.1 and Proposition 6.1, identifies the four rational divisor classes of degree two. Their theorem distinguishes the four degree-zero divisor classes from the eight rational points of the Picard scheme's identity component. We use the divisor classes in the following parametrization.

| Degree $d$ | All rational effective divisors, as sets of rational parameters |
|---|---|
| 0 | One empty divisor |
| Odd | Empty, since complex conjugation pairs the geometric points |
| 2 | Four individual divisors |
| 4 | One $\mathbb P^2(\mathbb Q)$ and three $\mathbb P^1(\mathbb Q)$ |
| Even $d\ge6$ | Four disjoint copies of $\mathbb P^{d-3}(\mathbb Q)$ |

Indeed for $d=2n\ge2$ the classes are $(n-1)D_0+D_i$. For $d=4$ the canonical class has three sections, and the other three classes each have two by duality. Above degree four, the Riemann–Roch dimension is $d-2$. The [parametrization](fermat_divisors.py) computes an explicit basis on each class, using $2D_i=L_i|_C$ to impose at most four residual linear conditions at any degree. A projective vector of rational coefficients produces $D+\operatorname{div}(\sum c_jN_j/H)$. The program gives the representative $D$, all $N_j$, and $H$. Scalar-equivalent vectors give the same effective divisor. The [replay](verify.py) compares these optimized bases' dimensions with a direct ideal calculation, checks all three nontrivial divisor differences and reproduces the small-degree classes and the map results. [Exact receipt](receipt.json)
