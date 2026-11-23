# From finite local map data to every rational hyperelliptic morphism

Given two smooth hyperelliptic equations over the rationals, with target genus at least two, the new compiler and reconstruction engine determine **every nonconstant rational morphism between their smooth projective curves**. The degree bound comes from the genera. The coefficient height bound comes from the exact finite incidence algebra. The algorithm accepts the two equations and requires neither a proposed map nor a coefficient cutoff. Its general elimination step can be expensive. The executed examples and the complete theorem have the respective scopes recorded below. [General compiler](rational_maps.py) · [Exact finite-algebra solver](finite.py)

For an elliptic target, a second theorem determines every map whose fibre over the selected point at infinity is the single source point at infinity. In genus two this specializes to a closed coefficient classification, reducing eight coefficient unknowns to two equations in one scale parameter. In higher genus a complete polynomial-decomposition reduction handles constant vertical multipliers. [Polynomial maps](maps.py) · [Genus-two classification](cubic-proof.md) · [Constant-multiplier reduction](composition.py)

The local data in these theorems are points of an explicitly constructed **map incidence scheme**. Prime-field counts of the source and target supply necessary restrictions on the existence of a map. Recovering a map uses the equations of its graph and exact identities. The reconstruction implements classical rational univariate representation and Hensel lifting. Its contribution here is a general interface, proved map-degree bounds, complete map templates, the small cubic classification and replayable certificates. [Rouillier, 1999](https://doi.org/10.1007/s002000050114) · [Demin, Rouillier and Ruiz, 2024](https://arxiv.org/abs/2402.07141)

## Every map to a hyperelliptic curve of genus at least two

Let

$$C:y^2=F(x),\qquad D:v^2=G(u),$$

where the rational polynomials are squarefree, of degrees $n,m\ge5$. Either degree may be odd or even. The genera are $g=\lfloor(n-1)/2\rfloor$ and $h=\lfloor(m-1)/2\rfloor$. A nonconstant map has degree

$$1\le d\le \left\lfloor\frac{g-1}{h-1}\right\rfloor.$$

This is the Riemann–Hurwitz bound. For every such $d$, use coprime polynomials $A,B$ with $\max(\deg A,\deg B)=d$, and put

$$u=A/B,\qquad v=y\,C_0/B^{h+1}.$$

Write $G^{\rm hom}(A,B)=\sum_{i=0}^m g_iA^iB^{m-i}$ and $\varepsilon=2h+2-m$, which is zero or one. The exact coefficient equation is

$$B^\varepsilon G^{\rm hom}(A,B)=F C_0^2,\qquad
\deg C_0\le\left\lfloor\frac{(2h+2)d-n}{2}\right\rfloor.$$

There are two disjoint affine charts. On the first, $A$ is monic of degree $d$. On the second, $\deg A<d$ and $B$ is monic of degree $d$. An additional variable $z$ and equation $z\operatorname{Res}(A,B)=1$ remove every common-factor solution. Scaling $(A,B)$ by $\lambda$ scales $C_0$ by $\lambda^{h+1}$, so each rational map has exactly one normalized coefficient point. The target is fixed, and opposite vertical signs count as different maps.

**Exhaustion proof.** The target regular differentials $du/v$ and $u\,du/v$ pull back to anti-invariant regular differentials on the hyperelliptic source. Their ratio proves that $u\circ f\in\mathbb Q(x)$. Taking the quotient of its differential by the pullback of $du/v$ then proves that $v\circ f$ is anti-invariant, hence is $y$ times an element of $\mathbb Q(x)$. Write its horizontal coordinate as a coprime $A/B$. At a zero of $B$ of order $k$, the homogeneous numerator is nonzero because $A$ is nonzero there. The identity for the square of $v$ shows that the order of $C_0$ is $(\varepsilon k-\operatorname{ord}F)/2$, an integer at least zero. Away from $B=0$, squarefreeness of $F$ likewise rules out poles of $C_0$. Thus $C_0$ is polynomial, with the displayed degree bound. Conversely, any coefficient solution gives a nonconstant map of the function fields, which extends to the smooth projective curves. The field degrees show that its degree is $d$.

The de Franchis finiteness theorem makes each coefficient incidence finite. This is also a standard effective construction through a bounded morphism scheme. Every finite-dimensional quotient algebra is treated by the solver, including nonreduced presentations. Therefore its Gröbner and local-reconstruction stages terminate. This proves a complete decision for the stated input class. Its runtime is governed by exact elimination and coefficient growth. [Riemann–Hurwitz](https://stacks.math.columbia.edu/tag/0C1B) · [Baker, González-Jiménez, González and Poonen, computable de Franchis–Severi, §5.3](https://math.mit.edu/~poonen/papers/finiteness.pdf)

At the unramified degree $g-1=d(h-1)$ the compiler uses a further exact simplification. Put $W=A'B-AB'$. The pullbacks of all target differentials are

$$\frac{A^jB^{h-1-j}W}{C_0}\frac{dx}{y},\qquad0\le j<h.$$

The cases $j=0,h-1$ and coprimality of $A,B$ force $W/C_0$ to be polynomial. Regularity bounds its degree by $g-1-d(h-1)=0$. It is a nonzero constant, so the implementation substitutes $C_0=sW$ and retains only the scalar $s$.

## All maps with one point above infinity, including elliptic targets

For odd-degree models with $n,m\ge3$, let the selected rational Weierstrass points be the points at infinity. Every map with a singleton fibre above target infinity is of the form

$$u=U(x),\qquad v=yV(x),\qquad G(U)=FV^2.$$

The pointed Abel maps show that such a map commutes with the hyperelliptic involutions. Regularity away from infinity then gives the polynomial form. Write $d=\deg U$, $e=\deg V$. Degree comparison gives $md=n+2e$. Differentiation and squarefreeness of $G$ imply $V\mid U'$. Therefore

$$0\le e\le d-1,\qquad d\le\left\lfloor\frac{n-2}{m-2}\right\rfloor.$$

This supplies a finite degree list even when $D$ is elliptic. For each degree, equate coefficients in $G(U)-FV^2$ and invert the leading coefficient of $U$. All geometric points of this incidence are reduced. Indeed a tangent vector satisfies $G'(U)\delta U=2FV\delta V$. The polynomials $G'(U)$ and $FV$ are coprime, so $FV$ divides $\delta U$. Since $\deg(FV)=n+e>d$, both variations vanish. The inverse variable is then fixed. This also gives an elementary proof of finiteness.

For $n=5,m=3$ the only permitted pair is $(d,e)=(3,2)$. The [complete coefficient classification](cubic-proof.md) proves the two source-coefficient obstructions and recovers all maps from two equations of degrees four and six in one scale. It includes every rational and geometric map in this infinity-fibre class, with arbitrary rational leading coefficients and translated coordinates.

When $e=0$, the source degree is $md$. Divide $F$ by its leading coefficient. Its formal $m$th root at infinity determines a unique monic degree-$d$ polynomial $H$ with zero constant term. Polynomial division either finds an exact $F/\operatorname{lc}(F)=K(H)$ or proves the absence of a constant-multiplier map. Every map must have $U=aH+b$ and $V=s$. The next-to-leading coefficient determines $b$ linearly in $a$. A univariate gcd determines the allowed nonzero $a$, and $\operatorname{lc}(F)s^2=\operatorname{lc}(G)a^m$ supplies the vertical scale. This is a specialization of classical polynomial decomposition, independently compared with SymPy's implementation in the solver study. [Kozen and Landau, 1989](https://www.cs.cornell.edu/kozen/Papers/poly.pdf) · [Comparison](../solver-comparison/README.md)

## A complete local reconstruction from a finite coefficient algebra

Let $I\subset\mathbb Q[c_1,\ldots,c_r]$ be any zero-dimensional coefficient ideal, including the inverse equations for its open conditions. Set $A=\mathbb Q[c]/I$ and $D=\dim A$. A reduced Gröbner basis provides a standard-monomial basis and multiplication matrices $M_i$.

The rank $N$ of the trace pairing $(a,b)\mapsto\operatorname{Tr}(m_{ab})$ is the number of distinct geometric points. At a geometric point $P$, the regular trace weights evaluation by the positive length $\mu_P$ of its local algebra. Try

$$\ell_k=c_1+kc_2+\cdots+k^{r-1}c_r.$$

Each pair of distinct points excludes at most $r-1$ integers. Consequently some $0\le k\le(r-1)\binom N2$ separates all points. The degree of the squarefree part of the characteristic polynomial of $M_{\ell_k}$ detects this event exactly. Let $P(T)\in\mathbb Z[T]$ be that squarefree polynomial in primitive integral normalization.

The moment matrix $H_{ij}=\operatorname{Tr}(M_\ell^{i+j})$, for $0\le i,j<N$, is invertible. Over a splitting field it is a Vandermonde matrix, the nonzero diagonal weights $\mu_P$, and the transpose Vandermonde matrix. Solving

$$\sum_{i=0}^{N-1}q_{ji}\operatorname{Tr}(M_\ell^{i+k})
 =\operatorname{Tr}(M_jM_\ell^k),\qquad0\le k<N$$

gives $c_j=q_j(T)$ on the entire reduced support. Thus a rational incidence point corresponds exactly to a rational root of $P$, and its coefficients are the evaluations $q_j(T)$.

If $a_N>0$ is the leading coefficient of primitive $P$ and $M=\max_{i<N}|a_i|$, every rational root $a/b$ in lowest terms has

$$|a|\le A_0=a_N+M,\qquad1\le b\le B_0=a_N.$$

The denominator assertion is the rational-root theorem and the numerator assertion is Cauchy's root bound. Choose a prime avoiding the leading coefficient, discriminant and coordinate-polynomial denominators. Enumerate every root of $P$ in that finite field and Hensel-lift each one to a modulus $Q>2A_0B_0$. Rational reconstruction then has at most one possible fraction on each branch. The exact equation $P(a/b)=0$ decides whether the branch gives a rational point. These operations recover every rational incidence point, so an empty result proves absence within the complete map scope.

The certificate records the original coefficient equations, the quotient basis, multiplication matrices, a primitive integral separating polynomial with positive leading coefficient, coordinate polynomials, excluded-prime conditions, lift modulus and every modular branch. Native discovery also records literal polynomial cofactors expressing every supplied basis element as a combination of the original equations. The verifier expands these identities, recomputes the supplied basis by Buchberger's algorithm, and reduces every original equation by it. These two ideal inclusions certify the exact original incidence, including its saturation equation. This route replays in Python without Singular and avoids repeating the original elimination. Certificates without cofactors instead trigger an independent exact elimination of the original equations. The verifier checks the multiplication matrices and the identities

$$\bigl(M_j-q_j(M_\ell)\bigr)^D=0.$$

It also substitutes the coordinates into every original equation and repeats the bounded local recovery. The nilpotence identities retain validity for nonreduced input ideals. Canonical integral normalization is checked before applying the rational-root bounds. This is an exact computational certificate, with polynomial arithmetic and the standard algebra algorithms forming its trusted mathematical implementation. [Implementation and verifier](finite.py)

## Executed examples and their scope

The [archive](receipt.json) contains thirty-eight independently replayable cases. Eighteen cover constant vertical multipliers across six degree/genus choices, each with a cover, a twist and a perturbation. Six exercise the complete genus-two cubic classification, six exercise the unrestricted rational-map compiler, and eight test the finite-algebra reconstruction itself. Every count in the archive is checked against its original equations. The finite-algebra controls include unequal local multiplicities, a separating form that succeeds only at the third trial, nonintegral rational coordinates, nonsplit points and a scheme that has points over every completion but no rational point.

| Executed map calculation | Geometric maps | Rational maps | Completeness scope |
|---|---:|---:|---|
| Generic odd genus-two curve to itself | 2 | 2 | Every nonconstant map to the supplied target |
| Generic even genus-two curve to itself | 2 | 2 | Every nonconstant map to the supplied target |
| $y^2=-(x^6+1)$ to itself | 24 | 8 | Every nonconstant map; the curve has no real or rational point |
| Rationally changed genus-two coordinates | 2 | 2 | Every nonconstant map; recovered horizontal coordinate $u=1/x$ |
| Generic genus-two quadratic twist | 2 | 0 | Every nonconstant map; distinguishes geometric from rational existence |
| Incompatible genus-two targets | 0 | 0 | Every nonconstant map to the supplied target |
| Generic genus-two to elliptic cubic | 2 | 2 | Every map with the singleton infinity fibre |
| Cubic with elliptic $j=1728$ | 4 | 2 | Every map with the singleton infinity fibre |
| Cubic with elliptic $j=0$ | 6 | 2 | Every map with the singleton infinity fibre |
| [Unramified genus-three to genus-two quotient](unramified-example.md) | 8 | 2 | Every nonconstant map, by an independent branch-automorphism exhaustion |

The first six rows were completed by the general coefficient compiler and finite-algebra reconstruction. The last row was completed by an independent exact search through all 336 possible images of a source branch triple and all 120 target branch triples. Its two rational maps are $u=x^2$, $v=\pm xy$. This example has two points above target infinity, so it also lies outside the singleton-infinity class. The unrestricted coefficient elimination on this degree-two example was interrupted after 120 seconds; the independent automorphism proof supplies its complete result. The theorem proves termination in principle, while the measured cost of general elimination remains an implementation issue.

The archive uses literal original-generator cofactors where compact. The empty Möbius chart instead uses independent Python elimination from its original equations, avoiding a redundant seventeen-megabyte cofactor expression. Both verification modes establish equality with the original coefficient ideal. The standalone genus-two example records use the same current equations and certificates as the archive.

The [Frobenius controls](frobenius-controls.json) additionally check two known covers and three rational exclusions, including an even-degree twist. The default replay validates these counts and the independent unramified example, and rejects eighteen invalid inputs or altered finite/map certificates plus one altered Frobenius receipt.

## Why local solvability alone is insufficient

The finite reduced equation

$$(T^2-13)(T^2-17)(T^2-221)=0$$

has a point over every completion of $\mathbb Q$ and none over $\mathbb Q$. At a prime away from $2,13,17$, either one of 13 and 17 is a square or their product is a square. At 13 use 17, at 17 use 13, and at 2 use $17\equiv1\pmod8$. A real root also exists. All three integers are nonsquares in $\mathbb Q$.

The local recovery engine rejects every branch after the derived precision threshold and exact substitution. It thereby tests the global compatibility of a branch, going beyond the existence of local points. The example is classical. [Colliot-Thélène, Arizona Winter School 2015 notes](https://swc-math.github.io/aws/2015/2015Colliot-TheleneNotes.pdf)

The general theorem above covers maps between the specified hyperelliptic curves. Recovering an unspecified target, implementing the analogous complete search on arbitrary curve presentations, and constructing algebraic correspondences between arbitrary motives remain broader questions. Unrestricted maps to elliptic targets require a different output structure: composing one nonconstant map with elliptic multiplication can already produce infinitely many maps of increasing degree. The singleton-infinity condition gives the finite class proved here. A common factor in local zeta functions alone does not construct a curve map. The optional [Frobenius filter](frobenius-proof.md) proves rational exclusions when a required target factor is absent, before coefficient elimination.

## Run and verify

The default engine uses Python 3.11+ and SymPy. The Frobenius controls also require NumPy. It also uses the established Singular solver when `Singular` is on the executable path or `CREST_SINGULAR` names its executable. Singular substantially reduces the general coefficient-elimination cost. Both routes use the same trace, reconstruction and identity checks. The archived cofactor certificates replay with Python and SymPy alone. For a certificate without cofactors, independent elimination uses Buchberger instead of F5B in the Python route, and `std` instead of `slimgb` in the native route.

From this directory run

```sh
python maps.py examples/cubic.json --output cubic-result.json
python maps.py examples/cubic.json --verify cubic-result.json
python rational_maps.py examples/mobius.json --output all-maps.json
python rational_maps.py examples/mobius.json --verify all-maps.json
python rational_maps.py examples/negative.json --screen
python verify.py
python verify.py --rebuild
```

Input arrays list coefficients in increasing powers. A general map record gives $u=A(x)/B(x)$ and $v=yC(x)/B(x)^{h+1}$. A polynomial-map record gives $u=U(x)$ and $v=yV(x)$. The `--degree` and `--chart` switches explicitly restrict the resulting claim. The unrestricted command searches every proved degree and every required chart. A resource interruption leaves a computation unfinished and supplies no absence conclusion.
