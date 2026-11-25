# Certified period searches and exact homology lattices

For nine genus-two models, this study checks the period search against the geometric endomorphism ring and the Néron–Severi lattice. The comparison identifies each candidate inside a common integral homology basis. It gives an independent check of the saturation of eight split lattices and determines the full quaternionic order in the ninth example.

The verifier starts from the equations and uses explicitly certified maps or correspondences and verified rank bounds to determine the full lattices. It then compares them with two independent bounded enumerations from the enclosed periods. Its algebraic inputs and degree limits are specified below.

## Period enclosures

The change of variables $x=z_0+1/t$ uses a Gaussian-rational $z_0$. Around one root of the moved sextic, the other five roots are joined to it by straight segments. Each loop traverses a segment and returns on the other sheet. An exact dyadic rotation sends the proposed angular cut to the negative real axis. Disjoint root balls and strict angular inequalities certify the order of the outgoing rays. For every other root along a segment, both endpoints must lie in an open half-plane through that root. Convexity gives a single square-root branch along the entire segment.

The square roots of the outgoing ray directions are taken against the same angular cut. In the local coordinate $\sqrt{t-e_1}$ at the shared branch point, the two rays of every lifted loop are antipodal and any two loop directions interleave. Comparing the regularized square roots at the common endpoint makes their signs coherent. With the resulting orientation the five loops have pairwise intersection $+1$ above the diagonal, up to an overall sign convention. The first four have unimodular intersection matrix, so they form a basis of the genus-two homology. The alternating relation of the five periods and Riemann’s first bilinear relation are checked by ball containment as consistency checks.

The substitution $t=e_1+(e_j-e_1)\sin^2(\pi s/2)$ cancels the square-root endpoint singularities. FLINT/Arb integrates on dyadic subintervals of $0\le s\le1$, with the derivative bound check using analytic square roots. From the enclosed periods, ball matrix operations give the complex structure $J$ on real homology and the positive Gram matrix $G$ associated to the principal polarization. Midpoints propose the diagram. The certified topology and quadrature determine the enclosures.

## Exact enumeration in a ball-defined metric

Let $M$ denote a real symmetric form given by a ball matrix. At an integer scale $S$, round its entries to a symmetric integral matrix $B$. Bound the residuals $|(SM-B)_{ij}|$ by symmetric integers $r_{ij}$. With $A=B-\operatorname{diag}(\sum_jr_{ij})$, the inequality $2|v_iv_j|\le v_i^2+v_j^2$ proves $A\preceq SM$. The verifier requires an exact positive rational LDL decomposition of $A$.

Given $A=LDL^{\mathsf T}$, the enumerator fixes the coordinates from last to first. For fixed later coordinates, the remaining contribution at index $i$ is $D_i(v_i+a/b)^2$ with rational shift $a/b$. If the remaining budget is $H$, the exact integer condition is $|bv_i+a|\le\lfloor\sqrt{Hb^2/D_i}\rfloor$. This yields every vector including boundary vectors, one representative per sign. The optional reduction of an ambient Gram matrix changes the search basis by a matrix with determinant $\pm1$ checked in exact arithmetic.

For elliptic subvarieties, Hermite’s rank-two bound supplies a nonzero vector $u$ of squared length at most $2d/\sqrt3$ in any saturated homology plane of polarization degree $d$. The program therefore enumerates all possible witnesses through degree $X$ using $A\preceq SG$ and the exact integer bound $\lfloor 2XS/\sqrt3\rfloor$. For a witness from a true plane the vector
$$
\frac{d(u\wedge Ju)}{u^{\mathsf T}Gu}
$$
is its primitive oriented Plücker vector. Dyadic interval endpoints and exact integer inequalities reject every degree-coordinate pair that cannot give an integer vector. The remaining vectors must satisfy primitivity, the Plücker relation and the polarization-degree identity. Ball containment gives a shortlist. Every plane in the shortlist must subsequently agree with the plane of an exact isotropic class.

For endomorphisms, the form
$$
Q(R)=\operatorname{tr}(R^{\mathsf T} G R G^{-1})+2^{20}\lVert JR-RJ\rVert_F^2
$$
is positive definite on the real matrices. Its value on a genuine endomorphism is the real homology trace of $R^\dagger R$. The same integral lower-bound method exhausts all matrices of possible trace at most $T$. For each enumerated matrix outside the independently certified rational endomorphism space, the replay requires a commutator entry whose enclosure excludes zero. For those inside the space, exact trace arithmetic applies the actual cutoff. With the default penalty the ambient norm excludes all extraneous matrices on these cases. A low-penalty control on $C$ admits 181 additional matrices to the shortlist and rejects each by an enclosed commutator entry.

At fixed cutoffs, a sequence of enclosures shrinking to the same period matrix in a fixed basis eventually excludes every false candidate. To see this, keep one positive definite integral lower form throughout refinement. Its ellipsoid supplies a fixed finite set of vectors. For elliptic candidates, each real Plücker coordinate vector has finitely many integer targets in the initial enclosure. A false target differs from the real vector in at least one coordinate, so a sufficiently narrow ball excludes it. For endomorphisms, every non-complex-linear matrix has a nonzero commutator entry, and an over-bound norm is eventually separated from the cutoff. Exact tests remain necessary for acceptance at a finite stage. In these replays the independent lattices supply that test.

## Exact acceptance and saturation

For a certified algebraic map $f:C\to E_0$ of degree $d$, its pullback differential evaluated on the periods gives a real matrix $F$ for $f_*$. The endomorphism $f^*f_*$ is $d$ times the orthogonal projector off $\ker F$, so its integral homology matrix is enclosed by
$$
d\,G^{-1}F^{\mathsf T}(FG^{-1}F^{\mathsf T})^{-1}F.
$$
The algebraic map is checked before this calculation. Because the true matrix is integral, an enclosure with a unique integer in each entry identifies it exactly. Products of these verified matrices fill the rank bound, checked using elliptic reduction fields or the split CM bound. The degree-three auxiliary isogenies are checked by substitution and a degree calculation. The supplementary CM maps use the exact isomorphism between the two models and their order-four automorphism.

Let $V$ be the rational span of the verified matrices once it reaches the rank bound. Then $V=\operatorname{End}^0(J)$. An element of $V$ that is integral on homology defines a holomorphic torus endomorphism. Algebraicity for complex abelian varieties and descent of their homomorphisms through algebraically closed extensions yield
$$
\operatorname{End}(J_{\overline{\mathbb Q}})=V\cap M_4(\mathbb Z).
$$
An integral double kernel computes the saturated intersection. The Rosati-symmetric part gives the Néron–Severi lattice, with exact form
$$
D\cdot D'=\frac{\operatorname{tr}(R)\operatorname{tr}(R')}{4}
-\frac{\operatorname{tr}(RR')}{2}.
$$
For each split example the torsion-certified basis maps into this lattice by an integral unimodular coordinate matrix. For each primitive isotropic class of positive degree, the exact integral kernel supplies the corresponding Plücker vector. The verifier compares the full set of these vectors with the period shortlist.

## Results

| Model | End rank | NS discriminant | Degree bound | Matched elliptic planes |
|---|---:|---:|---:|---:|
| $C$ | 2 | $-4$ | 150 | 2 |
| $C'$ | 2 | $-16$ | 150 | 2 |
| $T_1$ | 4 | $32$ | 150 | 52 |
| $T_2$ | 4 | $32$ | 150 | 52 |
| Degree-three gluing, identity | 4 | $6$ | 150 | 126 |
| Degree-three gluing, transposition | 4 | $24$ | 150 | 86 |
| Degree-three gluing, three-cycle | 4 | $24$ | 150 | 62 |
| CM gluing, $j=1728$ | 8 | $-16$ | 40 | 822 |
| Quaternionic model | 4 | $12$ | 60 | 0 |

For $y^2=x^5-x^4+4x^3-8x^2+5x-1$, the verifier replays the [exact correspondences](../qm-certificate/README.md) for $\alpha,\gamma$ and the two ordinary Frobenius fields excluding a CM square. The tangent actions under all four embeddings of $\mathbb Q(\zeta_8)$ recover the integral homology matrices. The order generated by $1,\alpha,\gamma,\alpha\gamma,\gamma\alpha,\alpha^2,\gamma^2$ is the full geometric ring. Its reduced trace pairing has determinant $-36$, yielding the maximal order of reduced discriminant $6$. The Néron–Severi form has no vector primitive at $3$ of square zero modulo $9$. Any nonzero rational isotropic vector could be scaled to such a vector, so this congruence also rules out elliptic quotients in every degree.

## Replay and scope

Run `python research/elliptic-quotients/rigorous-search/verify.py` from the repository root. The replay needs Python with `python-flint`, `sympy` and `cypari2` and uses the adjacent algebraic input records. It compares with [the recorded results](results.json). Select one or more cases by name, or pass `--lattices-only` for the exact lattice and spectrum replay without the two ambient comparisons.

The bounded searches prove that no true candidate in their range is omitted. Exact membership and the full lattices here use separately verified algebraic inputs and rank witnesses. Finding those inputs from every genus-two equation is a separate search problem. Costa–Mascot–Sijsling–Voight and Lombardo provide prior algorithms for geometric endomorphism rings, and Kani supplies the correspondence for elliptic subcovers. This study establishes the explicit certificate and comparisons on its stated inputs. [Costa–Mascot–Sijsling–Voight](https://doi.org/10.1090/mcom/3373) · [Lombardo](https://doi.org/10.1090/mcom/3358) · [Kani](https://doi.org/10.1016/j.jnt.2018.05.011)
