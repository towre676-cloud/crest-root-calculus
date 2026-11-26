# Curve models for the level-three twins

The discriminant −19 pair comes from two principal polarizations on one unpolarized CM abelian surface. The archived equations had coefficients in a degree-24 field. Both models now descend explicitly to
\[
F=\mathbb Q(\theta),\qquad \theta^4-\theta^3-3\theta^2+2\theta+4=0.
\]
Exact arithmetic verifies the descent, the shared elliptic quotient pair and distinct Galois orbits of the models. The connection between the explicit equations and the twin lattices is now certified without recognition: Arb enclosures and exact lattice algebra recover each equation's refined Humbert form, as described in the section on the certified link below. The opposite-secant equality of the two lattice spectra has a separate exact proof. [Visual companion](../../../docs/visuals/README.md)

## From the lattice to the archived equations

Shioda–Mitani gives products of CM elliptic tori as unpolarized presentations. The search finds both principal polarizations on \(\mathbb C/\mathcal O_3\times\mathbb C/\mathfrak a_5\), with \(\mathcal O_3\) of conductor three and \(\mathfrak a_5=\langle5,(3+3\sqrt{-19})/2\rangle\). The refined Humbert matrices are
\[
Q_1=\begin{pmatrix}4&-2&0\\-2&12&-4\\0&-4&17\end{pmatrix},
\qquad Q_2=\begin{pmatrix}4&0&2\\0&5&-2\\2&-2&36\end{pmatrix}.
\]
A symplectic basis gives the period matrices. Siegel reduction, theta constants and Rosenhain coordinates give numerical moduli values. The chain was cross-checked on three known curves using Arb period enclosures, with relative agreement to \(10^{-54}\). The recognition step uses the residual and coefficient height, with an additional check on the size of the reduced field polynomial.

The archived equations have the form \(y^2=x^6-s_1x^4+s_2x^2-1\), with dihedral invariants \(u=s_1^3+s_2^3,\ v=s_1s_2\). Their exact arithmetic comparison with the recognized Igusa invariants takes place over the field with defining polynomial \(\eta^{24}+\eta^{18}+27\eta^{12}+163\eta^6+64\). [Lattice records](twin_lattices_-19.json) · [Recognized invariants](moduli_-19.json) · [Degree-24 coefficients](equations_-19_reduced.json)

## A model over the quartic field

Put \(\Delta=u^2-4v^3\). For \(v\Delta\ne0\), the archived sextic is geometrically isomorphic to
\[
Y^2=aX^6+bX^5+cX^4+dX^3+\Delta cX^2+\Delta^2 bX+\Delta^3a ,
\]
with
\[
a=u-v^2,\quad b=3\Delta+3u^2-2uv^2,\quad
c=\Delta(15u+v^2),\quad d=2\Delta(5\Delta+5u^2+2uv^2).
\]
The exact coefficients of \(u,v\) in the common quartic basis are in [dihedral_-19.json](dihedral_-19.json). The [expanded sextics and field embeddings](descent_-19.json) are constant-first.

For the change of variables set \(s=s_1^3,\ t=s_2^3,\ \delta=s-t,\ P=X+\delta,\ Q=X-\delta\). Substitute \(x=(s_1^2/v)P/Q\) and \(Y=tQ^3 y/\sqrt{\delta}\). Direct multiplication gives the right-hand side as
\[
\frac{s^2P^6-sv^2P^4Q^2+tv^2P^2Q^4-t^2Q^6}{\delta}.
\]
Expanding with \(st=v^3\) and \(\delta^2=\Delta\) gives the quartic-field equation. The two substitutions are checked coefficient by coefficient in the original degree-24 fields, after clearing the square-root factor. The existence of a field-of-moduli model in the extra-involution case is classical. [Cardona–Quer](https://arxiv.org/abs/math/0207015)

## Certified identification of the two models

Put \(g_1(U)=U^3-s_1U^2+s_2U-1\) and \(g_2(U)=U^3-s_2U^2+s_1U-1\). The degree-two maps \((x,y)\mapsto(x^2,y)\) and \((x,y)\mapsto(x^{-2},y/x^3)\) land on \(Y^2=g_1(U)\) and \(Y^2=-g_2(U)\), respectively. If \(\rho_i\) is a root of \(g_1\), the two torsion points \((\rho_i,0)\) and \((\rho_i^{-1},0)\) pull back to the same pair of Weierstrass points. The first pullback subtracts the fibre at \(x=\infty\), and the second subtracts the fibre at \(x=0\). Their difference is \(\operatorname{div}(x)\). Thus the kernel of \(E_1\times E_2\to\operatorname{Jac}(C)\) is the graph of \(\psi(\rho_i,0)=(\rho_i^{-1},0)\), and the polarization pulls back to \(2(\theta_1+\theta_2)\). [Howe–Leprévost–Poonen](https://arxiv.org/abs/math/9809210)

The certificate verifies exactly that the quotient \(j\)-invariants satisfy the ring class polynomial \(H_{-171}\). This proves that each is the modular \(j\) of an ideal lattice of \(\mathcal O_3\). Arb encloses the four possible modular \(j\)-values and singles out one for each quotient. After centring each elliptic cubic and putting it in the form \(Y^2=4X^3-g_2 X-g_3\), the invariants determine the \(\wp\)-coordinate scaling \(\lambda^{-2}\) by the ratio \((g_3/g_3^\tau)(g_2^\tau/g_2)\). Here \(g_2^\tau,g_3^\tau\) are the invariants of the selected ideal lattice. Arb excludes five of the six bijections between the cubic roots and the three scaled half-period values. It also certifies the reciprocal root pairing. Existence of both matchings follows from the equality of \(j\)-invariants and the reciprocal cubic identity, respectively.

The matched half-periods generate the graph subgroup. Adjoining their lifts to the product lattice gives the Jacobian lattice, on which the doubled product polarization is integral and unimodular. The code computes its refined Humbert form in exact arithmetic. It checks an integral isometry to the corresponding opposite-secant lattice for each model, excludes an isometry to the other, and checks the first discrepancy in quotient counts. The replay performs these checks at all 24 isolated embeddings of the original equation field. The exact descent transfers the outcome to the quartic models. [Certificate](certified_link_-19.json) · [Verifier](certify_link.py)

## Exact comparisons

The first degree at which the two Humbert lattices differ is three, with **zero versus four** optimal classes. At every degree prime to three their ordinary square counts coincide by the shared internal residue line of the two opposite secants. Möbius inversion gives the same equality for primitive counts. [Incidence proof and exact replay](../opposite-secant-spectra/README.md)

Both models have the same unordered pair of degree-two elliptic quotient \(j\)-invariants. If \(d_0=v^2-4u+18v-27\), the two elementary symmetric functions in those \(j\)-values are
\[
\frac{256(u^2-2v^3-9uv+54v^2-27u)}{d_0},\qquad
\frac{256^2(v^2-3u+9v)^3}{d_0^2}.
\]
They agree exactly across the two curves. The resulting quadratic divides \(H_{-171}\) over \(F\). For an independent moduli check, \(I_2^5/I_{10}=16^5(v+15)^5/(64d_0^2)\) in the stated Clebsch normalization. Its two minimal polynomials are irreducible quartics and differ. Thus the two algebraic models have quartic moduli fields and lie in distinct Galois orbits.


## Discriminant −43

The pair \(\mathtt{p3}\) has explicit models over \(F=\mathbb Q(\theta)\), with \(\theta^4-\theta^3+6\theta^2-\theta+7=0\). Their geometric optimal elliptic quotient counts agree in **every degree prime to 3**. Both have two classes in degree 2. Their first difference is in degree **6**, with **zero** classes for \(\mathtt{p3\_q1}\) and **four** for \(\mathtt{p3\_q2}\). [Quartic equations](models_certified_-43.json) · [Lattice certificate](certified_link_-43_p3.json) · [Opposite-secant theorem](../opposite-secant-spectra/README.md)

The certificate evaluates the descent formula and all three absolute Igusa invariants of each sextic in exact number-field arithmetic. The first invariant generates \(F\). The two elliptic quotient \(j\)-invariants satisfy the conductor-three ring class polynomial \(H_{-387}\) exactly, so the four CM ideal classes exhaust the possibilities. Arb excludes the wrong ideal classes and half-period bijections at each of the four embeddings. The resulting two-torsion gluing gives an integral Humbert form isometric to the intended archive form, and the verifier records its change of basis to the opposite-secant lattice. Primitive square representations give the quotient counts. [Exact verifier](verify_d43.py)

### Quartic descent for the six Mestre models

The six models in \(\mathtt{p0},\mathtt{p1},\mathtt{p2}\) have integral sextics over the same quartic field \(F\). The arrays list coefficients in ascending order and use \(w=\theta\). Each of the six stored Mestre conic points is checked by exact substitution. A pencil of lines through that point parametrises the conic and pulls the Mestre cubic back to a sextic. The replay recomputes its absolute Igusa invariants in \(F\). [Quartic sextics and points](mestre_quartic_models_-43.json) · [Descent replay](verify_mestre_descent_d43.py)

For the construction, each integral conic was minimised prime by prime using generators of the prime ideals, which exist because \(F\) has class number one. The reduction's radical permits a change of lattice and rescaling that lower the determinant valuation. In the rank-two case with valuation one, an isotropic line of the nonsingular binary block and the radical give the plane used for that rescaling. The six reduced matrices have unit determinants. Short vectors supply points and removing the fractional coefficient ideal gives the displayed integral models up to quadratic twist.

An independent test checks local solubility at 333 candidate finite places across the six conics. For an integral Gram matrix \(G\), the rational primes to test divide \(2N_{F/\mathbb Q}(\det G)\). At every other place a nonsingular point exists over the residue field and lifts by Hensel's lemma. For a diagonalisation \(\operatorname{diag}(d_0,d_1,d_2)\), the Hilbert symbol of \((-d_1/d_0,-d_2/d_0)\) is \(+1\) at every candidate place. At dyadic places the unit parts are reduced modulo \(\mathfrak p^{2e+1}\), where \(e=v_{\mathfrak p}(2)\). The field is totally imaginary. [Local certificate](mestre_obstruction_-43.json) · [Local replay](mestre_obstruction_d43.py)

### Exact lattice identification via Richelot paths

Two chains of Richelot isogenies connect the six quartic equations to the certified bielliptic models.

| Starting model | One step | Two steps | Three steps |
|---|---|---|---|
| $\mathtt{p3\_q1}$ | $\mathtt{p2\_q2}$ | $\mathtt{p0\_q2}$ | $\mathtt{p1\_q1}$ |
| $\mathtt{p3\_q2}$ | $\mathtt{p2\_q1}$ | $\mathtt{p0\_q1}$ | $\mathtt{p1\_q2}$ |

For each edge, the verifier multiplies three quadratics with distinct roots and checks its absolute Igusa triple against the source. It checks that the three quadratics are linearly independent, then checks the triple of the product of their Richelot brackets against the target. In each match the first invariant is a root of an irreducible quartic, and the other two are evaluated as exact polynomials in that same root. This certifies geometric isomorphism under a single field embedding for each triple. The edge witnesses therefore give isogenies up to Galois conjugacy, which suffices for scalar complex multiplication. [Exact Richelot witnesses](richelot_paths_-43.json) · [Richelot verifier](verify_richelot_d43.py)

Both factors of each bielliptic starting model have CM by the order of conductor three in \(\mathbb Q(\sqrt{-43})\). Scalar multiplication by \(3i\sqrt{43}\) thus preserves their product lattice. If an isogeny has kernel killed by \(e\), multiplication by \(e\) factors back through the isogeny. Composing with the scalar endomorphism and the original isogeny gives that scalar multiplied by \(e\) on the codomain. The bielliptic gluing and each Richelot step have kernels killed by two. Consequently, after \(s\) steps, \(6\cdot 2^s\sqrt{43}\,J\) is an **integral** homology matrix, where \(J\) represents multiplication by \(i\) on real homology.

Galois conjugation changes the tangent scalar at most by the sign of \(\sqrt{-43}\), so the integrality bound holds at every embedding. Arb encloses the matrix from each quartic equation at all four embeddings. Each entry's enclosure isolates one integer. The integrality theorem makes that integer the exact entry. The full geometric endomorphism ring is the integral centralizer of this matrix, computed by an integer kernel. Taking its Rosati-symmetric part gives the Néron–Severi lattice and hence the refined Humbert form. Every computed form matches exactly its intended archive form. An optional bounded Rosati search independently reproduces the same order. [Equation-to-lattice verifier](mestre_link.py) · [Exact results](mestre_link_quartic_-43.json)

All four pairs agree in every degree prime to \(3\) by the opposite-secant theorem. Their first differing geometric counts are:

| Pair | First distinguishing degree | \(\mathtt{q1}\) | \(\mathtt{q2}\) |
|---|---:|---:|---:|
| \(\mathtt{p0}\) | 3 | 0 | 2 |
| \(\mathtt{p1}\) | 3 | 2 | 0 |
| \(\mathtt{p2}\) | 3 | 2 | 0 |
| \(\mathtt{p3}\) | 6 | 0 | 4 |

The archived octic models and the initial recognition computation record how the equations were discovered. [Octic models](mestre_models_-43.json) · [Orbit discovery](orbit.py) · [Recognition](recognise_F.py)

## Discriminant −67

**All seven pairs have exact equations over the same quartic field and certified Humbert forms.** The field is
\[
F_{-67}=\mathbb Q(\theta),\qquad \theta^4-\theta^3-6\theta^2+5\theta+13=0 .
\]
The discovery computation locates the fourteen polarizations on products of conductor-three ideals and recovers their Igusa invariants from the CM orbits. The following checks work with the exact sextics. [Polarized lattices](twin_lattices_-67.json) · [Discovered invariants](invariants_F_-67.json) · [Recognition code](recognise_direct.py)

For the two bielliptic models in \(\mathtt{p6}\), the quotient \(j\)-invariants are roots of \(H_{-603}\) by exact substitution. The labelled two-torsion gluing yields the intended Humbert form at all four embeddings. The curves agree at every degree below \(9\), where their counts are zero and four. [Bielliptic equations and certificate](certified_bielliptic_-67.json)

Each of the twelve other models is an integral sextic over \(F_{-67}\). To construct the equations, conic minimization reduces the integral Gram determinants to units. A short-vector search yields a point on each conic, checked by substitution in the original form. A pencil through that point parametrizes the conic and gives the Mestre sextic. The replay also computes all three absolute Igusa invariants of every output over \(F_{-67}\). The first invariant generates the field in every case. [Quartic models and conic points](mestre_quartic_models_-67.json)

The two bielliptic models anchor twelve Richelot edges. All the edges have branch-root witnesses in a common explicit degree-24 field. The verifier checks the nonzero Richelot determinant and both invariant triples exactly at each edge. Thus the resulting paths give \(6\cdot 2^s\sqrt{67}\,J\) as an integral homology matrix for a model \(s\) steps from its anchor, with \(s\le5\). At all four embeddings of every quartic model, Arb isolates the integer entries. The integral centralizer gives the full endomorphism ring and the Rosati-symmetric lattice. Each resulting Humbert form matches precisely its own archive form. [Richelot witnesses](richelot_paths_-67.json) · [Exact verifier](verify_d67.py) · [Lattice and spectrum receipt](mestre_link_quartic_-67.json)

The exact lattices match the opposite-secant pairs. Within each pair the geometric optimal quotient counts agree in every degree prime to three. Three pairs first differ at degree \(3\), three at degree \(6\), and the bielliptic pair at degree \(9\).

## Discriminant −163

**All nineteen opposite-secant pairs have exact equations and certified Humbert forms.** In each pair the geometric optimal quotient counts agree in every degree prime to three. The pair that first differs at degree 21 has zero versus four optimal classes at that degree. Both curves have sextics over the common quartic field
\[
\mathbb Q(\theta),\qquad \theta^4-2\theta^3-9\theta^2+10\theta+28=0 .
\]
[Quartic sextics](quartic_models_-163.json) · [Exact descent](descend_d163.py) · [All 19 spectrum certificates](richelot_certified_-163.json)

### Starting equations from elliptic torsion

Take \(E:y^2=x^3-2174420x+1234136692+\tfrac14\). The certificate gives all roots of its two-division cubic and its three-division quartic in a single degree-24 field. Vélu's formulas for the four cyclic order-three kernels give cubics whose \(j\)-invariants satisfy \(H_{-1467}\) exactly. Their three roots are the images of the original two-torsion. For two such cubics and a chosen root bijection, write \(x\mapsto(ax+b)/(cx+d)\) for the corresponding projective transformation. With \(\alpha_j\) the roots of the first cubic, set \(u_j=\alpha_j+d/c\). Then
\[
C:\quad y^2=\prod_{j=1}^3(x^2-u_j).
\]
Its elliptic factors have branch roots \(u_j\) and \(1/u_j\), which are affine transforms of the matched cubic roots. The replay checks the reciprocal identity exactly. The two-torsion graph yields the principal polarization by the same argument used at discriminant −19.

The two starting equations in the twin census are \(\mathtt{p16\_q1},\mathtt{p16\_q2}\). Their glued lattices match the prescribed archive at all 24 embeddings. Exact embeddings place their dihedral invariants in the common quartic field. The descended equations have the same absolute Igusa invariants. [Torsion field](d163_torsion_field.json) · [Starting equations](models_certified_-163.json) · [Gluing replay](verify_d163.py)

### Algebraic equations along Richelot paths

A record in [the equation file](richelot_models_-163.json) specifies a parent, a Möbius chart and a matching of six branch roots. It gives each of the 36 other equations in factored form.

1. At a starting sextic take the roots in the order specified by the rational isolating rectangles. The auxiliary starting curve is another even sextic constructed from its recorded elliptic indices and torsion permutation. Its roots are the ordered pairs of square roots of its three \(u_j\).
2. The three indices of the chart select roots mapped to \(0,1,\infty\). Apply \(t\mapsto1/(t-2-3i)\), then list the images of \(0,1\), the three remaining roots in inherited order, and \(\infty\). Each matched pair gives a monic quadratic \(A,B,C\).
3. The child equation has right-hand side \([A,B][B,C][C,A]\) divided by its leading coefficient, with \([A,B]=A'B-AB'\). Its ordered roots come from the quadratic formula on each bracket in turn. The replay specifies the radical branches.

The exact equations use at most five Richelot steps and quadratic radicals. The verifier checks a nonzero Richelot determinant and six distinct roots at every step. Each starting Jacobian has scalar multiplication by \(6\sqrt{-163}\). After \(s\) Richelot steps the integral matrix \(6\cdot2^s\sqrt{163}\,J\) exists on homology, where \(J\) is the complex structure. Arb isolates its integer entries. Its integral centralizer gives the full endomorphism ring and its Rosati-symmetric part gives the Néron–Severi lattice. All 19 computed Humbert form pairs are matched to the opposite secants by exact integral isometries. [Richelot formula](https://msp.org/ant/2011/5-4/ant-v5-n4-s.pdf) · [Equation-to-lattice replay](verify_d163_paths.py)

| First differing degree | Pairs |
|---:|---:|
| 3 | 3 |
| 6 | 6 |
| 9 | 6 |
| 12 | 3 |
| 21 | 1 |

## Reproduce

From this directory, with SymPy, python-flint and cypari2 installed:

~~~sh
python -B descend.py --verify
python -B certify_link.py --verify
python -B verify_d43.py
python -B verify_mestre_d43.py
python -B mestre_obstruction_d43.py
python -B verify_mestre_descent_d43.py
python -B mestre_link.py
python -B verify_d67.py
python -B verify_d163.py
python -B descend_d163.py
python -B verify_d163_paths.py
~~~

The descent replay checks both change-of-variables identities, the common CM quotient pair and the absolute invariant that separates the moduli orbits. The gluing replay rebuilds the equation-to-lattice certificate. The original theta and recognition chain is in [pipeline.py](pipeline.py), with controls in [validate.py](validate.py). It uses PARI through cypari2, mpmath, NumPy and python-flint. The archived [degree-24 worksheet](twins_-19.md) gives the handoff calculation.
