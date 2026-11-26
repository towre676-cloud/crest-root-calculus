# Curve models for the level-three twins

The discriminant −19 pair comes from two principal polarizations on one unpolarized CM abelian surface. The archived equations had coefficients in a degree-24 field. Both models now descend explicitly to
\[
F=\mathbb Q(\theta),\qquad \theta^4-\theta^3-3\theta^2+2\theta+4=0.
\]
Exact arithmetic verifies the descent, the shared elliptic quotient pair and distinct Galois orbits of the models. The connection between the recognized moduli and the analytic periods uses high-precision recognition. The opposite-secant equality of the two lattice spectra has a separate exact proof. [Visual companion](../../../docs/visuals/README.md)

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

## Exact comparisons

The first degree at which the two Humbert lattices differ is three, with **zero versus four** optimal classes. At every degree prime to three their ordinary square counts coincide by the shared internal residue line of the two opposite secants. Möbius inversion gives the same equality for primitive counts. [Incidence proof and exact replay](../opposite-secant-spectra/README.md)

Both models have the same unordered pair of degree-two elliptic quotient \(j\)-invariants. If \(d_0=v^2-4u+18v-27\), the two elementary symmetric functions in those \(j\)-values are
\[
\frac{256(u^2-2v^3-9uv+54v^2-27u)}{d_0},\qquad
\frac{256^2(v^2-3u+9v)^3}{d_0^2}.
\]
They agree exactly across the two curves. The resulting quadratic divides \(H_{-171}\) over \(F\). For an independent moduli check, \(I_2^5/I_{10}=16^5(v+15)^5/(64d_0^2)\) in the stated Clebsch normalization. Its two minimal polynomials are irreducible quartics and differ. Thus the two algebraic models have quartic moduli fields and lie in distinct Galois orbits.

## Reproduce

From this directory, with SymPy installed:

~~~sh
python -B descend.py --verify
~~~

The verifier checks both change-of-variables identities, the common CM quotient pair and the absolute invariant whose minimal polynomials separate the moduli orbits. The original theta and recognition chain is in [pipeline.py](pipeline.py), with controls in [validate.py](validate.py). It uses PARI through cypari2, mpmath, NumPy and python-flint. The archived [degree-24 worksheet](twins_-19.md) gives the handoff calculation.

All eight \(D=-43\) forms have candidate polarizations on explicit products of CM elliptic curves, recorded in [twin_lattices_-43.json](twin_lattices_-43.json). The historical search selected no moduli polynomial. Its Python floating-point cutoff underflowed at high precision. The updated driver compares residuals in PARI throughout. The field degrees and exact curve models are open.
