# Uniform equations for the opposite-secant family

Let \(f>157\) satisfy \(f\equiv181\bmod300\). Put
\[
 k=\frac{4f+1}{25},\qquad \tau=3\sqrt{-f},\quad \operatorname{Im}\tau>0.
\]
Use the CM elliptic curves
\[
 E=\mathbb C/\langle1,\tau\rangle,\qquad
 E_h=\mathbb C/\langle225,\tau-h\rangle,\quad h\in\{114,39\}.
\]
The first model glues \(E\) to \(E_{114}\) by matching the two ordered half-period bases. The second glues \(E\) to \(E_{39}\) by sending the first half-period to the first and the second to the sum. In each case the following construction gives an explicit sextic over \(\overline{\mathbb Q}\). Its lattice is exactly the corresponding opposite-secant lattice of the [delayed-difference family](README.md#an-infinite-family-with-an-unbounded-first-difference).

## The equations

For \(\sigma=\tau\) and \(\sigma=(\tau-h)/225\), put
\[
 j_\sigma=j(\sigma),\qquad v_\sigma=\frac{j_\sigma}{j_\sigma-1728},\qquad
 E_\sigma:\quad Y^2=X^3-3v_\sigma X-2v_\sigma.
\]
The two \(j_\sigma\) are singular moduli of the order of discriminant \(-36f\). Thus they are algebraic roots of the same ring class polynomial, specified by their displayed CM lattices. Each cubic has \(j\)-invariant \(j_\sigma\). Both Weierstrass invariants of the lattices are nonzero since the order differs from the exceptional orders of discriminants \(-3,-4\).

Label the cubic roots by the half-periods \(1/2,\sigma/2,(1+\sigma)/2\). More explicitly, if the normalized lattice has invariants \(g_{2,\sigma},g_{3,\sigma}\), its labelled root at \(P\) is
\[
 \frac{2g_{2,\sigma}}{3g_{3,\sigma}}\,\wp_\sigma(P).
\]
These are specified roots of the algebraic cubic. Indeed the change of abscissa between the two short Weierstrass equations scales both invariants by the corresponding powers. Dividing the identities for \(g_3\) and \(g_2\) gives the displayed multiplier.

Let \(\rho_i\) be the three labelled roots for \(E\) and \(\eta_i\) those for \(E_h\), with \(i=1,2,3\). Use respectively the permutation \(\pi=(1,2,3)\) in one-line notation and \(\pi=(1,3,2)\). Let
\[
 \frac{aX+b}{cX+d}
\]
be the unique projectivity carrying \(\rho_i\) to \(\eta_{\pi(i)}\). The two elliptic curves are nonisomorphic by the degree calculation below, so \(c\ne0\). An affine match would also carry their origins at infinity to each other. Put \(u_i=\rho_i+d/c\). The exact factored equation is
\[
                 C_h:\qquad y^2=\prod_{i=1}^3(x^2-u_i).
\]
All \(u_i\) are distinct and nonzero. The two elliptic quotient maps use \(x^2\) and \(1/x^2\), with vertical coordinates \(y\) and \(y/x^3\). Their cubic root sets are \(u_i\) and \(1/u_i\), which are affine transforms of \(\rho_i\) and \(\eta_{\pi(i)}\), respectively. For the second transform use
\[
 \eta_{\pi(i)}=\frac ac+\frac{bc-ad}{c^2u_i}.
\]
Over \(\overline{\mathbb Q}\) the scalar factors in the cubic equations are absorbed by scaling the vertical coordinates. The two-torsion pullbacks give the graph of the indicated match, and the polarization pulls back to twice the product polarization. This is the [same bielliptic gluing calculation](../twin-curves/README.md#certified-identification-of-the-two-models) used for the finite twin models.

## Integral lattice calculation

Write \(O=\langle1,\tau\rangle\) and \(\mathfrak a_h=\langle225,\tau-h\rangle\). With
\[
 c_h=\frac{h^2+9f}{225}
 =\begin{cases}(k+231)/4 & h=114,\\ (k+27)/4 & h=39,\end{cases}
\]
these are proper invertible ideals of \(O\). Indeed the forms \((225,-2h,c_h)\) are integral, have discriminant \(-36f\) and are primitive. The only potential common prime is 3, and \(c_h\equiv2\bmod3\).

Every homomorphism \(E\to E_h\) is multiplication by
\(\alpha=225A+(\tau-h)B\) with \(A,B\in\mathbb Z\). Its integral matrix in the two displayed bases is
\[
 T=\begin{pmatrix}A&hA-c_hB\\ B&225A-hB\end{pmatrix},\qquad
 \deg\alpha=225A^2-2hAB+c_hB^2.
\]
The degree form is \(225(A-hB/225)^2+fB^2/25\). No nonzero pair gives degree 1, establishing the nonisomorphism used in the equation construction.

On the product, express a Néron–Severi class by its two restriction degrees \(p,q\) and its off-diagonal homomorphism \(\alpha\). In the product homology basis its alternating matrix is
\[
 F=\begin{pmatrix}pI_2&T^{\mathsf T}I_2\\ I_2T&qI_2\end{pmatrix},
 \qquad I_2=\begin{pmatrix}0&1\\-1&0\end{pmatrix}.
\]
The glued lattice adjoins half the columns of \(W=\binom{1_2}{S}\), where
\[
 S=\begin{pmatrix}1&0\\0&1\end{pmatrix}\ (h=114),
 \qquad
 S=\begin{pmatrix}1&1\\0&1\end{pmatrix}\ (h=39).
\]
Integrality on the enlarged lattice requires \(FW\equiv0\bmod2\) and \(W^{\mathsf T}FW\equiv0\bmod4\). In both cases the first condition is \(B\equiv0\bmod2\) and \(p\equiv q\equiv A\bmod2\). The upper-right entry of \(W^{\mathsf T}FW\) is \(p+q+226A-(h+S_{12})B\). The second condition therefore gives \(p\equiv q\bmod4\). These conditions are sufficient since the product generators and added halves generate the lattice.

The descended principal polarization has coordinates \((2,2,0)\) in the product. Intersections on the product are four times those on the Jacobian. Thus the refined Humbert norm of the descended class is
\[
 (D\cdot\theta)^2-2D^2=\frac{(p-q)^2}{4}+\deg\alpha .
\]
Modulo the principal class, integral coordinates are \((y,A,B')\) with \(y=(p-q)/4\) and \(B=2B'\). In these coordinates the norm is \(4y^2+225A^2-4hAB'+4c_hB'^2\). The unimodular changes \((A,B')=(x-z,-z)\) for \(h=114\) and \((A,B')=(x,-z)\) for \(h=39\) give
\[
 225x^2+4y^2+6xz+kz^2,\qquad
 225x^2+4y^2+156xz+(k+27)z^2.
\]
These are the two lattice norms in the basis used for the secants. The construction therefore gives exact genus-two equations for every parameter of the delayed-difference family and its fixed-field Pell subsequence.

## Humbert's thesis and replay

In *Sur les courbes de genre un* (1885), §§73 and 88, Humbert describes fibres of elliptic functions by the sum of their points in the elliptic group. In degree two the two points have fixed sum. Choosing the origin at the branch point over infinity makes them \(P\) and \(-P\), giving the three labelled nonzero two-torsion points used here. Humbert also constructs a plane equation from relations among theta functions in §20. The uniform gluing and the comparison of quotient spectra use the additional polarization and lattice calculations above. [Thesis, printed p. 88](https://patrimoine.sorbonne-universite.fr/medias/67/5d/47/22/675d4722-bf7e-425b-be95-b0a559a93643/files/TH_000349_001_pdf.pdf#page=95)

Run from the repository root with cypari2 and python-flint installed.

~~~sh
python -B research/elliptic-quotients/opposite-secant-spectra/bielliptic.py
~~~

The replay checks the congruences modulo four exhaustively, reconstructs both Humbert forms independently for \(M=9,33,57,81\), and verifies the first difference at each \(M\) by exact representation counts. Arb verifies the elliptic root labels and encloses both sextics at every test parameter. At \(M=9\) it also checks containment of zero on evaluating the ring class polynomial at the three singular moduli. The uniform proof of algebraicity and the lattice identities is given above.
