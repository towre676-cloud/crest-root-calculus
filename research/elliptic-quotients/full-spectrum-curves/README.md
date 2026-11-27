# Equations for the determinant-1,440 full-spectrum pair

The smallest pair in the [full-spectrum census](../opposite-secant-spectra/delay/README.md) has inequivalent refined Humbert forms \(A\) and \(B\) of Gram determinant \(1440\). Their elliptic quotient counts agree in **every** degree. This directory gives an explicit genus-two curve for each form. Each curve has two exact elliptic covers of degree three. Its full refined Humbert form is certified at all \(16\) complex embeddings of its coefficient field.

## The lattices

Both forms occur on the same product \(E_1\times E_2\), with \(E_1=\mathbb C/\mathbb Z[3\sqrt{-10}]\) and \(E_2=\mathbb C/\langle 9,\,3+3\sqrt{-10}\rangle\). Both factors have CM by \(\mathcal O_3=\mathbb Z[3\sqrt{-10}]\), the order of discriminant \(-360\) with class number \(8\). The two forms are two different principal polarisations on this product. [Search](find_d1440.py) · [Lattices](lattices_1440.json)

The exact moduli computation below identifies the degree-eight orbit of each model.

## The equations

Let \(F\) be the polredabs model of the ring class polynomial \(H_{-360}\), and let \(\theta\) be a root of it, so that \(j_1=\phi(\theta)\). The coefficients lie in the ring class field \(\mathcal H=\mathbb Q(\theta,\sqrt{-10})\) of \(\mathcal O_3\), of degree \(16\). In [models_1440.json](models_1440.json) they are written as polynomials in one generator of \(\mathcal H\), with the defining polynomial `field_polynomial_w`. For each of \(A\) and \(B\), the following data are exact.

The curve is
\[
C:\ y^2=Q(x)\,R(x),\qquad Q=x^3+x^2+q_1x+q_0 .
\]
The target of the first cover is \(E_1:\ Y^2=X^3+3j_1(1728-j_1)X+2j_1(1728-j_1)^2\), and the cover is
\[
\varphi_1=\Bigl(v+\frac{\kappa x}{Q(x)},\ \sqrt{c_1}\,\frac{S_1(x)\,y}{Q(x)^2}\Bigr).
\]
This is a degree-three map, certified by the polynomial identity \(Q^3g_1(v+\kappa x/Q)=c_1RS_1^2\).

The second cover is \(\varphi_2=\bigl(P_2/R,\ \sqrt{c_2}\,S_2y/R^2\bigr)\), with \(P_2=x^2+p_1x+p_0\). Its target is \(E_2:\ Y^2=g_2(X)\), where \(g_2=\operatorname{Res}_x(Q,XR-P_2)/\operatorname{Res}_x(Q,R)\). It is certified by the identity \(R^3g_2(P_2/R)=c_2QS_2^2\).

The two covers have complementary triple fibres. The fibre of \(\varphi_1\) over infinity consists of the Weierstrass points on \(Q=0\), and the fibre of \(\varphi_2\) consists of those on \(R=0\). In both curves \(p_1=2/3\).

The verifier checks in exact arithmetic that \(H_{-360}(j_1)=H_{-360}(j(E_2))=0\) and that \(j_1\ne j(E_2)\). It also checks that \(QR\) is squarefree of degree six, so the curve has genus two, and that both maps have degree three.

## Lattice certification

Both quotients have CM by exactly \(\mathcal O_3\). Put \(\chi=\varphi_1^*+\varphi_2^*:E_1\times E_2\to J(C)\) and \(\psi=(\varphi_{1*},\varphi_{2*})\). Then
\[
\psi\chi=\begin{pmatrix}3&\beta\\ \beta^\vee&3\end{pmatrix},\qquad
\begin{pmatrix}3&\beta\\ \beta^\vee&3\end{pmatrix}\begin{pmatrix}3&-\beta\\ -\beta^\vee&3\end{pmatrix}=(9-\deg\beta)\,I .
\]
So \(\ker\chi\) is killed by some \(e\in\{1,\dots,9\}\), and \(e\mid2520\). Each map is optimal since a degree-one map from a genus-two curve to an elliptic curve is impossible. Thus its pullback embeds its elliptic target. The distinct \(j\)-invariants imply distinct elliptic subvarieties, so \(\chi\) is an isogeny. Positivity of its pulled-back polarization gives \(9-\deg\beta>0\).

Let \(\chi'\) satisfy \(\chi\chi'=[e]\). Transporting the scalar \(3\sqrt{-10}\in\mathcal O_3\) by \(\chi\) and \(\chi'\) gives an endomorphism of \(J(C)\) whose tangent action is \(3e\sqrt{-10}\). Hence the homology matrix \(7560\sqrt{10}\,J\) is integral up to sign. Here \(J\) is the complex structure on real homology.

Arb encloses the period matrix of each equation at every embedding of \(\mathcal H\). Every entry of \(7560\sqrt{10}\,J\) contains exactly one integer, and its square is checked to be \(-10\cdot7560^2\). The integral centralizer of this matrix has rank \(8\), so it is the full geometric endomorphism ring. Its Rosati-symmetric part gives the Néron–Severi lattice and the refined Humbert form. The Rosati, centralizer and Humbert steps come from [mestre_link.py](../twin-curves/mestre_link.py).

At all \(16\) embeddings the resulting form is isometric to its target and not to the other form. [Certificate](certificate_1440.json)

| Curve | Census child Gram | Certified reduced form | Degree-3 classes | Degree-2 classes |
|---|---|---|---:|---:|
| \(C_A\) | \([45,18,21;\,18,12,10;\,21,10,17]\) | \([9,2,2;\,2,12,-4;\,2,-4,16]\) | 2 | 0 |
| \(C_B\) | \([45,3,36;\,3,5,4;\,36,4,36]\) | \([5,-1,1;\,-1,9,-1;\,1,-1,33]\) | 2 | 0 |

## Scope

Both sextics are defined over the degree-16 ring class field. Each has an exact octic field of moduli. The verifier evaluates the absolute Igusa invariants in the ring class field, checks that \(J_2^5/J_{10}\) has degree eight, and expresses the other two invariants in its power basis by rational linear algebra. Both fields of moduli are isomorphic to the field defined by \(x^8-5x^6+12x^4-5x^2+1\). Descent of either sextic to its field of moduli remains open.

The covers are defined over the degree-16 field if the targets are written as the twists \(Y_i^2=c_i g_i(X_i)\). In those coordinates the second component of each map is \(c_iS_i y/D_i^2\) with \(D_1=Q\), \(D_2=R\).

The discovery step solves the degree-three cover equations numerically on all \(8\) Shimura conjugates at 1000 digits. It then recognises six normalised coordinates of each curve in \(\mathcal H\). In this run every coefficient was stable between rational-approximation bounds of \(10^{333}\) and \(10^{250}\). The verifier checks the exact models and lattice identifications independently of that recognition. [Discovery](gluing33.py) · [Exact verifier](exact_models.py)

~~~sh
python -B find_d1440.py 25
TWIN_DPS=1000 python -B gluing33.py A   # and B; discovery, minutes each
python -B exact_models.py --verify     # exact identities and 16-embedding certificates per curve
~~~
