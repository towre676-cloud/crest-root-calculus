# Certified absence for a quaternionic-multiplication Jacobian

**25 November 2026.** A geometrically simple abelian surface with quaternionic multiplication becomes geometrically non-simple at every good prime. The quaternion algebra survives specialization, and the endomorphism algebra of a geometrically simple abelian surface in finite characteristic cannot contain it. Frobenius therefore cannot certify that $C$ has no elliptic map, and the absence certificate has to be global. This study gives one for
$$
C:\ y^2=x^5-x^4+4x^3-8x^2+5x-1,
$$
one of the two QM curves in the LMFDB genus-two database (Booker–Sijsling–Sutherland–Voight–Yasaki, §6.4).

## The certificate

**1. Two endomorphisms as explicit morphisms.**
- *Numerical step.* Periods give integer matrices for the endomorphisms, with tangent matrices $M$ satisfying $M\Pi=\Pi R$. Two are used: $\alpha$ with $\alpha^2=2$ (Rosati trace 20), and $\gamma$ with $\gamma^2+\gamma-1=0$ (Rosati trace 6). For sample points $P$, numerical Jacobi inversion finds the divisor $D_P$ with $\mathrm{AJ}(D_P)=M\,\mathrm{AJ}(P)$. The inversion uses a double-precision Newton search with routed Abel–Jacobi paths, refined at 50 digits, and every inverse is verified at high precision.
- *Interpolation.* Since $\psi\circ\iota=-\psi$, the Mumford form of $\psi(P)=\alpha(P-W_0)$ is $u=T^2-s_1(x)T+s_0(x)$ and $v=y\,(g_1(x)T+g_0(x))$. These four rational functions are interpolated. For $\alpha$ they have degrees 8, 8, 14 and 14, with singular-value gaps of 40 orders of magnitude or more. Every coefficient is recognized in $L=\mathbb Q(\zeta_8)$, with denominators powers of 17.

**2. Exact verification.** Reduction modulo $u$ follows $T^k\equiv a_kT+b_k$, so the check is an identity in $L(x)$:
$$
f(x)\,(g_1T+g_0)^2-f(T)\equiv0\pmod u .
$$
It holds, so $\psi$ is a morphism $C\to\mathrm{Sym}^2C$. By the Albanese property, $\psi-\psi(W_0)$ extends to an endomorphism of $J_L$. Its action on invariant differentials is computed exactly by the trace formula, with $dT_i=(T_i\,ds_1-ds_0)/u'(T_i)$:
$$
\psi^*\Big(\sum_i T_i^k\,\frac{dT_i}{V_i}\Big)=H_1\,\frac{dx}{y},\qquad H_1T+H_0=T^k(Ts_1'-s_0')(g_1T+g_0)^{-1}\bmod u .
$$
$H_1$ collapses to $c_0+c_1x$ with constant coefficients. For $\alpha$ the result is
$$
M_\alpha=\begin{pmatrix}t^3+t^2+t&-t^3-2t^2\\2t^2+t&-t^3-t^2-t\end{pmatrix},\qquad t=\zeta_8,\qquad M_\alpha^2=2I \text{ exactly}.
$$

For $\gamma$ the four functions have degrees $2,2,5,5$. The second generator was chosen by Rosati trace because function degrees track it. The direct anticommuting generator $\beta$, with $\beta^2=-3$ and trace 36, was not resolved by 40 samples at any degree up to 18.

**3. Relations.** The following hold exactly in $L$, with $\beta=\alpha\gamma-\gamma\alpha$, which anticommutes with $\alpha$ because $\alpha^2$ is scalar:
$$
M_\alpha^2=2,\qquad M_\gamma^2+M_\gamma-1=0,\qquad M_\beta^2=-6,\qquad M_\alpha M_\beta+M_\beta M_\alpha=0,\qquad M_\beta\ne0 .
$$
Tangent representations are faithful in characteristic 0, so $\operatorname{End}^0(J_{\overline{\mathbb Q}})\supseteq(2,-6)_{\mathbb Q}$.

**4. Division.** The Hilbert symbols $(2,-6)_2=(2,-6)_3=-1$, so $(2,-6)_{\mathbb Q}$ is the division algebra of discriminant 6, as the literature records for this curve.

**5. Not a CM square.** If $J\sim E_K^2$, every ordinary reduction would have geometric Frobenius field $K$. Here the fields differ already at $p=7$ and $p=11$: $\mathbb Q(\sqrt{-10})$ and $\mathbb Q(\sqrt{-57})$.

**Conclusion.** A division quaternion algebra over $\mathbb Q$ embeds neither in a commutative algebra nor in $M_2(\mathbb Q)$, and step 5 excludes $M_2(K)$. Hence $\operatorname{End}^0(J_{\overline{\mathbb Q}})=(2,-6)_{\mathbb Q}$, so $J$ is geometrically simple and $C$ has no nonconstant map to an elliptic curve, in any degree.

[`certificate.json`](certificate.json) records every relation as `true` and the final verdict `no_elliptic_quotient: true`. [`functions_alpha.json`](functions_alpha.json) and [`functions_gamma.json`](functions_gamma.json) hold the exact Mumford functions, and `verified_*.json` the exact tangent matrices.

## Replay

~~~sh
python -B certificate.py
~~~

This starts from the exact rational functions of both correspondences. It rechecks membership in the symmetric square and recomputes the tangent matrices before checking the algebra relations, Hilbert symbols and the two fields. The result must agree with the receipt. The interpolation programs and samples record how the candidates were discovered. The finite-field distinction between QM and a CM square is established prior art. See [Lombardo, Algorithm 4.5](https://arxiv.org/pdf/1610.09674) and the correspondence methods in [Costa, Mascot, Sijsling and Voight](https://doi.org/10.1090/mcom/3373).
