# The first elliptic quotient from the equation alone

A repeated elliptic factor can give elliptic maps even when the curve has no elliptic involution. This study searches for the first quotient directly from the coefficients using periods and then certifies the resulting map by exact substitution. The two repeated-factor Richelot neighbours \(T_1,T_2\) each have four optimal maps of degree four. The exact Richelot theorem below proves minimality and determines their spectra in every degree. For arbitrary inputs, periods supply candidates and the completeness of the computed spectrum remains numerical.

## Method

**Homology and polarization.** After the Möbius change $x=z_0+1/\tilde x$ with complex $z_0$ maximizing clearance, the six branch points are in general position. Fix a branch point $e_1$ and the star of straight segments to the others. In the local coordinate $z=\sqrt{x-e_1}$, the loop $\gamma_j$ lifting $[e_1,e_j]$ is a straight line through $0$ oriented along $\arg(e_j-e_1)/2$. With angles measured from a cut missing every segment, this gives $\gamma_j\cdot\gamma_k=+1$ whenever $j$ precedes $k$. Its period is $2\int_{e_1}^{e_j}\omega_+$. The parametrization $x=e_1+(e_j-e_1)\sin^2(t/2)$ makes the integrand analytic, and a nearest-point branch convention makes every square root continuous. Two identities test this independently of each other. The five loops satisfy the alternating relation $\gamma_1-\gamma_2+\gamma_3-\gamma_4+\gamma_5=0$, and the resulting form satisfies Riemann's first relation $\Pi E^{-1}\Pi^{\mathsf T}=0$ with $\det E=1$. Both hold to between $10^{-71}$ and $10^{-73}$ on every input.

**Endomorphisms and subvarieties.** In lattice coordinates the complex structure is $J_\Lambda=P^{-1}J_0P$. Endomorphisms are the integer $R$ with $J_\Lambda R=RJ_\Lambda$, and the Néron–Severi group consists of the integral alternating $E'$ with $J_\Lambda^{\mathsf T}E'J_\Lambda=E'$. Both are found by LLL at 70 digits, and the candidates are rechecked at the working precision. The lattice $L$ of an elliptic subvariety of degree $d$ has $E|_L=d$ times a generator. For any $u\in L$, $H(u,u)=d|u|^2/\operatorname{covol}L$, and Hermite's bound for the shortest $u$ gives $H(u,u)\le2d/\sqrt3$. For the exact form this proves that a shortest-vector search suffices. The implementation uses the computed form and `qfminim` to propose the directions. A vector $u$ lies in one exactly when the Plücker vector of $u\wedge J_\Lambda u$ is rational, and then $d=|\langle E,\pi\rangle|$ is an exact integer.

**Maps.** Let $L^\perp=\{v:E(v,L)=0\}$ be the lattice of the complementary subvariety. The optimal quotient is $z\mapsto\ell(z)$, where $\ell=a\,dx/y+b\,x\,dx/y$ vanishes on $\mathbb C\,\Pi(L^\perp)$, and its target lattice is $\ell(\Pi\mathbb Z^4)$. With a Weierstrass base point, $X=\wp(\ell(\mathrm{AJ}(P)))$ is a rational function $R(x)$ of degree $d$. It is interpolated from Abel–Jacobi values at sample points chosen for path clearance. The coefficients, $(a{:}b)$, $g_2$ and $g_3$ are recognized in $\mathbb Q(i)$. Write \(R=N/D\) and \(L=a+bx\). The certificate checks
\[
f(N'D-ND')^2=L^2D(4N^3-g_2ND^2-g_3D^3).
\]
It gives the vertical coordinate \(Y=y(N'D-ND')/(LD^2)\) over the claimed field. The calculation verifies the differential as well as both map coordinates. It checks coprimality, nonconstancy and nonsingularity of the elliptic target, so the map degree is \(\deg R=d\).

## Results

| Curve | Period candidate rank for $\operatorname{End}$ | Minimal degree | Classes at minimum | Spectrum through $100$ | Exact certificates |
|---|---:|---:|---:|---|---|
| $C$ (torsion study) | 2 | 2 | 2 | 2 (2) | both, over $\mathbb Q$ |
| $C'$ (torsion study) | 2 | 4 | 2 | 4 (2) | both, over $\mathbb Q$ |
| $C_2$ (PR #17) | 4 | 2 | 4 | equals the PR #17 law through $150$ | — |
| $T_1$, no involution | 4 | 4 | 4 | 4 (4), 20 (8), 52 (8), 68 (8), 100 (8) | all four: two over $\mathbb Q$, two over $\mathbb Q(i)$ |
| $T_2$, no involution | 4 | 4 | 4 | 4 (4), 20 (8), 52 (8), 68 (8), 100 (8) | all four: two over $\mathbb Q$, two over $\mathbb Q(i)$ |

The two targets are
\[
\begin{aligned}
T_1:\quad y^2 &=154x^6-90x^5-1605x^4-900x^3+5268x^2+7200x-1792,\\
T_2:\quad y^2 &=28x^6-450x^5-1317x^4+900x^3+6420x^2+1440x-9856 .
\end{aligned}
\]
Each has exact maps with pulled-back differentials \(dx/y\), \(x\,dx/y\), and \(dx/y\pm(i/2)x\,dx/y\). The first two maps are defined over \(\mathbb Q\). The other two are exchanged by conjugation, as seen from their distinct differential lines. The independent exact replay checks all twelve stored maps for the four curves \(C,C',T_1,T_2\). For \(C,C'\) the first degrees also follow from the separate torsion study. The period enumeration on \(C_2\) reproduces the proved PR #17 spectrum class by class through degree 150.

## Exact spectrum of every smooth Richelot neighbour of \(C_2\)

Let \(\psi:J(C_2)\to J(T)\) be a Richelot isogeny to a smooth genus-two Jacobian, with kernel \(K\subset J(C_2)[2]\). Work over the algebraic closure, so the result also applies to quadratic twists. If \(B\subset J(C_2)\) is an elliptic subvariety of optimal degree \(d\), its image has optimal degree \(d'=2d/|K\cap B[2]|\), by restricting the identity \(\psi^*\Theta_T=2\Theta_{C_2}\). The correspondence between elliptic subvarieties is bijective by taking the connected component of the inverse image.

**The intersection is determined by one of four degree-two directions.** Use the presentation \(\Psi:E^2\to J(C_2)\) from [the isotypic-gluing theorem](../isotypic-gluing/README.md). Its kernel is the graph of a transposition on \(E[2]\). If the primitive direction \((a,b)\) has opposite parity, its map \(h(P)=(aP,bP)\) into the product retains trivial kernel under \(\Psi\). Consequently \(\Psi h(E[2])\) gives its image subgroup \(B[2]\), which is the same as for one of the two axes depending on which coordinate is odd.

If \(a,b\) are both odd, \(\Psi h\) has kernel \(H\), the fixed subgroup of order two of the transposition. Its image \(B[2]\) is
\[
\Psi h\bigl(\{P\in E : 2P\in H\}\bigr).
\]
The input uses only \(E[4]\), so it depends only on \((a,b)\) modulo 4. Multiplying both entries by \(-1\) does not change the group. The two possibilities are represented by \((1,1)\) and \((1,-1)\). Each of the four representative directions has degree two on \(C_2\).

The four representatives therefore determine all degree corrections. Their images on a smooth target cannot have degree one by the genus inequality, so each correction gives either degree two or four. Let \(q\) be the number of representatives with image of degree two. These are all the degree-two classes on \(T\), since every other source degree is at least ten and can shrink by at most a factor of two. A reduced elliptic involution of \(T\) gives two degree-two classes, one for each lift. The exact Gröbner stabilizer computation in [the bielliptic study](../bielliptic-reconstruction/README.md) therefore determines \(q\) from the target equation.

**All multiplicities.** Let \(m\ge1\) have every prime divisor congruent to \(1\bmod4\), and write \(\omega(m)\) for its number of distinct prime divisors. The PR #17 formula gives \(2^{\omega(m)+2}\) source classes at degree \(2m\). The two opposite-parity residue types exchange under swapping coordinates. The two odd types modulo 4 are carried to the two opposite-parity types by \((a,b)\mapsto((a+b)/2,(a-b)/2)\), which preserves this parametrization of \(2m\). Thus each of the four groups supplies \(2^{\omega(m)}\) classes. The full target multiplicities are
\[
\#\{d'=2m\}=q\,2^{\omega(m)},\qquad
\#\{d'=4m\}=(4-q)\,2^{\omega(m)} .
\]
All other degrees have multiplicity zero. The count through degree \(X\) has leading coefficient \((4+q)/(4\pi)\) and error \(O(\sqrt X\log(2+X))\), by applying the PR #17 count to the four equal parts with their respective scale factors.

The exact replay enumerates all fifteen pairings of the six branch points. The two with zero Richelot determinant are excluded. Among the thirteen smooth results, four pairings give \(q=0\), eight give \(q=2\), and one gives \(q=4\). These are counts of pairings. The four listed in [targets.json](targets.json), corresponding to \(T_1,T_2\) and their \(x\mapsto-x\) versions, all have \(q=0\). Hence their spectra are exactly doubled in degree, and all four stored maps at the minimum are optimal. This conclusion uses no approximate periods.

## What the computations certify

The polynomial identity certifies the recognized maps and their degrees. Exact Richelot transport establishes the complete spectra for these specific neighbours. For other inputs, the real period matrices, endomorphism and Néron–Severi candidates, and lattice membership tests are still computed at 70 digits. Hermite's bound is rigorous for the exact positive form. The numerical search does not by itself supply a generally certified absence test or multiplicity count. General equation-only certification would require certified period enclosures and a way to certify the lattice identifications, or independent algebraic bounds. Costa, Mascot, Sijsling and Voight already give rigorous endomorphism-ring algorithms, including numerical proposals followed by algebraic certification. The specific result here is the four-class finite-torsion transport law and its exact spectrum for this Richelot family. Literature priority for the law remains open. [Richelot construction](https://www.math.uzh.ch/sepp/magma-2.20.4-cr/html/text1441.htm) · [Costa, Mascot, Sijsling and Voight](https://doi.org/10.1090/mcom/3373)

## Replay

~~~sh
pip install sympy
python -B verify_exact.py

pip install cypari2 mpmath numpy
python -B replay.py --output period-replay.json
~~~

The exact verifier uses the archived [map receipts](replay.json) and [spectra](spectra100.json). The numerical replay rebuilds the periods and quotient candidates independently. Map coefficient recognition currently works over \(\mathbb Q(i)\). The complete equation-only solver for general genus-two inputs, including maps over other number fields, remains open.
