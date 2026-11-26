# Rigorous rank-two zeros and pairwise projection bounds

The three twists identified in the [exact residual study](../exact-residual-spectrum/README.md) have **analytic rank two and Mordell–Weil rank two**. The corresponding untwisted curves have both ranks zero. The analytic result combines an exact central modular symbol with rigorous real intervals for the ordinary second derivatives.

| Cremona curve | Twist | Certified interval for $L''(1)$ of the twist | Certified interval for untwisted $L(1)$ |
| --- | ---: | --- | --- |
| $338\mathrm c1$ | $-8$ | $[11.699390116033,11.699390116034]$ | $[2.573914226693,2.573914226694]$ |
| $676\mathrm b1$ | $-8$ | $[7.288897137561,7.288897137562]$ | $[1.496541884049,1.496541884050]$ |
| $867\mathrm e1$ | $-3$ | $[3.736897421981,3.736897421982]$ | $[4.428336432042,4.428336432043]$ |

The two conductor-10816 twists have exact $a_3$ values $1$ and $-2$, respectively. They are distinct rational isogeny classes. The [central-order receipt](central.json) also gives the minimal twisted models, two rational points for each twist, and the coinciding descent bounds. The descent uses PARI's unconditional rank bounds with fixed points and search effort zero. [PARI elliptic-curve documentation](https://pari.math.u-bordeaux.fr/dochtml/html/Elliptic_curves.html#ellrank)

## The analytic order

For a weight-two elliptic newform of conductor $N$, set $\alpha=2\pi/\sqrt N$ and write
$$
 g(y)=\sum_{n\ge1}a_n e^{-\alpha n y},\qquad
 \Lambda(s)=\alpha^{-s}\Gamma(s)L(s).
$$
All six curves have sign $+1$. The Fricke relation is $g(1/y)=y^2g(y)$. Splitting the Mellin integral at $1$ gives
$$
 \Lambda(s)=\int_1^\infty g(y)(y^s+y^{2-s})\,\frac{dy}{y},
 \qquad
 \Lambda''(1)=2\int_1^\infty g(y)(\log y)^2\,dy .
$$
The exponential decay justifies differentiation under the integral. For each twist the exact central symbol gives $L(1)=0$, and the functional equation gives $\Lambda'(1)=0$. Therefore the ordinary $L'(1)$ vanishes and $L''(1)=\alpha\Lambda''(1)$. For each untwisted curve $L(1)=\alpha\Lambda(1)$. [Mellin method for special values](https://doi.org/10.1080/10586458.2004.10504528)

The verifier evaluates the first $K=1024$ Fourier terms over $1\le y\le Y=2048$ with complex ball integration. It passes the analyticity flag to the logarithm. For tail estimates, $|a_n|\le d(n)\sqrt n\le 2n$ and both $1$ and $(\log y)^2$ are at most $y^2$ on this interval and beyond. The discarded Fourier tail is bounded by
$$
 \frac{4e^{-\alpha(K+1)}}{1-e^{-\alpha}}
 \left(\frac1\alpha+\frac{2}{\alpha^2(K+1)}
                    +\frac{2}{\alpha^3(K+1)^2}\right).
$$
This comes from integrating $4n y^2e^{-\alpha n y}$ over $[1,\infty)$ for each omitted term. The integral tail beyond $Y$ is at most
$$
 \frac{4e^{-\alpha Y}}{(1-e^{-\alpha Y})^2}
 \left(\frac{Y^2}{\alpha}+\frac{2Y}{\alpha^2}+\frac{2}{\alpha^3}\right),
$$
using $\sum_{n\ge1}n e^{-\alpha n y}=e^{-\alpha y}/(1-e^{-\alpha y})^2$. For each twisted curve the sum of the displayed error bounds is below $10^{-20}$ before scaling by $\alpha$. The verifier adds both tails outward before rounding the reported bounds. [Rigorous ball integration](https://flintlib.org/doc/tutorials/numerics/integrals.html)

## The pairwise test

Work in the inverse-mass inner product on coefficient vectors. For a pair of gluing classes $i,j$ in one genus, set $c=e_i-e_j$. For a subspace $U$ define its energy on the pair as
$$
 \mathcal E_U(i,j)=\frac{\|P_Uc\|^2}{\|c\|^2}.
$$
The coprime spectra coincide exactly when $\mathcal E_V(i,j)=1$, where $V$ is the coprime square kernel. The visible fraction of squared norm is $1-\mathcal E_V(i,j)$. In the split genus, with the forced space $F$ and its orthogonal residual $R$ inside $V$, one has $\mathcal E_V=\mathcal E_F+\mathcal E_R$. For any finite coefficient prefix the kernel contains $V$, so its projector energy bounds the all-degree energy from above.

The verifier constructs the incidence basis of $F$ explicitly. Choose a base isotropic vertex orbit $O_0$ in each ambient lattice. For each other orbit $O_j$ take the function with value $|O_0|$ on $O_j$, value $-|O_j|$ on $O_0$, and zero elsewhere. On a split plane orbit $\mathcal E$ with $k_i$ endpoints in $O_i$ it gives the coefficient
$$
 |\mathcal E|\bigl(|O_0| k_j-|O_j| k_0\bigr).
$$
These are the sum-zero vertex functions of the [forced incidence construction](../forced-square-kernel/README.md#ambient-lattices-and-incidence-relations), converted to relation coefficients by the class masses, up to a common factor on each ambient block. The rational projector from the resulting basis agrees with the previous incidence diagonal formula. The verifier computes the inner products on actual pair differences, including off-diagonal terms.

The replay checks **223,958 same-genus gluing pairs across 28 levels**. It uses exact kernels at 18 levels and finite-count upper bounds at the other ten. Across the entire sweep the greatest upper bound on the hidden fraction is $59/120$ in the split genus at $(-8,13)$, for scalar parameters $6\omega$ and $1$. Its decomposition is
$$
 \frac{59}{120}=\frac16+\frac{13}{40}.
$$
The fraction $61/120$ is attained on that pair and bounds the visible fraction from below across the whole sweep. The largest upper bound on the residual-only fraction is $1/3$, attained in the nonsplit genus at $(-3,17)$. This measures the energy on a pair and can differ from the largest single-class diagonal.

For fixed $D$, the forced bound $\mathcal E_F(i,j)\le \rho_F(i)+\rho_F(j)\ll_D 1/\ell$ yields a pairwise uniform criterion. A bound $\mathcal E_R(i,j)\le 1-\epsilon_D$ for all gluing pairs at sufficiently large $\ell$, with $\epsilon_D>0$, suffices for eventual separation. The finite data establish no such asymptotic bound.

## Replay

~~~sh
python -m pip install cypari2 numpy python-flint
study=research/elliptic-quotients/residual-rank-certificate
python -B "$study/central.py" --check-symbols --verify "$study/central.json"
python -B "$study/pair_margin.py" --include-prefix --verify "$study/pairs.json"
~~~

The default central replay uses the exact central-symbol results already recorded in the [packet receipt](../exact-residual-spectrum/packets.json). The supplied independent numerical comparison is commit `78fb04ceb95dad6ab546334a401b61ec2c4c6dd4`.
