# Torsion denominators: proofs

**Epistemic status.** Theorem 1 and Corollaries 2–4 are proved here from standard facts that are cited by name: Tate's isogeny theorem over finite fields, the Poincaré complete reducibility theorem, the bijection between abelian subvarieties and symmetric idempotents (Birkenhake–Lange, *Complex Abelian Varieties*, §5.3, whose proof is characteristic-free), and the equal-exponent property of complementary elliptic subvarieties in a principally polarized abelian surface (Kani, *The number of curves of genus two with elliptic differentials*, J. reine angew. Math. 485, 1997). The identity $d=\min\{m\ge1:me\in\operatorname{End}\}$ is proved in [polarized-reduction](../elliptic-quotients/polarized-reduction/README.md#why-polarization-constrains-degree). Proposition 5 is proved here, and its instance is checked by exact symbolic computation. The numerical certificates are exact finite-field computations. The involution oracle used in the sweep is high-precision numerical and is labelled that way. Novelty against the literature has not been assessed. The divisibility $d\mid|b(a)|$ is the real-Weil-polynomial form of the Howe–Lauter gluing exponent bound. What is new here is the exact equality and the finite certificate that realizes it.

## Setting

Let $X/\mathbb F_p$ be a smooth projective curve of genus two, let $A=J_X$ carry its principal polarization, and let $k=\overline{\mathbb F}_p$. Write the Weil polynomial as $T^2h(T+p/T)$ with real Weil polynomial $h(s)=(s-a)(s-b)$. Here $a\ne b$ are integers, and the ratio of a Frobenius root in the $a$-factor to one in the $b$-factor is not a root of unity. Over $k$ this is equivalent to $\operatorname{tr}\pi^{120}$ differing on the two factors. A root-of-unity ratio of two Weil numbers in a quartic CM field has order $n$ with $\varphi(n)\le4$, and every such $n$ divides $120$. By Tate's theorem, $A_k$ is isogenous to $E_a\times E_b$ with $\operatorname{Hom}_k(E_a,E_b)=0$. Hence $A_k$ has exactly two elliptic subvarieties $B_a,B_b$, and $\operatorname{End}^0(A_k)=\operatorname{End}^0(E_a)\times\operatorname{End}^0(E_b)$. Let $F$ be the $p$-power Frobenius and $V$ the Verschiebung, so $FV=VF=[p]$. Both are integral, and the Rosati involution exchanges them. Put
$$
\varphi_a=F+V-b\in\operatorname{End}(A),\qquad N=|a-b|,\qquad e_a=\varphi_a/(a-b)\in\operatorname{End}^0(A).
$$

## Theorem 1 (exact torsion denominator)

*Let $d_a$ be the degree of the optimal map $X_k\to B_a$, equivalently the exponent of the restricted polarization. For every prime $\ell\ne p$,*
$$
v_\ell(d_a)=v_\ell(N)-s_\ell,\qquad s_\ell=\max\{\,s\le v_\ell(N):\ \varphi_a(A[\ell^s])=0\,\}.
$$
*The set of such $s$ is an initial segment of $\{0,\dots,v_\ell(N)\}$. Consequently $s_\ell$ is certified by $\varphi_a$ vanishing on generators of $A[\ell^{s_\ell}]$, plus, unless $s_\ell=v_\ell(N)$, one point $W\in A[\ell^{s_\ell+1}](k)$ with $\varphi_a(W)\ne0$.*

*Proof.* On $E_a$ the endomorphism $F+V=\pi+\bar\pi$ acts as $a$, and on $E_b$ as $b$. So $\varphi_a$ acts as $a-b$ on the first factor of $\operatorname{End}^0$ and as $0$ on the second. Thus $e_a$ is an idempotent with image $B_a$. It is symmetric because $(F+V)^\dagger=V+F$. A symmetric idempotent is determined by its image, so $e_a$ is the norm idempotent of $B_a$. By the cited identity, $d_a$ is the least $m\ge1$ with $me_a\in\operatorname{End}(A_k)$. Since $\operatorname{End}(A_k)$ is a lattice in $\operatorname{End}^0(A_k)$, membership is decided prime by prime in $\operatorname{End}\otimes\mathbb Z_{(\ell)}$. For $\ell\nmid N$, the element $e_a$ is already $\ell$-integral.

Fix $\ell\mid N$ with $\ell\ne p$, and let $\psi\in\operatorname{End}(A_k)$. Then $\psi\in\ell^s\operatorname{End}(A_k)$ if and only if $\psi(A[\ell^s])=0$. Necessity is clear. Conversely, $[\ell^s]$ is a separable isogeny with kernel $A[\ell^s]\subseteq\ker\psi$, so $\psi$ factors as $\chi\circ[\ell^s]$ with $\chi\in\operatorname{End}(A_k)$. The same holds in $\operatorname{End}\otimes\mathbb Z_{(\ell)}$: if $u\psi\in\ell^s\operatorname{End}$ with $\ell\nmid u$, then $u\psi$ kills $A[\ell^s]$, and $u$ acts invertibly there. Hence $me_a$ is $\ell$-integral exactly when $\varphi_a$ kills $A[\ell^{v_\ell(N)-v_\ell(m)}]$, and the least admissible $v_\ell(m)$ is $v_\ell(N)-s_\ell$. Vanishing on $A[\ell^{s}]$ implies vanishing on $A[\ell^{s'}]$ for $s'\le s$, which gives the initial-segment property. Finally, $A[\ell^s]$ is étale, so testing on $k$-points suffices. On $A(\mathbb F_{p^r})$ one has $F^r=1$, so $V=[p]\,F^{r-1}$ there. $\square$

**Remark (the residue characteristic).** The $p$-part of $d_a$ is invisible to étale torsion. It is governed by the action of $\varphi_a$ on the $p$-divisible group. Corollary 3 sidesteps it with a second prime.

**Remark (degree at least two).** A degree-one map between smooth projective curves of genus $\ge1$ is an isomorphism, so $d_a\ge2$. If $N$ is a power of $\ell$, a kill at level $v_\ell(N)$ is therefore impossible. A certified kill at level $v_\ell(N)-1$ then determines $d_a=\ell$ with no witness. This closure is used once, for $C$ at $p=23$, and it is recorded explicitly in the certificate.

## Corollary 2 (both blocks agree)

In a principally polarized abelian surface, an elliptic subvariety and its complement have the same exponent (Kani). Hence $d_a=d_b$. In every computed certificate the two blocks are certified independently and agree. This is a nontrivial consistency check on the implementation.

## Corollary 3 (characteristic zero, exact degree)

*Let $X/\mathbb Q$ have genus two, and suppose $J_{\overline{\mathbb Q}}$ has an elliptic subvariety. Let $p\ne q$ be primes of good reduction at which the reduction is multiplicity-one in the sense above. Then $J_{\overline{\mathbb Q}}\sim E\times E'$ with $E\not\sim E'$. The common optimal degree $n$ of $X$ onto its elliptic quotients satisfies*
$$
n/p^{v_p(n)}=d^{(p)}/p^{v_p(d^{(p)})},\qquad n/q^{v_q(n)}=d^{(q)}/q^{v_q(d^{(q)})},
$$
*so $n$ is determined. Every nonconstant map from $X$ to an elliptic curve has degree a positive multiple of $n$, so $n$ is the minimum geometric degree.*

*Proof.* By Poincaré reducibility $J_{\overline{\mathbb Q}}\sim E\times E'$. If $E\sim E'$, the reduction at $p$ would have a repeated geometric elliptic factor. The transfer of the prime-to-$p$ part through Tate modules and smooth proper base change is the argument in [polarized-reduction](../elliptic-quotients/polarized-reduction/README.md#why-polarization-constrains-degree). It applies because the specialized symmetric idempotent is the unique one onto its block. Corollary 2 removes any ambiguity about which characteristic-zero factor reduces to which block. A map $X\to E''$ factors, after translation, as an optimal map followed by an isogeny. $\square$

## Corollary 4 (Richelot pairs in one prime)

Suppose $X$ has an involution with elliptic quotient $E_1$ of degree two, and $X'$ is a non-product Richelot neighbour. Proposition 5 gives a map $X'\to E_1$ of degree four, so the optimal degree $n'$ divides $4$. Since $n'\ge2$, it is a power of two. One odd multiplicity-one prime therefore determines $n'\in\{2,4\}$ completely.

## Proposition 5 (the degree-four map)

*Let $\hat\psi:J_{X'}\to J_X$ be the dual Richelot isogeny, and let $q:J_X\to E_1$ be the optimal quotient of degree two. The composite $X'\to J_{X'}\to J_X\to E_1$ has degree four.*

*Proof.* For $h:J_{X'}\to E$, the degree of $X'\to E$ is the integer $hh^\dagger$. Richelot isogenies respect principal polarizations up to $2$, so $\hat\psi\hat\psi^\dagger=[2]$. Also $qq^\dagger=[2]$. Hence $hh^\dagger=q\hat\psi\hat\psi^\dagger q^\dagger=2qq^\dagger=4$. $\square$

Concretely, the correspondence $G_1(x)H_1(z)+G_2(x)H_2(z)=0$, $yv=G_1(x)H_1(z)(x-z)$, sends $(z,v)$ to two points $(x_i,y_i)$ of $X$. The map is $(z,v)\mapsto(x_1^2,y_1)+(x_2^2,y_2)$ on $E_1$. Every coordinate of the sum is symmetric in $x_1,x_2$ up to a factor $1/v$. For the pair in the README it gives
$$
u=\frac{289z^4-120z^3-3446z^2+8797}{(5z+13)^2(z^2-10z+13)},\qquad
Y=-\frac{528(z+7)(3z-1)(z^2-5)(2z^2+2z-7)(7z^2+18z+3)}{w\,(5z+13)^3(z^2-10z+13)}.
$$
It satisfies $Y^2=(u-1)(u-4)(u-9)$ identically on $C'$. Here $u$ has degree $4$ in $z$, and $z$ has degree $2$ on $C'$, so the map has degree $8/2=4$.

## Effective bounds

The field degree needed for $A[\ell^s]$ to be rational divides the order of $F$ in $\operatorname{GL}_4(\mathbb Z/\ell^s)$. That order is bounded by $|\operatorname{GL}_4(\mathbb Z/\ell^s)|$. The torsion level needed is at most $v_\ell(N)\le\log_\ell(4\sqrt p)$, because $0<|a-b|\le4\sqrt p$. These bounds are crude but finite. In the recorded certificates every kill level needed $r\le4$. Witnesses were found with $r\le8$.
