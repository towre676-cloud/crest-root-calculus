# Lifting the fixed-curve point certificates to variable nonhyperelliptic curves

Let $K$ be a number field with a certified list on the smooth projective genus-two curve $C_6:v^2=u^6-1$. The maps below are defined over $K$. Thus one recovers every $K$-point on their sources by factoring over the original field. This is a direct finite-map application of the point certificate.

## Power covers

For positive integers $m,n$ let $D_{m,n}$ be the smooth projective normalization of
\[
 y^{2m}=x^{6n}-1.
\]
Put $a=2m,\ b=6n,\ h=\gcd(a,b)$. The map $(x,y)\mapsto(x^n,y^m)$ lands on $C_6$. Over a finite affine input point $(u,v)$, the full affine fibre is the product of the sets of $K$-roots of $X^n-u$ and $Y^m-v$. Every remaining point is an infinity branch.

Over the algebraic closure the $x$-map has degree $a$, since each simple zero of $x^b-1$ gives a totally ramified place of degree $a$. There are $b$ such places and no other finite branch points. At $z=1/x=0$, divide $y$ by the formal unit $(1-z^b)^{1/a}$ whose constant term is one. The equation becomes $(y')^a=z^{-b}$. Its branches are indexed by the value $\eta=(y')^{a/h}z^{b/h}$ in $\mu_h$. For fixed $\eta$, the exponents $a/h,b/h$ are coprime, giving one branch of index $a/h$. The indexing is Galois equivariant, so a branch is $K$-rational precisely when $\eta\in\mu_h(K)$. Its coordinate is the limiting value of $y^{a/h}/x^{b/h}$. Riemann–Hurwitz gives
\[
             g(D_{m,n})=1+\frac{ab-a-b-h}{2}.
\]
When $m=1$, the curve is hyperelliptic. For $m\ge2$, the order-$2m$ automorphism $y\mapsto\zeta_{2m}y$ fixes all $6n$ points with $y=0$. On a hyperelliptic curve of genus at least two, an automorphism descends through the unique degree-two map to the line. The kernel has order two. An automorphism of order greater than two therefore induces a nonidentity projective linear transformation, with at most two fixed values and at most four fixed points upstairs. This proves geometric nonhyperellipticity.

Suppose the complete $C_6(K)$ consists of its two infinities, the branch points $(z,0)$ for $z^6=1$, and $(0,z)$ for $z^2=-1$ in $K$. This applies throughout the original multiquadratic field tower. The full list on $D_{m,n}$ consists of the branch labels $\mu_h(K)$, points $(x,0)$ with $x^b=1$, and points $(0,y)$ with $y^a=-1$. If $W=\#\mu(K)$ and $R=\gcd(a,W)$, its size is
\[
     \gcd(h,W)+\gcd(b,W)+R\,\mathbf1_{2\mid W/R}.
\]

Over each of $\mathbb Q(\sqrt7)$ and $\mathbb Q(\sqrt{-7})$, all the curves $D_{m,n}$ except $D_{1,1}$ have exactly four points, consisting of $x=\pm1,y=0$ and the two rational infinity branches. For positive $7$ the extra base points have $u=\pm2,v=\pm3\sqrt7$, and for negative $7$ they have $u=\pm1/2,v=\pm3\sqrt{-7}/8$. In both fields each prime over 3 is unramified, with valuation of $v$ equal to one, excluding $m>1$. For $m=1,n>1$ in the real field the prime above 2 gives $n\mid2$. The number $-2$ is negative and $\sqrt2\notin\mathbb Q(\sqrt7)$, excluding $n=2$. In the imaginary field 2 splits, so the valuation of each proposed $u$ is $-1$, excluding $n>1$.

## Smooth plane curves with independently varying coordinate polynomials

For every $n\ge2$ and rational $a,b$ with $|a|,|b|\le1/2$, put
\[
 P(X)=X^n+aX,\qquad Q(Y)=Y^{3n}+bY,\qquad A_{n,a,b}: P(X)^6-Q(Y)^2=1.
\]
Its projective closure is a smooth plane curve of degree $6n$. To check the affine chart, a stationary value of $P(X)^6$ arises at $P=0$ or at $X^{n-1}=-a/n$. In the latter case $P(X)=a(n-1)X/n$, hence $|P(X)|\le1/2$. The same argument applies to $Q(Y)^2$ with $3n$ in place of $n$. If both affine partial derivatives vanish, $|P^6-Q^2|\le 1/64+1/4<1$. The leading form at infinity is $X^{6n}-Y^{6n}$, with $6n$ simple slopes, all smooth. Smoothness and Bézout give geometric irreducibility. Adjunction gives genus $(6n-1)(6n-2)/2$. The canonical line bundle is $\mathcal O(6n-3)$, very ample, proving geometric nonhyperellipticity.

For each affine input $(u,v)\in C_6(K)$ take all $K$-roots of $P(X)-u$ and $Q(Y)-v$. The Cartesian products give precisely the affine $K$-points of $A_{n,a,b}$. Its $K$-points at infinity have slopes $Y/X\in\mu_{6n}(K)$. Thus a finite point algorithm for $C_6$ supplies one on every member, with both the genus and the field degree unbounded on the proved tower. A partial input list produces verified points with partial status.

Even at fixed $n$ and $b$, these plane curves give infinitely many geometric isomorphism classes as $a$ varies over rational $(0,1/2]$. For a critical point of $P$ its image $u_c$ satisfies
\[
       u_c^{n-1}=-\frac{(n-1)^{n-1}}{n^n}a^n.
\]
Choose any $v$ above $u_c$ and a root $y$ of $Q(y)=v$. Since $|u_c|<1$, we have $v\ne0$, and $u-u_c$ is a local parameter at $(u_c,v)$ on $C_6$. The preimage on the source is a ramification point since $P'(x)=0$, while $Q'(y)\ne0$ by smoothness. Infinitely many target branch values occur as $a$ varies. By de Franchis, a fixed smooth complex curve admits only finitely many nonconstant maps to fixed $C_6$. Finitely many source isomorphism classes would allow only finitely many branch values, a contradiction.

For $n=2,\ a=b=1/2$, this genus-55 curve has exactly its two infinities over $\mathbb Q$. Over $\mathbb Q(\sqrt{17})$ the two additional points have $x=(-1\pm\sqrt{17})/4,\ y=0$. Indeed the field has affine base points only $(\pm1,0)$, and the two discriminants of the horizontal quadratic are $17/4,-15/4$. The only $K$-root of $Y^6+Y/2=0$ is zero, since every other root has fifth power $-1/2$, impossible at a prime over 2 in this quadratic field. Over $\mathbb Q(\sqrt7)$ the additional candidate horizontal values $\pm2$ lead to quadratic discriminants $33/4,-31/4$, neither square. This field again has only the two infinities.

The accompanying script checks 45 complete field–curve lists and 50 independent genus-formula comparisons, and preserves partial status on an incomplete input. Run the script with the verify subcommand and the fixed-curve classifier available. For the genus-55 example use the plane subcommand with arguments '[-17,0,1]' 2 1/2 1/2.

Standard references include [Molin and Neurohr, Computing period matrices and the Abel–Jacobi map of superelliptic curves, §3](https://webusers.imj-prg.fr/~pascal.molin/pdf/periods.pdf) and [Riemann surfaces of second kind and effective finiteness theorems, de Franchis' theorem](https://doi.org/10.1007/s00209-022-03018-3).
