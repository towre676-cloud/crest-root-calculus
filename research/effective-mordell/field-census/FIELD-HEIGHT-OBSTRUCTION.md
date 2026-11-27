# How fast a field-dependent bound must grow

Write $C:y^2=x^6-1$ and $H(a/b)=\max(|a|,|b|)$ for a reduced rational abscissa. Let $D_K=|\operatorname{disc}(K)|$. There are infinitely many distinct real quadratic fields $K$ and points of $C(K)$ with rational abscissa for which
$$
        H(x)\ \gg\ \sqrt{D_K}\,(\log D_K)^{1/4}.
$$
The implied constant is absolute for this construction. In particular, no uniform upper bound $A\sqrt{D_K}(\log D_K)^\beta$ can hold for any fixed $A>0$ and $\beta<1/4$. This is a necessary growth condition on possible answers to Effective Mordell for one fixed curve. Larger field-dependent bounds remain possible.

## The auxiliary elliptic group

On $E':Y^2=X(X^2+21X+441)$, let $P=(84,882)$. Its double has coordinates $(225/16,7335/64)$. Shifting the abscissa to $X+7$ gives an integral short Weierstrass equation. Its double still has a nonintegral abscissa, so Nagell–Lutz proves that $P$ is nontorsion. The cubic has exactly one real root. The compact connected Lie group $E'(\mathbb R)$ is a circle.

The usual rational descent class is $X\bmod\mathbb Q^{\times 2}$ for nonexceptional points. Its homomorphism property follows from intersecting the cubic with a line: the product of the abscissas of the three intersections is the square of the line's constant coefficient. The vertical and zero-abscissa cases give the usual extensions to the identity and rational two-torsion. Hence the abscissa of every odd multiple $mP$ has class 21. Take $x=\sqrt{X/21}>0$ and $z=Y/(21x)$. Both are rational, and $z^2=21(x^4+x^2+1)$.

Take odd multiples approaching one of the real points of abscissa $21$, from the side $X>21$. Then $x=a/b>1$ approaches one, with coprime positive integers $a,b$. Its field is $K=\mathbb Q(\sqrt{21(a^2-b^2)})$, since $x^6-1=(x^2-1)z^2/21$. This is a genuine quadratic field because the complete rational point list has no $x>1$. The quadratic discriminant formula gives
$$
     D_K\le 84(a^2-b^2),\qquad
     \frac{D_K}{H(x)^2}\le84(1-x^{-2}).
$$
Thus already mere density of the odd multiples proves the unboundedness of $H(x)/\sqrt{D_K}$. Faltings' theorem guarantees that infinitely many fields occur: each fixed field has finitely many points on this genus-two curve.

## A quantitative subsequence

Choose an orientation of the real circle and use additive coordinates in $\mathbb R/\mathbb Z$. The odd multiples are $P+j(2P)$; the rotation parameter $\alpha$ of $2P$ is irrational. Choose a circle coordinate for the target with $X=21$, where $Y\ne0$. The coordinate $X$ is smooth with nonzero derivative there. We may choose the increasing side to have $X>21$.

Let $p/q$ be any continued-fraction convergent to $\alpha$. For $0\le j<q$, each $j\alpha$ lies within $1/q$ on the circle of $jp/q$. The latter points form the entire $1/q$ grid. Relative to the starting point $P$, choose a grid point at oriented distance between $2/q$ and $3/q$ beyond the target. The corresponding odd multiple has distance between $1/q$ and $4/q$ beyond it. As the convergent denominators increase these give infinitely many distinct odd $m=2j+1\le2q$, and $0<x-1\ll1/q\ll1/m$. This argument supplies the side of approach as well as the approximation rate.

For this fixed elliptic curve the standard canonical-height comparison gives
$$
    h(X(mP))=2m^2\widehat h(P)+O_{E'}(1).
$$
Since $X=21x^2$, it follows that $\log H(x)=m^2\widehat h(P)+O(1)$ with positive leading coefficient. The preceding discriminant estimate now gives $D_K\ll H(x)^2/m\ll H(x)^2/\sqrt{\log H(x)}$. Rearranging gives $H(x)/\sqrt{D_K}\gg (\log H(x))^{1/4}\gg(\log D_K)^{1/4}$ on an infinite subsequence with unbounded discriminant. As a finite arithmetic check, $213P$ gives $H(x)/\sqrt{D_K}>1.53$ and a rational abscissa of height 42,237 bits. The script proves this using integer arithmetic without factoring the large radicand. The infinite sequence follows from the proof.

## An independent rank obstruction

Write $k=\mathbb Q(\zeta_{12})$, $E:Y^2=X^3-1$. Given previously chosen squarefree radicands, take a prime $q>3$ dividing none of them. Put $n=q+1$ and $n^6-1=d m^2$ with $d>0$ squarefree. Because $v_q(n^6-1)=1$, the new squareclass $d$ is independent modulo all preceding radicands and the classes $-1,3$. Now $(n,m\sqrt d)\in C(\mathbb Q(\sqrt d))$. On the rational twists $E_t:Y^2=X^3-t^3$ we obtain
$$
   (dn^2,d^2m)\in E_d(\mathbb Q),\qquad
   (-d/n^2,d^2m/n^3)\in E_{-d}(\mathbb Q).
$$
Transport their points to $E$ over $k(\sqrt d)$. Their abscissas have height $2\log n$, so the duplication lower bound in the main note gives strictly positive canonical heights. Their two Galois characters correspond to $d,-d$. Applying the nonzero CM endomorphism $1+2[\zeta_3]$, conjugated to its negative by the nontrivial automorphism of $\mathbb Q(\sqrt{-3})$, gives two further independent character directions, with classes $-3d,3d$.

After $r$ steps the field $k(\sqrt{d_1},\ldots,\sqrt{d_r})$ has degree $4\cdot2^r$, and the $4r$ characters are distinct. Therefore the elliptic rank is at least $4r$. This supplies explicit field sequences testing methods that require a rank ceiling. The six-radical arithmetic receipt gives degree 256 and rank at least 24.
