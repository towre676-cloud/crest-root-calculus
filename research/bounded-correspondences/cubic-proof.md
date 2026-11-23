# Every cubic map totally ramified over elliptic infinity

Let the source be the smooth projective curve $y^2=F(x)$ with $\deg F=5$, and let the target be $v^2=G(u)$ with $\deg G=3$. Both branch polynomials have rational coefficients and are squarefree. The two odd-degree models specify rational Weierstrass points at infinity. The following reduction computes every rational or geometric map whose fibre above target infinity consists of source infinity. It reduces the complete coefficient incidence to two univariate equations of degrees four and six. [Implementation](cubic.py)

## Polynomial form and the degree bound

For a morphism $f$ carrying source infinity to target infinity, the pointed Abel maps satisfy $i_Df=f_*i_C$. The hyperelliptic involutions act as negation on the pointed Abel images, since $P+\iota(P)$ is linearly equivalent to twice infinity. The Abel map of the positive-genus target is injective. Consequently $f\iota_C=\iota_Df$. If the entire infinity fibre is a singleton, both coordinate pullbacks are regular away from infinity. The invariant and anti-invariant parts of $\mathbb Q[x,y]/(y^2-F)$ give

$$u=U(x),\qquad v=yV(x),\qquad G(U)=FV^2$$

with polynomials $U,V$ and $U$ nonconstant. Conversely, any such identity defines the required morphism on the smooth projective curves. Its degree is $d=\deg U$. These facts follow from the standard universal property of the Jacobian and extension of maps from smooth curves to proper curves. [Milne, Jacobian varieties, §§2 and 6](https://www.jmilne.org/math/articles/1986c.pdf)

Write $n=\deg F$, $m=\deg G$, and $e=\deg V$. Comparison of degrees gives $md=n+2e$. Differentiation gives

$$G'(U)U'=F'V^2+2FVV'.$$

Squarefreeness of $G$ gives $\gcd(G'(U),V)=1$, by substituting $U$ in a Bézout identity for $G,G'$. Thus $V\mid U'$, so $e\le d-1$ and

$$d\le\left\lfloor\frac{n-2}{m-2}\right\rfloor.$$

Geometrically this is Riemann–Hurwitz with ramification index $d$ at infinity. For $n=5,m=3$, the only possible pair is $d=3,e=2$. [Riemann–Hurwitz, Stacks Project 0C1B](https://stacks.math.columbia.edu/tag/0C1B)

## The centered equations

Translate the source coordinate and divide its branch polynomial by its leading coefficient $L$. Write the resulting equation as

$$y^2=L(x^5+Ax^3+Bx^2+Cx+D).$$

Make the rational target change

$$u'=g_3u+g_2/3,\qquad v'=g_3v,$$

where $G(u)=g_3u^3+g_2u^2+g_1u+g_0$. Its equation becomes

$$v'^2=u'^3+Pu'+Q,$$

$$P=g_3g_1-g_2^2/3,\qquad Q=2g_2^3/27-g_3g_2g_1/3+g_3^2g_0.$$

Every map exists only if

$$C=12A^2/49,\qquad D=4AB/7.$$

Under those two conditions define

$$U_0=x^3+(3A/7)x+B/3,\qquad V_0=x^2+A/7,$$

$$P_0=4A^3/343-B^2/3,\qquad Q_0=8A^3B/1029+2B^3/27.$$

There is an exact polynomial identity

$$U_0^3+P_0U_0+Q_0=(x^5+Ax^3+Bx^2+Cx+D)V_0^2.$$

All maps in the normalized coordinates are

$$u'=Lq^2U_0,\qquad v'=yLq^3V_0,$$

where $q$ runs through the solutions of

$$P=L^2P_0q^4,\qquad Q=L^3Q_0q^6.$$

Rational values of $q$ give exactly the rational maps. Every algebraic value gives one geometric map. Smoothness of the target excludes $q=0$, and different nonzero values give different signed maps. In particular, the maps with $q$ and $-q$ have opposite vertical coordinate.

To prove exhaustion, write $U=ax^3+bx^2+cx+d$. Since $V\mid U'$ and both have degree two, $V=kU'$ for a nonzero scalar $k$. Comparing leading coefficients gives $a=9Lk^2$. The coefficients in degrees eight, seven and six then give $b=0$, $c=3Aa/7$ and $d=Ba/3$. The degree-five and degree-four equations force the two conditions on $C,D$. The remaining coefficients give $P=a^2P_0$ and $Q=a^3Q_0$. Taking $q=3k$ gives the stated formulas. Equivalently, write $U=aU_0,V=sV_0$, so $s^2=a^3/L$. Then $q=s/a$ gives $a=Lq^2$ and $s=Lq^3$. This also proves completeness for arbitrary nonsquare rational $L$.

The source and target discriminants in the monic centered identity are

$$\operatorname{disc}(F/L)=\frac{432A(4A^3+343B^2)^3}{7^{10}},\qquad -4P_0^3-27Q_0^2=-\frac{16A^3(4A^3+343B^2)^2}{7^9}.$$

Their nonvanishing agrees. The scale equations have at most six geometric solutions and at most two rational solutions. When maps exist, the generic geometric count is two, rising to four or six at the two elliptic automorphism loci when the fixed target has the corresponding invariant.

## A transparent family

For $c\ne\pm2$ take

$$F_c=(x^2-4)(x^3-3x-c),\qquad G_c=(u^2-4)(u-c).$$

The identity

$$G_c(x^3-3x)=F_c(x)(x^2-1)^2$$

gives the two rational maps $u=x^3-3x$, $v=\pm y(x^2-1)$. The preceding classification proves they exhaust the rational maps with the stated infinity fibre. At $c=0$ there are four geometric maps to the fixed target, while the rational count stays two. The example with $A=21,B=18$ instead has target $v^2=u^3+1728$ and six geometric maps, again two rational maps. A quadratic twist of a fixed target can leave geometric maps while eliminating rational ones. The implementation handles these cases through the same scale equations.

The result uses classical Jacobian and ramification theory and explicit polynomial elimination. Totally ramified elliptic covers have an established literature, including Rubinstein–Salzedo's construction for every genus. The implemented reduction supplies an exact coefficient classification and a small independently checkable interface for this program's reconstruction engine. [Rubinstein–Salzedo, *Covers of elliptic curves with unique, totally ramified branch points*, 2013](https://doi.org/10.1002/mana.201200285)
