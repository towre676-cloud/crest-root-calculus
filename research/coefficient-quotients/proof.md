# Polynomial quotients with an unknown target

Let $P\in\mathbb Q[x]\setminus\{0\}$ be an expanded polynomial. Write $P=A^2F$ with $A$ monic and $F$ squarefree, keeping the entire leading coefficient in $F$. The identity is unique by factor multiplicities. For a positive-genus normalization the vertical coordinate is $y=y_P/A$, where $y_P^2=P$. This first change has degree one.

Write $n=\deg F$ and $g(F)=\lfloor(n-1)/2\rfloor$ for $n\ge1$. Constants are handled as genus-zero controls. All positive-genus curves below mean their smooth projective models.

## Constant vertical multiplier

Search for $F=G(U)$, with $U$ monic, $U(0)=0$, $d=\deg U\ge2$ and $1\le g(G)<g(F)$. Since the target is unknown, absorb any nonzero constant multiplier into $G$. The degrees satisfy $n=md$, where $m=\deg G\ge3$. For every such divisor put $t=1/x$, $L=\operatorname{lc}(F)$, and use
$$
 U^\vee=1+u_1t+\cdots+u_{d-1}t^{d-1},\qquad
 t^nF(1/t)/L\equiv (U^\vee)^m\pmod {t^d}.
$$
The lower terms of $G$ begin at order $t^d$. The formal root with constant one is unique in characteristic zero and determines every coefficient of $U$. Successive Euclidean divisions by monic $U$ express $F$ in base $U$, with digits of degree less than $d$. All digits must be constant. Their constants determine $G$, including its leading unit, and the squarefreeness of $F$ forces that of $G$. This is a fixed-degree form of classical polynomial decomposition. [Kozen and Landau](https://www.cs.cornell.edu/kozen/Papers/poly.pdf)

A sieve prime $p>\max(d,m)$ must preserve the source degree and integrality of the coefficients. The formal root lifts coefficient by coefficient in $\mathbb Z_{(p)}$, since each equation divides by $m$, a unit. Monic long division stays in the same ring. A nonconstant digit modulo $p$ therefore rules out a rational composition. Survivors undergo the exact rational test. There is at most one normalized answer per divisor.

Over the common horizontal field $\mathbb Q(U)$, the source function field has degree $2d$ and the target function field has degree two. Thus $(x,y)\mapsto(U(x),y)$ has degree $d$, also geometrically.

## Odd-degree models and nonconstant multipliers

Let $F$ and the unknown $G$ have odd degree, with $n\ge5$ and target genus at least one. Consider only maps of degree at least two for which source infinity is the sole point above target infinity. The pointed Abel maps intertwine the pushforward on Jacobians. The hyperelliptic involutions act by negation on the pointed Abel images, and the target Abel map is injective. Consequently the map commutes with the involutions. Its coordinates have no poles on the affine source, so the invariant and anti-invariant parts of its affine coordinate ring give
$$
 u=U(x),\qquad v=yV(x),\qquad G(U)=F V^2.
$$
Conversely such an identity with nonconstant $U$ gives the required morphism. The function-field degree calculation above still applies. Pointed target coordinate changes allow $U$ monic with constant zero and $V$ monic. These changes are exactly $u\mapsto au+b,\ v\mapsto cv$ between pointed odd hyperelliptic presentations: the invariant functions of pole order at most two are spanned by $1,u$, and the anti-invariant functions of pole order at most $\deg G$ are spanned by $v$. Thus the three normalizations identify a unique presentation per pointed target-isomorphism class. In that presentation $\operatorname{lc}(G)=L$.

Put $e=\deg V$, $m=\deg G$. Degrees give $md=n+2e$. Differentiation gives
$$
 G'(U)U'=F'V^2+2FVV'.
$$
Since $G$ is squarefree, $G'(U)$ and $V$ are coprime. Hence $V\mid U'$, $0\le e\le d-1$, and $d\le\lfloor(n-2)/(m-2)\rfloor$. Both $m$ and $d$ are odd. The degree list is finite. The stratum $e=0$ is covered by the constant method.

For a fixed $(d,m,e)$, first treat the $e$ lower coefficients of $V$ as parameters. The formal root formula applied to $FV^2$ recovers $U$, and repeated monic division gives polynomial equations on those parameters for all digits to be constant. These operations commute with specialization because they divide only by nonzero rational constants. Target squarefreeness is automatic among solutions: if $a$ were a repeated root of $G$, a root of $U-a$ of multiplicity $k$ would satisfy $2\operatorname{ord}(V)\ge2k-1$. Thus $U-a\mid V$, contradicting $e<d$. The constant digits recover $G$ without an extra discriminant inequality.

## Reduce the parameter count using the derivative

Write $U'=dVW$ with $W$ monic of degree $w=d-1-e$. For $t=1/x$, set
$f(t)=t^nF(1/t)/L$, $v(t)=t^eV(1/t)$ and $b(t)=t^wW(1/t)$.
Then $U^\vee\equiv (fv^2)^{1/m}\pmod {t^d}$. For $1\le j\le e$, let $v_{<j}=1+v_1t+\cdots+v_{j-1}t^{j-1}$ and put
$$
 r_j=[t^j](fv_{<j}^{\,2})^{1/m},\qquad s_j=[t^j](v_{<j}b).
$$
The coefficient $u_j$ is $r_j+2v_j/m$. Comparing the reversed coefficient of order $j$ in $U'=dVW$ gives
$$
 \boxed{\displaystyle
 v_j=\frac{(d-j)r_j-ds_j}{d-2(d-j)/m}.}
$$
The denominator is positive for $m\ge3$ and $j\le e<d$. Thus every coefficient of $V$ is a polynomial over $\mathbb Q$ in the coefficients of $W$. The remaining derivative equations and the constant-digit conditions recover the original incidence. Conversely long division by monic $V$ recovers $W$ uniquely from every solution, so there is no loss of normalized maps.

Choose the smaller parameter list and work in $\min(e,d-1-e)$ variables. In the direct chart use the remainder of $U'$ modulo $V$. In the complementary chart first solve the top $e$ equations by the displayed recurrence. In either chart one can reduce later equations modulo the derivative conditions. When $e=d-1$, $W=1$, and the entire stratum is decided by univariate polynomial arithmetic without Gröbner elimination. When $w=1$, only one scalar remains. This bounds the number of variables. Polynomial degrees and coefficient sizes can still grow.

## Finiteness with an unknown elliptic target

Over $\mathbb C$, de Franchis–Severi gives finitely many nonconstant maps from a fixed smooth source to varying targets of genus at least two, up to target isomorphism. For elliptic targets we use the bound on $d$ and the chosen image of infinity.

Let $\Lambda=H^1(C,\mathbb Z)$, with the Hodge metric on $\Lambda_{\mathbb R}$. For $f:C\to E$, pull back an oriented symplectic integral basis $\alpha,\beta$ of $H^1(E,\mathbb Z)$. The pullbacks span a Hodge-stable real plane. On that plane the Hodge area agrees with the absolute cup product. So their exterior product is an integral vector of norm $\left|\int_C f^*\alpha\wedge f^*\beta\right|=d$ in the induced metric on $\bigwedge^2\Lambda_{\mathbb R}$. Bounded $d$ leaves finitely many lattice vectors. Each vector determines its plane and the index of the pullback lattice in the saturated lattice of that plane. There are finitely many sublattices of that finite index.

The lattice with its induced Hodge structure determines the target elliptic curve as a complex torus. Integrating the pullback holomorphic form from source infinity fixes the pointed map up to target isomorphism. Thus the bounded-degree pointed elliptic maps are finite too. An isomorphism intertwining two maps in this problem must take the image of source infinity to the image of source infinity, so the pointed normalizations preserve this finiteness. It follows that each coefficient incidence has finite geometric support, and its coordinate algebra is finite-dimensional even for a nonreduced presentation. The existing exact rational-point solver can reconstruct or exclude every normalized rational candidate.

The finiteness ingredients are classical. The explicit unknown-target incidence and complementary recurrence give the effective search in this scope. Maps with arbitrary rational horizontal coordinates and general infinity fibres have a wider scope. [Effective de Franchis–Severi background](https://math.mit.edu/~poonen/papers/finiteness.pdf) · [Pointed cubic predecessor](../bounded-correspondences/cubic-proof.md) · [Finite-algebra engine](../bounded-correspondences/README.md)

## An example with a proper derivative factor

Take
$$
 U=x^5+\tfrac54x^4-\tfrac{10}{3}x^3,\qquad
 V=x(x-1)(x+2),\qquad W=x,\qquad
 G(u)=u(u-\tfrac{44}{3})(u+\tfrac{13}{12}).
$$
Here $U'=5VW$, and $F=G(U)/V^2$ is squarefree of degree nine. Its curve has genus four, and the quotient has genus one. Given just the coefficients of $F$, the $(d,m,e)=(5,3,3)$ search recovers the map in one variable. The saturated genus-two cubic case instead has $W=1$.
