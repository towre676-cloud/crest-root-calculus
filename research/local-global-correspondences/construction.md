# Constructing the graph and its minimal fields

Write $C:y^2=f(x)$ for a smooth genus-two curve over $\mathbb Q$, with $f$ an integral quintic or sextic. The six branch points are the binary sextic's projective roots. A maximal isotropic subgroup of $J(C)[2]$ amounts to an unordered pairing of these points. The three pairs represent its three nonzero elements. A rational subgroup corresponds exactly to a Galois-stable pairing, even if Galois acts transitively on the three elements. The classical Richelot construction associates a principally polarized quotient to each pairing. A zero determinant in the quadratic-factor construction gives a product over the algebraic closure. [Bruin and Doerksen, §4](https://doi.org/10.4153/CJM-2011-039-3) · [Kohel and Smith, explicit graph](https://www.i2m.univ-amu.fr/perso/david.kohel/pub/ladic.pdf)

## Recovering pairings from a single split prime

For a quintic the program first makes the rational change $x=(nw+1)/w$ with $f(n)\ne0$, recording the matrix. Write $a$ for the leading coefficient after this change, and let $F(u)=a^5 f(u/a)$. This is monic integral. Its roots $u_i$ have absolute value at most $R=1+\max_{j<6}|[u^j]F|$ at every complex embedding. Choose a prime $p$ at which all six roots are distinct and rational. For each unordered pair let $s=u_i+u_j$, $t=u_i u_j$, and $\theta=s+kt$. Select $k$ so the fifteen $\theta$ values differ modulo $p$. Such a choice exists for sufficiently large split primes, since any two distinct unordered pairs exclude at most one value of $k$. Split primes exist by Chebotarev.

For a matching $M$ write
$$
m_M(T)=\prod_{i=1}^3(T-\theta_i),\qquad
h_M(T)=\sum_{i=1}^3 t_i\prod_{j\ne i}(T-\theta_j).
$$
The coefficients are algebraic integers. For a rational matching they are ordinary integers. With $H=2R+|k|R^2$, their absolute values are bounded by $B=\max((1+H)^3,3R^2(1+H)^2)$. Lift all six roots by Newton's method to a modulus greater than $2B$. Center-lifting the coefficients must recover every rational matching exactly. To reject accidental lifts, form $A=\mathbb Q[T]/m_M$ and recover $t=h_M/m_M'$, $s=T-kt$ in $A$. The program verifies the polynomial identity
$$
\operatorname{Norm}_{A/\mathbb Q}(x^2-(s/a)x+t/a^2)=f(x)/a .
$$
Squarefreeness of $f$ makes every successful decomposition a genuine matching. Distinct local matchings at the separating split prime cannot encode the same global one. This proves termination and completeness independently of a guessed coefficient cutoff. [Rational pairing recovery](partitions.py)

For *all* geometric matchings choose an integer $b$ such that the fifteen values $u_M=m_M(b)$ are distinct modulo $p$. Differences of the monic cubics have degree at most two, so this is possible for sufficiently large separating primes. The polynomial $R_{15}(Y)=\prod_M(Y-u_M)$ is integral by symmetric invariance. Its irreducible factors give the exact fields of definition of the fifteen kernels. A separate interpolation recovers each pairing and its target over the corresponding field. For any required algebraic-integer coefficient $d_M$ compute
$$
I(Y)=\sum_M d_M\prod_{N\ne M}(Y-u_N),\qquad d_M=I(u_M)/R_{15}'(u_M).
$$
The implementation lifts to more than twice an archimedean coefficient bound before reducing the answer in each number field. In detail $|u_M|\le L=(|b|+H)^3$. With $C=(1+R)^2$, all branch-quadratic coefficients are bounded by $C$, their determinant by $6C^3$, and their Wronskian coefficients by $4C^2$. The seven coefficients of $-aD\prod H_i$ are therefore bounded by $10,368|a|C^9$. Together with the preceding pairing bounds, $15\max(B,10,368|a|C^9)(1+L)^{14}$ bounds coefficients of every interpolation numerator. Division by the interpolated $D^2$ recovers the correct rational twist. The monic branch coordinate $u=ax$ corresponds to the rational curve model $v^2=aF(u)$, with $v=a^3y$. [All-orbit interpolation](kernel_fields.py)

The independent verifier works inside each exact field $\mathbb Q[u_M]$. It reconstructs the cubic pairing algebra and uses determinants of multiplication matrices for norms. It checks both source and target at seven arguments, sufficient for exact equality of their degree-six polynomials. For a product quotient it verifies the vanished determinant. [Independent field calculation](check_orbits.py)

## Global graph in the cubic pairing algebra

In a splitting field write $g_i(x)=x^2-s_ix+t_i$ for the *original-coordinate* quadratics, $f=c\prod_i g_i$, and $D=\det([1,-s_i,t_i]_{i=0}^2)$. Put
$$
H_i=g_j'g_k-g_jg_k',\quad K_i=H_i/D,\quad
R(z)=-cD^2\prod_i K_i(z)
$$
with cyclic indices. For $D\ne0$ one has $\sum_i g_i(x)H_i(z)=D(x-z)^2$ and $[H_i,H_j]=2Dg_k$. The target is $y'^2=R(z)$. An exact correspondence is
$$
g_0(x)H_0(z)+g_1(x)H_1(z)=0,\qquad yy'=c(x-z)g_0(x)H_0(z).
$$
The first equation gives $g_1H_1=-g_0H_0$ and $g_2H_2=D(x-z)^2$, verifying the square of the second equation on the two curves. A generic point has two images in either direction. Taking projective closures gives the cycle on the smooth projective curves.

For the induced maps on regular differentials, put $Q=g_0H_0+g_1H_1$. On a generic fibre with roots $z_i$ use $dz_i/dx=-Q_x/Q_z$. Partial fractions give $\sum_i z_i^j/((x-z_i)Q_z(x,z_i))=x^j/Q(x,x)$ for $j=0,1$. Here $Q(x,x)=-g_2(x)H_2(x)$ and $Q_x/(g_0H_0)=H_2/(g_0g_1)$ on $Q=0$. Substitution gives pullback $z^jdz/y'\mapsto x^jdx/y$. Reversing the projections gives $x^jdx/y\mapsto 2z^jdz/y'$ by the double-bracket identity. The transpose correspondence is the adjoint for the canonical polarizations. The composition acts by 2 on the Jacobian since differentials determine homomorphisms in characteristic zero. Thus the isogeny has degree four and pullback polarization twice the source polarization. The paired branch-point differences give its kernel by the Richelot correspondence theorem.

The normalized $K_i$ belong to the cubic algebra $A$. They can be recovered without choosing their order by $\operatorname{Tr}(g(x)K(z))=(x-z)^2$. For each omitted index $k$ put $e=g(x)K(z)\in A$, reduced to $e_0+e_1T+e_2T^2$. Then the first graph equation divided by $D$ is $Q_k=(x-z)^2-e(T)$. Let $d_0$ be the determinant of the coefficient columns $1,-s,t$ in the basis $1,T,T^2$. A Vandermonde evaluation gives $D^2=d_0^2\operatorname{disc}(m)$. On $Q_k$, the second graph equation becomes
$$
yy'= -\frac{c(x-z)d_0\operatorname{disc}(m)}{2m'(T)}
\big(e_1+e_2(\operatorname{Tr}(T)-T)\big).
$$
Indeed $D(e_i-e_j)=-d_0\operatorname{disc}(m)(e_1+e_2(\operatorname{Tr}(T)-T_k))/m'(T_k)$ for cyclic $i,j,k$. This formula defines each graph over the field of its omitted pair. The three conjugates all induce the same map, as their differential matrices are identical. Hence the homomorphism descends to the field of the unordered matching. Their rational orbit sum induces three times the descended homomorphism. The verifier checks the graph identity over **every** component of $A$. [Trace-dual algorithm and graph check](partitions.py)

## Identifying another sextic as the quotient

In genus two every curve isomorphism induces a Möbius transformation on the hyperelliptic base. Given target equation $v^2=g(w)$, the program finds all rational matrices $M=\left(\begin{smallmatrix}a&b\\c&d\end{smallmatrix}\right)$ with
$$
(cw+d)^6R((aw+b)/(cw+d))=s\,g(w),\quad s\in\mathbb Q^\times.
$$
A rational square $s$ gives the two signed rational maps. Otherwise this frame gives a map over $\mathbb Q(\sqrt s)$. The graph uses $\sqrt{s}\,yv=B(T,x,w)$, with the denominator of the transported cubic cleared in $B$. Strong Torelli for genus two identifies all polarized isomorphisms with curve isomorphisms, including the sign from the hyperelliptic involution. Thus exhaustive rational kernels and frames exhaust the requested rational $(2,2)$ isogenies.

The fast method constructs quadratic covariants by binary transvectants. Starting from the sextic $F$, set $I=(F,F)_4$, $L=(F,I)_4$, $M_2=(I,L)_2$, $N=(I,M_2)_2$. If the three quadratics form a basis and each has nonzero discriminant, their individual lines are equivariant under every projective identification. The polar form $B(q,r)=b_qb_r-2(a_qc_r+c_qa_r)$ scales by the square of the determinant. Consequently the squared ratios of the determinant-normalized scale factors are the three rational discriminant ratios. Exact square tests and the three cross-pairings leave finitely many sign choices. The induced linear map on $(x^2,x,1)$ gives the Möbius matrix, which is then tested on the entire sextic. [Covariant frames](frames.py)

The fallback has a finite height bound and also handles extra automorphisms. Work in integral degree-six charts with leading coefficients $A,B$. The scaled branch roots of the two sextics have bounds $R,S$. A matrix carrying a triple $(w_i)$ to $(z_i)$ lies in the kernel of rows $(w_i,1,-w_i z_i,-z_i)$. Multiply by $AB$ to make the entries algebraic integers. Their archimedean sizes are bounded everywhere by $H_0=\max(|A|S,|AB|,RS,|B|R)$. The four signed minors have absolute values at most $H=6H_0^3$. For a rational matrix in primitive integral coordinates, these minors equal one common algebraic integer times its entries, by Bézout's identity. Some complex conjugate of that nonzero integer has modulus at least one. Thus every primitive rational coordinate is bounded by $H$. At a common split prime with distinct finite branches, match a fixed triple against all 120 ordered source triples, then lift the six branches to modulus greater than $2H^2$. Rational reconstruction uniquely recovers each possible coordinate ratio. All retained frames satisfy the original characteristic-zero polynomial identity. [Height-certified branch matching](lift_frames.py)
