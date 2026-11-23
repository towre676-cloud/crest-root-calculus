# Chart I at a given mark: four-by-four cofactor certificates

Given a smooth ternary quartic over the rationals and a rational projective point on it, the [solver](solve.py) determines **every geometric transverse degree-three elliptic pencil at that mark**. It returns a triangular description over all coefficient fields, the conic parameter and the target equation on every surviving branch. The optional certificate gives polynomial cofactors for the six original chart-I equations. The calculation uses characteristic-zero univariate elimination and matrices of size four in place of a three-variable Gröbner lift. Its theoretical elimination degree is at most **340**, independently of the input coefficients. For a mark over any supplied number field of degree $d$, the [number-field solver](number_field.py) gives the same complete geometric support, with rational norm degrees at most **340$d$**. Coefficient heights and runtime still depend on the input. Chart II and discovery of unspecified marks are treated by the [full cubic atlas](../../applications/atlas/original/README.md).

The [new example](quartic.json) has geometric chart parameters $(s,b,c)=(0,0,0)$ at $[0:0:1]$. Two norm polynomials of degrees 142 and 198 meet in an eliminant of degree 101. Its simple linear factor gives the map. Its degree-ten factor, raised to the tenth power, contributes only excluded points where the sextic leading coefficient vanishes. The explicit original-ideal identities have targets $h^{46}s$, $h^{46}c$ and $h^{48}b$. The full geometric and rational spectra of this example and an infinite arithmetic family are proved below. [Certificate](certificate.json) · [Independent cofactor checker](proof.py) · [Evaluation of the six original cofactors](cofactors.py)

## From six equations to a free module of rank four

Choose the first permitted affine permutation of the input axes. Translate the mark to $(L,H,N)=(0,0,1)$ and use its tangent as $H=0$. The command reports the permutation and the actual affine translation and shear. Put $f(x,y)=F(s-x,-1,y)$ in this adapted frame, $a_{ij}=[x^iv^j]f(x,c+v)$, $h=a_{00}$, $u=a_{01}$, $v=a_{02}$, $w=a_{03}$. Here $w$ is a nonzero constant in the field of the mark. The transverse residual sextic is
$$
G(t,z)=t^6 f(1/t,c+b/t+z/t^2).
$$
The square shape is $G=h(S^2-K(z))$, where $S=t^3+e_2t^2+(e_1+f_1z)t+e_0+f_0z$ and $K$ is cubic. The top coefficients of $G$ determine the five shape parameters on $h\ne0$. Write $N_{ij}$ for the remaining denominator-cleared coefficient at $t^iz^j$, with $i=1,2$ and $j=0,1,2$, using the Part XXIII normalization. The [implementation](chart.py) constructs them using the displayed square. The [checker](proof.py) independently substitutes into the original quartic and compares the coefficients of the full squared numerator $(16h^3 S)^2$.

The Part XXX transverse-line osculation test gives a quartic equation for $c$ and a forced value for $b$. Here these equations are used as identities **over the polynomial ring in the unknown line slope**. Set
$$
Q=4hv-u^2,\qquad q=Q/(3w^2),\qquad D=12h^2w,\qquad B=-(4h^2 a_{12}-2hu a_{11}+u^2 a_{10}),\qquad L=Db-B.
$$
The identities $N_{22}=16h^2 Q$ and $N_{12}=16h^2(L-buQ)$ put $h^2q$ and $h^2L$ in the original ideal. The monic quartic $q$ has derivative $12wh/(3w^2)$ in $c$. Its roots on $h\ne0$ are therefore simple even on exceptional lines with repeated residual intersection elsewhere.

Work in $A=\mathbb Q[s][c]/q$, which is free of rank four without any localization in $s$. For the other four coefficients set
$$
r_{ij}=D^{6-i-2j}N_{ij}(s,B/D,c)\bmod q .
$$
These are polynomials. Indeed, in the coefficient $[t^i z^j]G$ the power of $b$ is at most $6-i-2j$, and the weight is additive under square completion. The reduced implementation performs polynomial multiplication modulo the monic $q$ at each step. A second check uses the unmodified original generators before substitution.

Every residual solution on $h\ne0$ has a separable cubic target. Absorb $s,b,c$ into the linear substitution $(X,Y,Z)\mapsto(X-sY,Y,Z-bX-cY)$ in the adapted quartic and call the result $F'$. Then $t^2G(t,z)=F'(-t,-t^2,z)$. The mark $[0:0:1]$ has tangent $Y=0$, so $G(0,z)=z^2(\alpha+\beta z)$ with $\beta\ne0$, the negative coefficient of $YZ^3$ in $F'$. Since $S(0,z)$ has degree at most one, $K$ has degree three. If $K$ had a repeated root $z_0$, each root of the monic $S(t,z_0)$ would give a singular point of $G=h(S^2-K)$. For $t\ne0$ the displayed substitution is an isomorphism onto the plane open $XY\ne0$, with inverse $t=Y/X$, $z=-YZ/X^2$. Smoothness forces $S(t,z_0)=t^3$. Since $G(0,z_0)=G_z(0,z_0)=0$, the expression for $G(0,z)$ gives $z_0=0$. Now $S=t^3+z(f_0+f_1t)$, and the coefficients of $t^2$, $t^3$ and $tz$ in $G$ vanish. Up to sign these are the coefficients of $X^4$, $X^3Y$ and $X^3Z$ in $F'$. Their vanishing makes $[1:0:0]$ singular, a contradiction. The fixed-degree elliptic-quotient finiteness theorem therefore applies to every residual solution on this open chart.

For a linear combination $r$ of the four rows, let $M_r$ be multiplication by $r$ on $1,c,c^2,c^3$. Cramer's rule gives
$$
e(s)=\det M_r,\qquad u_r=(1,c,c^2,c^3)\operatorname{adj}(M_r)(1,0,0,0)^{\mathsf T},\qquad u_r r=e(s)\pmod q.
$$
Thus each nonzero norm comes with a polynomial identity. Typically two rows suffice. If individual norms vanish, try coefficients $(1,k,k^2,k^3)$ for $0\le k\le12$. On a smooth quartic $h$ has no common component with $q$ over $\mathbb Q(s)$, since $h$ is the generically separable residual cubic of a projection in characteristic zero. A common component of $q$ and all four residuals would then give infinitely many transverse pencils at a fixed mark. Finiteness follows from the degree-three elliptic quotient theorem and the finite marked fibres for each quotient. Each of the at most four irreducible components rules out at most three scalar values. Consequently one of the thirteen combinations has nonzero norm.

Let $E$ be the monic gcd of the selected norms. Exact polynomial division checks $E\mid e_i$. A good-prime coprimality check on the primitive quotients verifies the converse over $\mathbb Q[s]$. More explicitly, for two norms write $a=e_1/E$ and $b'=e_2/E$. The first column of the inverse Sylvester matrix, whose columns in ascending coefficient order are $s^ja$ for $j<\deg b'$ followed by $s^jb'$ for $j<\deg a$, gives coefficients $A,B'$ with $Aa+B'b'=1$. Constant quotient cases require one scalar division. This supplies a characteristic-zero cofactor for $E$. The certificate also verifies the two norm cofactors by direct multiplication.

## Fibres, denominators and the original polynomial identities

For each monic irreducible $p(s)$ in $E$ work over $K_p=\mathbb Q[s]/p$. The output supplies the monic $d_p=\gcd(q,r_{20},r_{21},r_{10},r_{11})$ in $K_p[c]$ and an explicit five-entry Bézout row for $d_p$. The checker verifies that $d_p$ divides all five inputs and that the cofactor combination reproduces $d_p$. Set
$$
g_p=d_p/\gcd(d_p,h^4).
$$
The fourth power removes all excluded roots counted with multiplicity. Every remaining root is simple because $q'$ is a nonzero scalar times $h$. Thus $\deg(p)\deg_c(g_p)$ is the exact geometric count on this fibre. The conic parameter is $B/D\bmod g_p$. The solver reconstructs the square, checks the full sextic identity in $K_p[c]/g_p$ and verifies that the elliptic target cubic is separable. A singular source is rejected.

For clarity, the following formulas specify the six original polynomial cofactors as arithmetic circuits, including on a chart with no pencils. Their inputs are the explicit rational polynomials in the certificate. If any fibre survives, use the two targets $T_s=\prod_{g_p\ne1}p$ and $T_c$, whose value modulo each surviving $p$ is $g_p$. The ordinary Chinese remainder formula gives $T_c$. An empty chart uses the target $T=1$. Let $m_p$ be the multiplicity of $p$ in $E$, $M=\sum_p m_p$ and $m=\sum_{g_p\ne1}m_p$. For one of these targets, take $T_p=T$ on surviving fibres and $T_p=1$ elsewhere. In $K_p[c]$, multiply the certified Bézout row for $d_p$ by $(h^4T_p)/d_p$, an exact polynomial quotient. Lift the entries to $\mathbb Q[s,c]$ with $s$-degree less than $\deg p$. This gives
$$
X_p=h^4T_p=I_p+pW_p,\qquad I_p=A_{p0}q+\sum_j A_{pj}r_j,\qquad W_p=(X_p-I_p)/p .
$$
The quotient defining $W_p$ is exact. Multiply the identities with multiplicities $m_p$. A telescoping product expresses $\prod X_p^{m_p}-E\prod W_p^{m_p}$ in $q,r_1,\ldots,r_4$. The previous norm cofactors express $E$ in the same generators. Thus $h^{4M}T^m$ lies in their ideal; for an empty chart the target is simply $h^{4M}$.

To transport this identity back to the original six generators, set $H=h^2$, $k=6-i-2j$ and $P_j=D^k N_{ij}(s,B/D,c)$. Then
$$
Hr_j=HD^k N_{ij}-\frac{P_j-r_j}{q}Hq-\frac{D^kN_{ij}-P_j}{L}HL.
$$
Both displayed polynomial divisions are exact. Substitute $Hq=N_{22}/(48w^2)$ and $HL=(N_{12}+bu N_{22})/16$ to obtain the six original cofactors. The final target is $h^{4M+2}T^m$. For a single rational triple with $m=1$, the additional target $B-Db_0$ and the identity for $L$ yield $h^{4M+4}(b-b_0)$. This proves literal first-power coordinate membership in these examples. General projected multiplicities are retained explicitly. A large exponent on an excluded factor does not contribute spurious pencils. The [circuit evaluator](cofactors.py) computes each of the six coefficients at a rational point or at a residue point where the displayed denominators are nonzero. The proof of the identities is polynomial. Numerical evaluations are extra controls.

For a uniform degree bound assign weight one to $s$ and two to $c$. The coefficients $a_{ij}$ have weight at most $6-i-2j$. The monic elimination $q$ has weight at most eight, so reducing in $c$ never increases weight. The weights of $D,B$ are at most 12 and 13. The weights of $(r_{20},r_{21},r_{10},r_{11})$ are bounded by $(68,44,85,61)$. For a row of weight $W$, entry $(i,j)$ of its multiplication matrix has weight at most $W+2j-2i$. Every determinant term therefore has $s$-degree at most $4W$. This proves the bound 340 including fallback combinations.

## Marks over a number field

Supply a rational quartic, a monic irreducible polynomial $m(\theta)$ over $\mathbb Q$ and a projective point with three coordinates in $K=\mathbb Q[\theta]/m$. Each coordinate is a list of rational coefficients in increasing powers of $\theta$. The [quadratic example](quadratic.json) encodes $[0:\sqrt2:1]$. The frame, rank-four adjugates and four residuals are computed over $K$. This treats rational curves with algebraic marks of any degree and permits conic and map parameters in further extensions.

For a nonzero mark-field norm $e_i(s)$ from the rank-four calculation, replace $s$ by $U-j\theta$, with $j$ initially zero. Take the resultant with $m$ in $\theta$ to obtain a rational polynomial $n_i(U)$. Its degree is at most $340[K:\mathbb Q]$. Factor the gcd of the $n_i$ over $\mathbb Q$. For each irreducible factor $p(U)$, compute
$$
\gcd\bigl(m(\theta),e_1(U-j\theta)\bigr)
\quad\text{in}\quad (\mathbb Q[U]/p)[\theta].
$$
A linear gcd recovers $\theta$ and $s=U-j\theta$ in the same number field. If a gcd has higher degree, increase $j$ and restart. There are finitely many geometric pairs $(\theta,s)$ with $m(\theta)=e_1(s)=0$. Distinct pairs with different $\theta$ have equal $s+j\theta$ for at most one $j$. Thus a separating integer exists and this loop terminates, including cases with repeated roots. If a candidate arises from different conjugates of the separate norms, the simultaneous fibre gcd of the five chart equations still tests it exactly.

At each resulting factor, the routine checks the full square identity, smoothness of its elliptic target and the five-term fibre Bézout identity over $\mathbb Q[U]/p$. The surviving $c$-polynomial is squarefree. Its contribution per embedding of the originally prescribed mark is $\deg(p)\deg_c(g_p)/[K:\mathbb Q]$. Galois acts transitively on these embeddings, so division yields an integer. The exported [quadratic certificate](quadratic-certificate.json) includes the mark-field norm cofactors and each specialized fibre cofactor. The arbitrary-degree routine provides this triangular certificate. The standalone six-original-generator circuit evaluator currently takes rational marks.

| Prescribed input | Rational mark cases | Algebraic mark cases | Geometric outcome |
|---|---:|---:|---|
| New degree-three family and coordinate transforms | 9 | 0 | One at each chosen mark |
| Historical $C^\star$ | 3 | 0 | Two positive marks and one negative |
| Cyclic quotient | 1 | 1 | One at each chosen mark |
| Diagonal, perturbed and catalogue controls | 21 | 0 | Zero at each chosen mark |
| Fermat, Klein and a second Klein-family curve | 0 | 3 | One per embedding at the displayed quadratic mark |
| Asymmetric conjugate-fibre curve | 0 | 1 | One per embedding at the displayed quadratic mark |
| $X^4+Y^3Z-2Z^4$ | 0 | 1 | Zero at $[0:\sqrt[3]2:1]$ |

The [five additional algebraic inputs](algebraic-marks.json) include a cubic mark with rational norm degrees 516 and 648. The reduced geometric chart and its maps are determined for every qualifying marked quartic. The powers in original-ideal membership certificates can be nonminimal. The calculation does not identify all nilpotent thicknesses of the saturated chart, or find the unspecified marks on the source.

## An infinite family with geometric degrees 3n and rational degrees 3n²

Put
$$
S=t^3+t^2+tz+1+z,\qquad K_k(z)=z^3+kz^2+2z+1 .
$$
Define the plane quartic by $t^{-2}F_k(-t,-t^2,z)=S^2-K_k$. Explicitly, each term $A t^r z^n$ on the right gives the monomial $(-1)^{i+j} A X^iY^jZ^n$ with $j=r+n-2$ and $i=4-j-n$. All exponents are nonnegative. The point $[0:0:1]$ is rational. The identity recovers a degree-three map to $E_k:w^2=K_k(z)$ over $\mathbb Q$. A good reduction at 5 proves smoothness for every integer $k\equiv1\pmod5$. Its counts over the first three extensions are 5, 49 and 131. Its reciprocal characteristic polynomial is
$$
(T^2+5)(T^4-T^3+7T^2-5T+25).
$$
The existing [root-ratio and specialization calculation](../elliptic-quotients/polarized-reduction/README.md) finds exactly one geometric elliptic block of multiplicity one over the algebraic closure of the residue field. Hence the global Jacobian has one elliptic direction and its complementary surface is geometrically simple. The exhibited map, of prime degree, is optimal. Every geometric elliptic map is its composite with an isogeny, so the degrees are precisely $3n$ for positive integers $n$. Up to translation and target isomorphism there are $\sigma_1(n)$ classes in degree $3n$, one for each subgroup of order $n$ in the elliptic curve.

For $k\equiv6\pmod{5005}$, the good-prime traces of $E_k$ at $5,7,11,13$ are $0,-3,-4,-6$. For each allowed odd rational prime-isogeny degree $\ell$ there is a different residue prime $p$ in this list with $a_p^2-4p$ nonsquare modulo $\ell$. At 7 the odd trace excludes a rational order-two subgroup. These tests exhaust [Mazur's possible prime degrees](https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf). A rational finite subgroup with elementary divisors $m\mid n$ contains all of $E_k[m]$. After division by this subgroup its cyclic remainder would have a rational subgroup of prime order unless $m=n$. Thus all rational isogeny kernels are exactly the $E_k[m]$. The rational source point permits normalization by translation, and there is exactly one rational class in each degree $3m^2$. The invariant $j=\big(16(k^2-6)\big)^3/(16\operatorname{disc}K_k)$ is nonconstant. The unique optimal target therefore distinguishes infinitely many geometric source curves in either progression. [Family and sieve verifier](family.py)

## Replay

With Python 3.11+ and SymPy and NumPy installed, run from this directory.

```sh
python solve.py quartic.json 0 0 1
python solve.py quartic.json 0 0 1 --cofactors --output my-proof.json
python cofactors.py quartic.json my-proof.json 0 0 1 --at 11 5 7 --prime 1009
python number_field.py quadratic.json --cofactors --output mark-proof.json
python number_field.py quadratic.json --verify quadratic-certificate.json
python verify.py --workers 3
python verify.py --panel --workers 3
```

A coefficient entry `[i,j,a]` denotes $a c^i s^j$. Lists for univariate factors and for target cubics use increasing powers. The certificate is relative to the primitive integral equation in the reported adapted frame. The six-generator order in the evaluator is $N_{20},N_{21},N_{22},N_{10},N_{11},N_{12}$. Omitting `--prime` evaluates the characteristic-zero rational coefficients. The [exact evaluator](bezout.py) uses modular Bézout updates and Chinese remaindering with an integral Hadamard bound to reconstruct the required Sylvester coefficients at the chosen rational point. The default verifier checks six rational cases, a stored quadratic certificate, a new cubic computation, the family reduction and isogeny tests, a case in which all four individual norms vanish, and five deliberately invalid inputs or certificates. The extended panel covers all 40 cases in the table. Timestamps and runtimes depend on machine and concurrency.
