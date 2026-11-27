# The octic cyclotomic family

Let $k=\mathbb Q(w)$ with $w^4-w^2+1=0$, so $w$ is a primitive twelfth root of unity. Write $L_d=k(\sqrt d)$ for squarefree $d>1$ coprime to $3$. This choice gives one representative per field. The curve $C:y^2=x^6-1$ has exactly ten points over $k$, including its two infinities. Among the 91 independent inputs $d\le200$ there are 49 for which the CM quotient has certified rank at most two. All 49 are resolved. Forty-eight have exactly the base ten points, and $L_2$ has exactly 46. Of the 49 fields, 23 have rank zero and 26 have rank two. The other 42 have certified rank at least four and carry partial status in the panel.

## Abscissas outside the quartic field

The Jacobian has rank zero over $k$ by its two elliptic factors. Good reduction at 13 and 37 bounds its rational torsion by 144. In the odd chart $z=1/(x-1)$, $v=z^3y$, the verified Cantor addition in octic_divisors.py generates all 144 reduced divisor classes using $(z,1),(z+1,w^3)$ and five Weierstrass points. A pair of conjugate points whose abscissa lies in $L_d\setminus k$ gives a $k$-rational divisor of degree two outside the hyperelliptic pencil, hence is its unique effective representative. It is fully specified by the reduced Mumford form $(u,v)$ with $u$ irreducible quadratic. There are 96 such divisors over $k$.

If $\delta=\mathrm{disc}(u)$, the pair splits in $L_d$ precisely when $\delta/d$ is a square in $k$. Only 2 and 3 ramify in $k/\mathbb Q$. Thus at each rational prime $q>3$, any prime over $q$ with odd valuation of $\delta$ forces $q\mid d$. Conversely $q\mid d$ outside 2 and 3 forces odd valuation at *every* prime above $q$. The product $d_0$ of rational primes detected in the numerator or denominator of this principal ideal leaves only $d_0$ and $2d_0$ to test. The exact field square tests give 12 pairs for $d=2$, 36 for $d=7$, 24 for $d=11$, and 24 pairs for which no rational radicand works. These counts hold for every $d$ without a size bound. Each pair gives precisely two affine points.

## Abscissas in the quartic field

If $x\in k$, then $y$ is in $k$ or $y\sqrt d\in k$. Write $v=y\sqrt d$ on the latter branch; then $v^2=d(x^6-1)$. The map $(x,v)\mapsto (d x^2,dv)$ to $E_d:Y^2=X^3-d^3$ has readily computable finite fibres. The eight affine base points on $C$ over $k$ are returned separately. Each twist point with $v\ne0$ contributes a new point with ordinate $v/\sqrt d$. The field for these coordinates is represented as $L_d=k[t]/(t^2-d)$.

The quotient $E:y^2=X^3-1$ has rank zero over $k$. The rank over $L_d$ equals the sum of the rational twist ranks at $d,-d,3d,-3d$. For the rank-two inputs, the characters $-d,3d$ have rank one each; their generators have zero height pairing by Galois invariance. The implementation supplies the second rational direction via the rational degree-three isogeny $(X,Y)\mapsto(X+4B/X^2,Y(1-8B/X^3))$ on $Y^2=X^3+B$. This works even when rational rank bounds are exact and the rank calculation returns no point. Its image is divided if necessary. The usual trace projectors for $\operatorname{Gal}(k/\mathbb Q)$ show that the lattice spanned by rational character generators contains four times the full free group of $E_d(k)$. Projecting any point by a character yields a point *exactly* in the indicated character subgroup, to which the rational character model applies. Thus testing halves after adding all torsion suffices to saturate the lattice. All 26 inputs have index two.

Write $\mathcal H=2\hat h$ for the doubled canonical height. Linear coordinate changes identify $E_t$ with $E$ over a quadratic extension by taking $X$ to $X/t$. Use the height bounds from the main note with $A=(\log132)/3,\ B=(\log9)/3$:
\[
        h(X/t)-A\le \mathcal H(P)\le h(X/t)+B .
\]
For each rational rank-one character in the panel we certify that every nontorsion rational point has $\mathcal H>1/4$. A counterexample would have ordinary multiplicative height of $X/t$ less than $e^{1/4}132^{1/3}<7$. Enumerate that rational window, omit torsion, and for every remaining candidate check that its fourth multiple has multiplicative height $>330$ in the scaled abscissa. Since $e^4<55$ and $132^{1/3}<6$, this forces $16\mathcal H(P)>4$.

Let the candidate rational generator be $G$, and let $M$ be the multiplicative height of $X(G)/t$. If the index of $G$ has prime divisor $\ell$, then $\ell^2 < 4(\log M+B) < 4(\mathrm{bitlength}(M)+1)$. At all candidate primes test for exact divisibility (including the rational torsion shift for $\ell=2$). A negative certificate can be faster by a good reduction for which $[\#E(\mathbb F_q)/\ell]G$ is nonzero. Divide and repeat when possible. This certifies each primitive rational generator without presupposing analytic generators or a saturated database entry.

The CM height bound for $\operatorname{rank}E(L_d)\le2$ in the main note gives $H(x)<5$, and our lattice enumeration uses the weaker bound $H(x)<14$. Thus $\mathcal H(P)<2\log14+B<7$. For each primitive character generator $G_j$, use the maximum of $1/4$ and an exact lower bound from its fourth multiple, taking $(2/3)(\mathrm{bitlength}(M_4)-1)-2\le \log(M_4)-A$. Orthogonality bounds the coefficients of $P$ in a finite ellipse in this character basis. Invert the exact rational basis change, enumerate all lattice points in this ellipse, add the full torsion, and invert the quadratic fibre. This gives exactly six zero-ordinate twist points in 48 cases and gives 12 further nonzero-ordinate points in the case $d=2$. The 12 nonzero twist points, 24 points with abscissa outside $k$, and ten base points yield 46 for $L_2$.

## A uniform upper limit on exceptional radicands

Put $d'=d/\gcd(d,2)$. Assume rank at most two and a point above $x\in k$ with $y\notin k$. For every rational $q\mid d'$, necessarily $q\ge5$. At every prime $\mathcal P\mid q$, the equation $v^2=d(x^6-1)$ forces $\operatorname{ord}_{\mathcal P}(x^6-1)$ odd. It is positive, since a negative valuation of $x$ would give the even value $6\operatorname{ord}_{\mathcal P}(x)$. All primes above $q$ are unramified. Therefore $(d')^4$ divides the norm of the numerator ideal of $x^6-1$. For nonzero $\alpha\in k$ this norm is at most $H(\alpha)^4$, by the product formula applied to $1/\alpha$. Hence
\[
 d'\ \le\ H(x^6-1)\ \le\ 2H(x)^6\ <\ 2\cdot5^6=31250.
\]
Combining this with the exact three exceptions for abscissas outside $k$ proves that every rank-at-most-two $L_d$ with $d'\ge31250$ has precisely the ten base points.

## Fields with higher rank

In $L_7=\mathbb Q(i,\sqrt3,\sqrt7)$ we certify at least 130 points, namely 10 base points, 48 more with abscissa in $k$ and nonbase ordinate, and 72 with abscissa outside $k$. The first 48 additional points arise by direct substitution in the twist; the archived coordinates have second entry $v$, corresponding to $y=v/\sqrt7$. In $L_{11}$ there are at least 58 points, from the base ten and 48 points with abscissa outside $k$. Both total lists are partial. The octic $L_7$ properly contains the separately studied degree-four field $\mathbb Q(i,\sqrt7)$.
