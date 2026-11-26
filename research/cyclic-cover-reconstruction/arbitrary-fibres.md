# Reconstructing cyclic quotients with arbitrary infinity fibres

An exact search now recovers horizontal poles away from the source branch locus and allows the source point at infinity to map to a finite target branch point. It also handles mixed target branching. The input is one rational cyclic-cover equation, with the same power removal, character and source-infinity hypotheses as in the [cyclic theorem](proof.md). The target equation, quotient coordinates and pole locations are unknown. Give a total degree cap when elliptic targets are included. For higher-genus targets Riemann–Hurwitz supplies a cap from the source alone.

## Scope and normalizations

In a character chart write \(y^k=B(x)\), with finite multiplicities \(0,\ldots,k-1\), and \(\gcd(\deg B,k)=1\). The original input may have larger exponent and removable powers. The character chart and its index account for those before the search. The target is \(z^k=G(u)\), with
\[
G(u)=c\prod_{\nu=1}^{k-1}G_\nu(u)^\nu,\qquad
m=\deg G,\qquad \gcd(m,k)=1.
\]
The \(G_\nu\) are monic, squarefree and pairwise coprime. Their degree profile determines the genus. This allows a mixed target with a rational fully ramified place at infinity. In particular, an equivariant map sends the fully ramified rational source infinity to a rational full fixed point of the target action. If the original target infinity is unramified or partly ramified, send that fixed point to infinity by a rational change of \(u\) and remove powers from its radicand.

Write \(u=A/Q\) in coprime polynomials. Normalize \(Q\) to be monic and scale \(u\) so that \(A\) is monic. If source infinity maps to target infinity, let \(a=\deg A>q=\deg Q\) and use translation to set the coefficient of \(x^q\) in \(A\) to zero. Otherwise its image is a finite root of \(G\), which can be translated to \(0\). In that chart \(a<q\) and the root multiplicity \(\nu_0\) must be coprime to \(k\). Both normalizations are over \(\mathbb Q\) for a rational map. Normalize the numerator of the vertical multiplier to be monic by scaling \(z\). The working map has degree \(d=\max(a,q)\), and the original degree multiplies by the character-chart index.

The implementation lists normalized maps in these marked target charts. Changing the marked target infinity can give more than one normalized model of the same underlying map.

## The two types of finite poles

Valuations at a pole of horizontal order \(b\) give \(mb+r=0\bmod k\), where \(r\) is its multiplicity in \(B\). At branch places, \(r\ne0\), so the denominator part \(Q_0\) is specified by known rational factors and finitely many admissible orders. Elsewhere \(k\mid b\). Thus
\[
Q=Q_0H^k,\qquad D=D_0H^m,\qquad
D_0^k=Q_0^m R,\qquad F=B/R.
\]
Here \(H\) is an unknown monic polynomial coprime to \(B\), \(R\) removes the selected branch factors with their multiplicities, and \(D\) is the vertical denominator. The search bounds \(\deg H\) by \((d-\deg Q_0)/k\) and recovers its coefficients from the incidence. Repeated roots of \(H\) are permitted. Away from the poles, regularity forces the vertical multiplier to have no denominator because every finite multiplicity of \(B\) is less than \(k\). The exact identity is
\[
 G_h(A,Q)=FV^k,\qquad G_h(A,Q)=Q^mG(A/Q),
\]
and the vertical coordinate is \(z=yV/D\) in the working chart.

## Derivative reduction for mixed branching

Set \(P=\lfloor(k-1)/\max\{\nu:G_\nu\ne1\}\rfloor\). If a finite root of \(V\) has order \(v>0\) and its horizontal image is a target root of multiplicity \(\nu\), its horizontal ramification index \(h_x\) satisfies
\[
\nu h_x=r+kv> \nu Pv ,
\]
so the Wronskian has order \(h_x-1\ge Pv\) there. The same inequality covers roots of \(V\) that are also roots of \(F\). At a pole of order \(b\) the Wronskian has order \(b-1\). With
\[
 T_0=\prod_{p\mid Q_0}p^{\operatorname{ord}_p(Q_0)-1},
\]
the resulting factorization is
\[
 A'Q-AQ'=(a-q)\,T_0 H^{k-1}V^P W
\]
for a monic \(W\). At a root of \(H\) of multiplicity \(b'\) this uses only the guaranteed contribution \((k-1)b'\le kb'-1\), leaving the rest to \(W\).

Let \(r_\ast=m\) in the infinity-pole chart and \(r_\ast=\nu_0\) in the finite-root chart. The highest-degree term of the homogenized target is the one with \(A^{r_\ast}Q^{m-r_\ast}\), so
\[
 k\deg V+\deg F=mq+r_\ast(a-q),\qquad
 \deg W=a+q-1-\deg(T_0H^{k-1})-P\deg V.
\]
Both are exact stratum constraints. A further local bound distributes the branch weights of \(F\) among all finite target roots. A root of multiplicity \(\nu\) receives weight at least the residue of \(\nu d\bmod k\), since the fibre polynomial has degree \(d\). For the root at \(u=0\) in the finite-root chart, use \(a\) in place of \(d\).

Reverse the polynomials in \(t=1/x\) and write \(\alpha,\gamma,v,\tau,w,f\) for the reversals of \(A,Q,V,T_0H^{k-1},W,F/\operatorname{lc}(F)\). Put \(\ell=a-q\). Then
\[
 \alpha=(fv^k\gamma^{r_\ast-m})^{1/r_\ast}\pmod {t^{|\ell|}},
 \qquad
 \ell\alpha\gamma-t\alpha'\gamma+t\alpha\gamma'
   =\ell\tau v^P w .
\]
In the complementary search the coefficient of \(v_j\) can be removed at every \(1\le j<|\ell|\) up to the degree of \(V\). If \(C_j\) is the coefficient of the difference of the Wronskian sides with \(v_j=0\), the recovered value is
\[
 v_j=\frac{r_\ast C_j}{r_\ast\ell P-k(\ell-j)} .
\]
There are **no resonant denominators** in admissible strata. From the pole congruence \(mq+\deg R=0\bmod k\) and the degree identity above, \(r_\ast\ell=\deg B\bmod k\). Therefore the denominator reduces modulo \(k\) to \(P\deg B\), which is nonzero since \(1\le P<k\) and \(\deg B\) is coprime to \(k\). The formula holds for negative \(\ell\) as well as positive \(\ell\). Squarefree targets give \(P=k-1\), recovering the existing pole-chart recurrence.

The complementary search solves the remaining linear Wronskian equation for \(A\). The alternative starts from the leading formal root and leaves the missing lower coefficients free. The solver picks the smaller coefficient list. It then reads the target coefficients by descending degrees of \(A^iQ^{m-i}\). In the finite-root chart the descending order starts with \(i=\nu_0\) and increases. For mixed targets the repeated root factors are unknown monic polynomials, and division by their product recovers the simple factor. Discriminants and pairwise resultants enforce the precise branching profile. Coprimality for \(A,Q\) and \(H,B\) is imposed by resultants. This exact incidence, with open conditions imposed by inversion, is the certificate input.

## Completeness and the degree cap

For any solution the two original equation identities are checked by polynomial arithmetic in the working chart, and the character calculation transports them to the original equation. The coprimality and degree conditions give a nonconstant map of the recorded degree. Conversely any equivariant map in the stated scope yields one of the profiles, charts and pole divisors above. Its coefficients satisfy the exact identities and every derivative reduction, so its normalized coefficients appear in the appropriate incidence. Positive-genus targets and bounded degree give finitely many maps up to pointed target isomorphism, by the [bounded pointed argument](../coefficient-quotients/proof.md#finiteness-with-an-unknown-elliptic-target). Each target action has finitely many full fixed points and only finitely many such marked models after normalization. The saturated incidences are therefore zero-dimensional. The solver enumerates the rational points, and its verifier also checks any empty degree strata.

When the target genus is at least two, the total degree is at most \(g(C)-1\), with sharper per-profile bounds. A degree cap is genuinely needed when elliptic targets are allowed. Multiplication by \(N\) on an elliptic target preserves equivariance and gives infinitely many degrees \(N^2d\). It can introduce poles away from the branch locus.

**Hyperelliptic consequence.** For odd-degree hyperelliptic source models, any morphism to a hyperelliptic curve of genus at least two respects the involutions. If \(P\) is the source Weierstrass point and \(f(P)=Q\), pushforward gives \(f(R)+f(\iota R)\sim 2Q\). The divisors move, so they are fibres of the target's unique degree-two pencil and \(Q\) is Weierstrass. The target therefore has a rational odd-degree model and the exponent-two search needs no degree cap in higher genus. For genus one, take \(Q\) as the origin, and apply the bounded exponent-two search. This yields a source-only search for arbitrary-fibre hyperelliptic quotients in that source scope.

**An optimal rational map that breaks the cyclic symmetry.** Consider the cubic cover
\[
 y^3=(1-2x)(x^2-1)^2.
\]
Set \(t=(1-2x)/(x^2-1)\), \(\xi=y/(x^2-1)\), and \(\eta=1+xt\). Then
\(\xi^3=t\), \(\eta^2=\xi^6+\xi^3+1\), and \(x=(\eta-1)/t\). The following map is defined over \(\mathbb Q\).
\[
 U=\xi+\xi^{-1},\qquad V=\eta(\xi+1)/\xi^2,\qquad
 E:\ V^2=(U+2)(U^3-3U+1).
\]
The coordinates are invariant under \((\xi,\eta)\mapsto(\xi^{-1},\eta/\xi^3)\). Over their function field, \(\xi\) satisfies \(\xi^2-U\xi+1=0\) and \(\eta\) is recovered from \(V,\xi\). Thus the map has degree two and is optimal, since a lower-degree elliptic map would have degree one. Moving the root \(U=-2\) to infinity gives the cubic \(Z^2=w^3+9w^2+6w+1\), with \(j=790272\). The order-three source action fixes source infinity. Equivariance would require an order-three target automorphism fixing its image, which this \(j\)-invariant excludes. A trivial target action would force the map to factor through the rational source quotient. This gives a rational, optimal counterexample to unconditional equivariance already in degree two. The replay checks the curve and map identities. General non-equivariant discovery remains open.

## Exact replay

The replay recovers \(A=x^4+x+2,\ Q=x^3,\ G=u^2-1\) from \(B=A^2-Q^2\), with the unbranched pole inferred from the source. It recovers \(u=1/[x(x-1)^3]\) from \(B=x(1+x(x-1)^3)\), a chart with source infinity over a finite branch root and with both kinds of finite pole. Further checks cover a repeated unbranched pole, mixed target powers, a nonconstant multiplier at a double target root, a lifted exponent-six source, an exhaustive higher-genus search with no supplied cap, exact exclusion on a perturbed source, and rejection of an incomplete receipt. Every exported map used in these checks is substituted in the original coordinates. The degree-four unbranched example has an independent Gröbner calculation of the raw coefficient identity. A degree-eight doubling example checks the complementary coefficient step at a repeated unbranched pole, where the formula reduces the coefficient list from twelve to six, including the unknown pole coefficients.

~~~sh
python -B research/cyclic-cover-reconstruction/verify_fibres.py
python -B research/cyclic-cover-reconstruction/fibres.py source.json --max-degree 4
python -B research/cyclic-cover-reconstruction/fibres.py source.json --higher-genus-only --mixed-targets
~~~

Two exponent-four strata check partial target ramification. The first recovers \(u=x^3+x\) with \(G(u)=(u-1)^2(u+1)\). The second recovers \(u=x/(x-1)^4\) with \(G(u)=u^3(u-1)^2\) from \(B=x^3\bigl(x-(x-1)^4\bigr)^2\). In the latter, the solver discovers \(x=1\) as an unbranched pole, and source infinity maps to the triple root at \(u=0\). Both certificates replay the original equation identities.

The JSON input has the keys **exponent** and **coefficients** in ascending order, as in [the earlier example](branch-poles.json). **--mixed-targets** enables mixed profiles, and **--verify** checks an existing report. The specialized [branch-supported solver](branch-poles.md) retains its faster fixed-denominator route. This general incidence proves completeness in its scope. Its running time can grow substantially with the cap and the number of profiles.
