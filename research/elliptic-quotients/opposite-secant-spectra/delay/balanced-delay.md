# Finite quotient delay beyond the third successive-minimum length

Write \(\lambda_1(G)\le\lambda_2(G)\le\lambda_3(G)\) for the **squared** successive minima of a positive integral ternary lattice \(G\). Write \(\delta(L_+,L_-)\) for the first degree at which the geometric optimal elliptic quotient counts on the two realized children differ.

**Theorem.** There is an effective infinite family of opposite-secant pairs of genus-two Jacobians of CM product type with finite first difference, satisfying
\[
 \frac{\lambda_3(G)}{\lambda_1(G)}\longrightarrow1,
 \qquad
 \delta(L_+,L_-)\gg\lambda_3(G)^{4/7}.
\]
In particular \(\delta/\sqrt{\lambda_3(G)}\) is unbounded even for parents with asymptotically equal successive minima. The counts agree in every degree prime to three. The construction gives their first common degree \(n\), their first differing degree \(ns\), and counts **two versus zero** there. The ratio \(\delta/\sqrt{\lambda_3(G)}\) is asymptotic to \(s\), where \(s\) runs through \(2k^2+1\) for \(k\ge2\), \(k\equiv2\pmod3\).

## Construction and proof

Put \(Y=2k\), \(Z=2k^2\), \(s=Z+1\). Then \(Y^2=2(s-1)\), and the vector \(v=(1,Y,Z)\) has Euclidean squared length \(s^2\). Also \(Y\equiv1\) and \(s\equiv0\pmod3\). For \(w=(x,y,z)\) use the integral coordinates
\[
 a=x,\qquad b=y-Yx,\qquad c=z-Zx,
\]
and define the integral quadratic perturbation
\[
 B(w)=\alpha b^2+6(ab+ac+bc),\qquad
 \alpha=6(3^e+1),
\]
where \(3^e\) is the least power of three strictly exceeding \(11s^2\). In the adopted coordinates the Gram matrix has off-diagonal entries three and diagonal \((0,\alpha,0)\), with determinant \(-54\cdot3^e\). The coordinate change has determinant one. In particular the determinant in the original coordinates has only the prime divisors two and three. All coefficients of the Gram matrix are divisible by three. Since \(Y,Z\) are even and \(4\mid\alpha\),
\[
 B(w)\equiv2(xy+xz+yz)\pmod4 .
\]

The perturbation vanishes on the axis \(\mathbb Ze_3\) and on \(\mathbb Zv\). These are its only zeros of Euclidean length at most \(s\). Indeed, when \(b=0\) it equals \(6ac\), giving exactly the two stated lines. If \(b\ne0\) and \(\|w\|\le s\), then
\[
 \frac{|a|}{|b|}\le\frac{s+1}{Y}\le2\sqrt s,\qquad
 \frac{|c|}{|b|}\le s+(s-1)\frac{s+1}{Y}\le3s^{3/2}.
\]
Thus \(|ab+ac+bc|/b^2\le11s^2<\alpha/6\), so \(B(w)>0\). For all \(w\) the estimates \(|a|\le\|w\|\), \(|b|\le(1+Y)\|w\|\), and \(|c|\le s\|w\|\) give
\[
 |B(w)|\le W\|w\|^2,\qquad
 W=\alpha(1+Y)^2+6(1+Y+s+(1+Y)s)=O(s^3).
\]

**Choosing the scale.** Set \(J=s^2\) and let \(N\) be the first integer congruent to one modulo six strictly greater than
\[
 2Ws(12Js+1).
\]
We choose \(n\) among \(N,N+6,\ldots,N+6J\) so that for every nonsquare integer \(1\le R\le s^2\), the integer \(n^2R\) has distance strictly greater than \(WR\) from every square. For a given nonsquare \(R\), at most one candidate fails. To see this, a failure yields
\[
 \|n\sqrt R\|_{\mathbb R/\mathbb Z}\le \frac{WR}{n\sqrt R}\le\frac{Ws}{N}.
\]
Two failures would give a difference \(h\) with \(1\le h\le6J\) and \(\|h\sqrt R\|_{\mathbb R/\mathbb Z}<1/(12Js+1)\). If \(m\) is the integer nearest to \(h\sqrt R\), then \(|h^2R-m^2|\ge1\) since \(R\) is nonsquare. Hence
\[
 |h\sqrt R-m|\ge \frac{1}{h\sqrt R+m}\ge\frac{1}{12Js+1},
\]
a contradiction. There are at most \(J\) forbidden values of \(R\) and \(J+1\) candidates. This gives an effective choice with \(n\ll s^7\). The finite tests use the two squares adjacent to \(n^2R\).

**The spectra.** Take
\[
 Q(w)=n^2\|w\|^2+B(w),\qquad
 L_+=\{x-y\equiv0\pmod3\},\qquad
 L_-=\{x+y\equiv0\pmod3\}
\]
inside the parent \(G=(\mathbb Z^3,Q)\). The parent Gram matrix reduces to \(I\) modulo three. The two planes cut complementary pairs of points on its isotropic conic, and among projective points of square unit norm they both contain exactly the axis of \(e_3\). Consequently the ordinary square counts agree in every degree prime to three, and divisor inversion gives the same assertion for primitive counts.

If \(Q(w)\le n^2s^2\), the norm bound forces \(\|w\|\le s\). In fact \(\|w\|^2\ge s^2+1\) implies
\[
 Q(w)\ge(n^2-W)(s^2+1)>n^2s^2
\]
by the choice of \(n\). Now set \(R=\|w\|^2\). If \(R\) is nonsquare, the choice of \(n\) excludes a square value of \(Q(w)\). If \(R=r^2>0\) and \(Q(w)=m^2\), then \(|m-nr|\le Wr/n<1\), so \(m=nr\) and \(B(w)=0\). Thus the only nonzero vectors of square norm at or below \(n^2s^2\) are \(\pm je_3\) for \(1\le j\le s\), together with \(\pm v\). On the axis only \(\pm e_3\) are primitive. The vectors \(\pm v\) are primitive since the first coordinate is a unit, and they belong only to \(L_+\) since \(Y\equiv1\pmod3\). Thus the only primitive counts below \(ns\) are two in degree \(n\) on each side. At \(ns\) the counts are two and zero. The finite difference also proves that the two lattices are inequivalent.

**Jacobian realization.** The congruences give \(Q(w)\equiv(x+y+z)^2\pmod4\), and \(Q(w)\ge(n^2-W)\|w\|^2\) excludes norm one. Both children represent \(n^2\) on \(e_3\). Their Gram determinant is nine times the determinant of \(n^2I+\mathrm{Gram}(B)\), which modulo any prime dividing \(n\) is nine times \(-54\cdot3^e\). Because \(n\equiv1\pmod6\), the root \(n\) is prime to the discriminant of either child. This also gives primitivity, since an odd prime dividing the form's content would divide both its determinant and \(n\), and the represented square is odd. Kani's realization criterion now supplies principally polarized abelian surfaces of CM product type over \(\overline{\mathbb Q}\). With norm one absent their polarizations are indecomposable, so they are Jacobians of smooth genus-two curves. The refined Humbert correspondence identifies the primitive counts with geometric optimal elliptic quotient classes. [CM realization criterion](https://www.queensu.ca/academia/kani/sites/acadkaniwww/files/uploaded_files/Preprints/refHum5.pdf) · [Subcover correspondence](https://doi.org/10.1007/s40316-018-0105-6)

Finally
\[
 n^2-W\le\lambda_1(G)\le\lambda_3(G)\le n^2+W,
\]
where the upper bound uses the coordinate basis. Since \(n>Ws\), all three squared minima are asymptotic to \(n^2\) and \(\delta/\sqrt{\lambda_3(G)}\sim s\). The effective choice \(n\ll s^7\) gives \(\delta=ns\gg n^{8/7}\gg\lambda_3(G)^{4/7}\), as claimed.

## Scope and replay

The determinant of the children is asymptotic to \(9n^6\), since the parent perturbation has size \(O(W)=o(n^2)\). The exponent of the maximal finite delay as a function of the **child determinant** remains open. The separate [elongated family](../README.md#an-infinite-family-with-an-unbounded-first-difference) still gives the stronger known determinant-scale lower bound of order the square root of the determinant. The construction here settles the proposed successive-minimum bound by a different, nearly round family.

The [receipt](balanced-delay.json) gives \(s=9,51,129,801\). It contains the integer parent Gram matrices, coprime witness roots and first distinguishing degrees. For the first three members the replay independently enumerates every vector inside the Euclidean cutoff and compares the two signed square counts. The last member uses the exact adjacent-square test and the proved description of the perturbation zeros.

~~~sh
study=research/elliptic-quotients/opposite-secant-spectra/delay
python -B "$study/balanced_delay.py" --verify "$study/balanced-delay.json"
~~~
