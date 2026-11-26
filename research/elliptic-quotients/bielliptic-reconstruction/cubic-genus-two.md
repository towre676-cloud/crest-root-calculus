# Elliptic quotient spectra of genus-two cyclic cubics

A normalized rational cyclic cubic of genus two has a complete geometric optimal elliptic degree law whenever its elliptic factors are non-CM. Write \(M(n)\) for its number of geometric quotient classes of degree \(n\), up to translation and isomorphism of the target. Then
\[
M(n)=
\begin{cases}
6\cdot 2^{\omega(m)},&
 n=2m,\ \text{all prime divisors of }m\text{ are }1\bmod3,\\
0,&\text{otherwise}.
\end{cases}
\]
Here \(m=1\) is allowed and \(\omega(m)\) counts distinct prime divisors. In particular the first degrees are 2, 14, 26 and 38, with multiplicities 6, 12, 12 and 12. Every elliptic map in this case breaks the order-three action. This identifies an entire family outside the equivariant search, including optimal maps.

## From the cyclic equation to a hyperelliptic model

Start with \(y^3=B(x)\) over \(\mathbb Q\). Remove cube powers and assume the resulting degree is coprime to 3, so source infinity is fully ramified. In genus two there are four full branch points. The sum of their weights is \(0\bmod3\), so two have weight one and two weight two. Replace the character by its inverse if needed so that the weight of infinity is one. The finite radicand is now \(A(x)C(x)^2\) with \(A=a_1x+a_0\), \(C=c_2x^2+c_1x+c_0\), and the factors are disjoint and squarefree. Divide the normalized \(y\) by \(C(x)\) and write the result as \(\xi\). Then \(t=\xi^3=A(x)/C(x)\), and the quadratic equation \(tC(x)-A(x)=0\) has discriminant
\[
 \eta^2=P\xi^6+Q\xi^3+R,\qquad
 \eta=2tc_2x+c_1t-a_1,\quad
 P=c_1^2-4c_2c_0,\quad Q=4c_2a_0-2c_1a_1,\quad R=a_1^2 .
\]
Solving the displayed linear formula for \(x\) gives the inverse on function fields. The program checks the cubic identity in the original coordinates, including the removed powers and inverse-character cases.

Geometrically rescale this sextic to \(Y^2=X^6+\lambda X^3+1\), with \(\lambda=Q/\sqrt{PR}\). In this model the complementary degree-two quotients have
\[
 U=X+X^{-1},\quad V_\pm=Y(X\pm1)/X^2,\quad
 E_\pm:\ V_\pm^2=(U\pm2)(U^3-3U+\lambda).
\]
Their \(j\)-invariants are \(j(\lambda)\) and \(j(-\lambda)\), where \(j(\lambda)\) is the invariant of the cubic with coefficients \(p=-6\lambda-15\), \(q=\lambda^2+14\lambda+22\) in the model \(z^2=w^3+pw+q\). Thus the coefficient-only solver only needs a quadratic field for the \(j\)-check, even when the marked involutions need larger fields.

## The isogeny and the two-torsion labels

On the Jacobian let \(\rho\) be the order-three rotation of \(X\), and let \(\iota\) be the lift of \(X\mapsto1/X\) preserving \(V_+\). Since the cubic quotient has genus zero, \(1+\rho+\rho^2=0\) on the Jacobian. The involution satisfies \(\iota\rho\iota=\rho^{-1}\). For the two elliptic adjoint inclusions \(i_\pm\) and projections \(q_\pm\), one has \(i_\pm q_\pm=1\pm\iota\) and \(q_\pm i_\pm=[2]\).

The map \(\tau=q_-\rho i_+\) is an isogeny of degree three. Indeed \(q_+\rho i_+=q_+\rho^{-1}i_+=[-1]\), by conjugating with \(\iota\) and adding the two expressions. Hence
\[
 \tau^\dagger\tau
 =q_+\rho^{-1}(1-\iota)\rho i_+
 =[2]-q_+\rho i_+=[3].
\]
For non-CM factors, Hom is a rank-one lattice, so an isogeny of prime degree generates it. The degree-three relation for the dihedral family is classical. [Shaska and Völklein, *Elliptic subfields and automorphisms of genus 2 function fields*](https://citeseerx.ist.psu.edu/document?doi=e3fa3ea10fd50e2c61ddd6570700daf5dac38fe0&repid=rep1&type=pdf)

On \(E_+[2]\) one has
\[
 i_-\tau=(1-\iota)\rho i_+=(\rho-\rho^{-1})i_+=i_+.
\]
The last equality uses \(1+\rho+\rho^{-1}=0\) and the fact that minus and plus coincide on two-torsion. Thus the isogeny induces exactly the labelled two-torsion gluing. In the degree law for non-CM bielliptic Jacobians, \(k=3\) and \(|E_+[2]^\pi|=4\). The argument gives the gluing type directly from the curve action, without an isogeny-class enumeration.

## From the exact degree law to the multiplicity formula

Use the [degree law](README.md#degree-law-for-every-cyclic-minimum-isogeny-degree) for primitive directions \((a,b)\) modulo simultaneous sign. Their degree is
\[
 d(a,b)=\frac{2(a^2+3b^2)}{s^2\,3^\epsilon},\qquad
 s=\begin{cases}2&a,b\text{ both odd},\\1&\text{otherwise},\end{cases}
 \qquad \epsilon=\mathbf 1_{3\mid a}.
\]
The assignment \((a,b)\mapsto(a+b\sqrt{-3})/s\) is a bijection with nonzero primitive Eisenstein integers modulo sign. Here primitive means indivisible by an integer greater than one. To invert it, write an Eisenstein integer with numerator coordinates of the same parity and denominator two. If the coordinates are even, cancel two. In that case the remaining coordinates cannot both be odd, since otherwise the integer would be divisible by two in the Eisenstein ring. Odd common factors are excluded by primitivity in either case.

The degree is twice the Eisenstein norm with its ramified factor at 3 removed. The ramified exponent is zero or one by primitivity and equals \(\epsilon\). A primitive Eisenstein integer has no inert prime in its norm and uses only one prime above each split rational prime. Unique factorization then gives two choices for each prime dividing \(m\), two choices for the ramified exponent, and three units modulo sign. This gives \(6\cdot2^{\omega(m)}\) classes. The [finite-kernel counting theorem](../isotypic-gluing/README.md) gives the cumulative leading coefficient \(3\sqrt3/(2\pi)\).

The two elliptic factors are isogenous to one another. If one is non-CM, any elliptic quotient of the Jacobian is non-CM. Equivariance would force a nontrivial order-three elliptic automorphism fixing the image of the fully ramified source point, or factor the map through the rational source quotient. Neither is possible. This explains why the equivariant solver misses the entire elliptic spectrum on these inputs.

## Exactly which optimal quotients descend to \(\mathbb Q\)

In the coefficient chart \(Y^2=P X^6+Q X^3+R\), the rational optimal spectrum consists of **two degree-two classes if \(R/P\) is a rational sixth power**, and is empty otherwise. The criterion can rule out every nonconstant rational elliptic map despite the six geometric degree-two classes.

Choose \(\delta^3=R/P\), \(s^2=R/P\) and write the reciprocal lift as \(\iota(X,Y)=(\delta/X,sY/X^3)\). Galois sends the rotation \(\rho\) to \(\rho^{\pm1}\). It sends \(\iota\) to \(\pm\rho^j\iota\), with \(j\) determined by the conjugate of \(\delta\) and the sign by that of \(s\). Over the algebraic closure the endomorphism algebra is \(M_2(\mathbb Q)\), since the elliptic factors are isogenous and non-CM. The elements \(1,\rho,\iota,\rho\iota\) span it. Indeed, \(1,\rho\) generate a quadratic field and conjugation by \(\iota\) acts on it nontrivially. A rational elliptic direction has a Galois-invariant symmetric projector, so its rational line in the two-dimensional matrix representation must be invariant under every matrix implementing the Galois conjugation.

If \(R/P\) is not a rational cube, the cubic for \(\delta\) remains irreducible over \(\mathbb Q(\zeta_3,s)\), whose degree is a power of two. Some Galois element fixes \(\rho,s\) and sends \(\iota\) to \(\rho\iota\). Its action on the matrix algebra is conjugation by a nontrivial power of \(\rho\), which has irreducible minimal polynomial \(T^2+T+1\) and hence no invariant rational line.

If \(\delta\) is rational but \(s\notin\mathbb Q(\zeta_3)\), some element fixes \(\rho\) and negates \(\iota\). It acts by conjugation by \(w=\rho-\rho^{-1}\), since \(w\) anticommutes with \(\iota\). Here \(w^2=-3\), so there is no invariant rational line. If \(s\in\mathbb Q(\zeta_3)\setminus\mathbb Q\), conjugation in that field both inverts \(\rho\) and negates \(\iota\). It acts by conjugation by \(w\iota\), with \((w\iota)^2=3\), which again has no invariant rational line.

Finally, if \(\delta,s\in\mathbb Q\), the ratio is the sixth power of \(s/\delta\). Galois fixes \(\iota\) and acts on the matrix algebra through conjugation by \(\iota\) when it inverts \(\rho\). The only invariant directions are the two eigenspaces of \(\iota\). The reciprocal lifts are rational, and their degree-two quotients have a rational point from the source point at infinity. Both directions therefore give rational optimal maps. For the excluded cases, a nonconstant rational map would still induce a rational elliptic direction by taking the image of its adjoint, so the same obstruction applies in every degree.

## Replay

~~~sh
python -B research/elliptic-quotients/bielliptic-reconstruction/verify_cubic.py
python -B research/elliptic-quotients/bielliptic-reconstruction/cubic_genus_two.py source.json --bound 150
~~~

Input JSON uses exponent 3 and ascending coefficients. The replay compares the multiplicity formula with primitive-direction enumeration at every degree through 10,000. It checks the exact change of model for both characters, a case with nonrational quotient \(j\)-invariants, and a CM input on which the non-CM spectrum abstains. The replay includes examples with and without rational quotients. With Python-FLINT installed, add **--finite-fields** to independently exclude a rational elliptic factor on three fixtures by irreducible Weil polynomials at good primes. A uniform CM spectrum remains open in this implementation.
