# An index-four obstruction and complete rational degree spectra

Let \(T/\mathbb Q\) be a smooth projective genus-one curve. An elliptic target has a rational origin. A genus-one target may have none. For infinitely many geometrically different smooth plane quartics the complete degree sets over \(\mathbb Q\), with variable target, are

| Target | Degrees for positive integers \(n\) | Smallest degree |
|---|---|---:|
| Genus-one curve | \(2n^2,\ 4n^2\) | 2 |
| Elliptic curve with a rational origin | \(4n^2,\ 8n^2\) | 4 |

There are three optimal genus-one quotient directions, all of degree two. On a diagonal member the sixteen geometric cubic directions form a single Galois orbit of size sixteen. The [solver](descent.py) starts from quartic coefficients and can recover rational sign coordinates if the presentation has been mixed. An independent [replay](verify.py) checks the maps symbolically, recomputes traces by counting the genus-one double covers, and checks the family specialization and incomplete-evidence controls. The monograph’s §132 treats the index relation in degree three. The following argument applies to all degrees. [Canonical text](../../monograph-publication/README.md)

## A degree-two line bundle on every genus-one target

For a proper smooth curve \(D/K\), the index \(\operatorname{ind}(D)\) is the gcd of the degrees of its closed points. Let \(f:C\to T\) be a nonconstant morphism over a field of characteristic zero, with \(C\) smooth, projective and geometrically connected of genus \(g>1\), and \(T\) of genus one. The map is finite flat. The actual vector bundle \(V=f_*\mathcal O_C\) has Euler characteristic \(1-g\). Riemann–Roch on \(T\) gives \(\deg V=1-g\), so \((\det V)^{-1}\) is a line bundle defined over \(K\) of degree \(g-1\). It has a nonzero section over \(K\) and yields a rational effective divisor of that degree. In particular
\[
\operatorname{ind}(T)\mid g-1,\qquad
\operatorname{ind}(C)\mid \deg(f)\operatorname{ind}(T).
\]
The second divisibility comes from pullback of closed-point divisors. For a smooth plane quartic of index four, every rational genus-one map therefore has even degree. A rational elliptic origin gives the stronger divisibility \(4\mid\deg(f)\).

## Certifying index four by valuations

Work in integral coordinates at an odd prime \(p\). Multiply the quartic by a scalar and rescale coordinates by powers of \(p\) as necessary. Write \(v_i\) for the three valuations of its pure fourth powers and \(u_i\) for their unit parts. Suppose the \(v_i\) are distinct in \(\{0,1,2,3\}\), and every nonzero mixed coefficient has valuation strictly larger than \(\max v_i\). Let \(i,j\) be the two indices with the same parity. If \(-u_i/u_j\) is a nonsquare in \(\mathbb F_p\), the smooth plane quartic has index four.

Indeed, over any finite extension of \(\mathbb Q_p\) let \(e,f\) be its ramification index and residue degree. Normalize a putative point so all coordinates are integral and at least one is a unit. The least pure-term valuation is at most \(e\max v_i\). Every mixed term has strictly greater valuation. If \(e\) is odd, the pure terms have distinct valuations modulo four. If \(e\equiv2\bmod4\) and \(f\) is odd, the only possible tie is the indicated pair. Upon division by the common power of the uniformizer, cancellation would force a nonsquare to be a square in \(\mathbb F_{p^f}\). The unit from writing \(p\) as a uniformizer power times a unit contributes a square since \(v_i-v_j\) is even. These cases exhaust \(4\nmid ef\). Thus every local point field has degree divisible by four. Completing a global point field at \(p\) and summing local degrees gives the same divisibility globally. A rational line cuts out a divisor of degree four, giving equality. The test applies to arbitrary sufficiently \(p\)-adically small mixed perturbations of the three pure terms. A failed test leaves the index undecided.

## Sign quotients, maps and their arithmetic

Suppose the quartic is even in all three coordinates. Write it with \(s=X_i,\ t=X_j,\ z=X_k\) as
\[
a s^4+s^2(mt^2+nz^2)+b t^4+o t^2z^2+c z^4=0 .
\]
The sign change in \(s\) has quotient \(T_i\) in weighted coordinates \([t:z:w]\), with
\[
w=2as^2+mt^2+nz^2,\qquad
w^2=A t^4+B t^2z^2+C z^4,\quad
(A,B,C)=(m^2-4ab,\ 2mn-4ao,\ n^2-4ac).
\]
The coefficients \(a,b,c\) are nonzero on a smooth even quartic. Each sign quotient is a double cover ramified at four fixed points. Its induced map on Jacobians has connected kernel, since a factorization through a nontrivial elliptic isogeny would leave a degree-one map from genus three to genus one. If the source has index four then \(T_i\) has index exactly two: it has a degree-two divisor, and a rational point would pull back to degree two on the source.

On \(z=1\), the map \((t,w)\mapsto(x,y)=(At^2,Atw)\) has degree two onto
\[
 E'_i:\quad y^2=x^3+B x^2+ACx .
\]
It extends over the points at infinity by properness. Its deck involution sends \((t,w)\) to \((-t,-w)\). The usual 2-isogeny with kernel generated by \((0,0)\) has target
\[
 E_i:\quad Y^2=X^3-2B X^2+(B^2-4AC)X .
\]
Indeed \(X=x+B+AC/x,\ Y=y(1-AC/x^2)\). The composite from \(T_i\) has degree four and gives \([2]\) on Jacobians after identifying \(\operatorname{Jac}(T_i)\simeq E_i\). For this identification, the cover \(T_i\to E'_i\) adjoins \(\sqrt{x/A}\). Over an algebraic closure it is the dual of the displayed 2-isogeny, seen by the divisor \(2[(0,0)]-2[O]\) of \(x\). Changing the constant \(A\) twists the covering by translation by its order-two deck point, which acts trivially on the Jacobian. Equivalently these are the binary-quartic Jacobian formulas. [Cremona, Fisher and Stoll, genus-one models](https://www.dpmms.cam.ac.uk/~taf1000/papers/minred-234.pdf)

The pullback differential lines of the three quotients are independent. On a canonical plane quartic sign change acts on differential numerators by the coordinate change times its determinant. Its invariant line is generated by \(X_i\). The three pullbacks therefore give a rational isogeny from the product of the three elliptic Jacobians to \(J(C)\). At one good prime, unequal traces exclude a rational isogeny between two factors. If all three pairs are separated, every elliptic quotient of \(J(C)\) over \(\mathbb Q\) has one of these three connected kernels. Via the Albanese torsor of \(C\), every rational map to a genus-one curve factors through its corresponding \(T_i\). More explicitly the factorization exists over the algebraic closure by connectedness and the Jacobian universal property. Its intermediate map descends because the original map and the surjective map \(C\to T_i\) are rational.

For each \(E_i\), at a good prime \(p\ne\ell\), a rational cyclic \(\ell\)-isogeny requires \(x^2-a_p x+p\) to have a root modulo \(\ell\). We eliminate each odd prime on Mazur’s rational-isogeny list by exhibiting a nonsquare discriminant. If \(16AC\) and \(B^2-4AC\) are both nonsquares over \(\mathbb Q\), \(E_i\) and \(E'_i\) each have only one rational subgroup of order two. The only rational order-two subgroup after the first isogeny is its dual kernel. A rational cyclic order-four subgroup is therefore impossible. Higher cyclic orders would have a characteristic subgroup of a forbidden order. For a finite rational isogeny kernel of type \((n,nr)\), divide out \(E_i[n]\). What remains is a rational cyclic subgroup of order \(r\), so \(r=1\) or \(2\). The possible isogenies, up to target isomorphism, are multiplication by \(n\), and multiplication by \(n\) followed by the displayed direction to \(E'_i\). [Mazur, rational isogenies of prime degree](https://www.math.columbia.edu/~goldfeld/Mazur-Goldfeld1978.pdf)

Push these isogenies out along the genus-one torsor. The class of \(T_i\) in \(H^1(\mathbb Q,E_i)\) is nonzero of order two: it has an actual degree-two divisor and has no rational point. Multiplication by \(n\) kills that class exactly for even \(n\). The quotient isogeny to \(E'_i\) already kills it, as witnessed by the explicit degree-two map to \(E'_i\). This proves both degree lists in the table. The attainable degrees to arbitrary genus-one targets include 2, 4, 8, 16 and 18. The degree 18 maps require an unpointed target. Degree 12 is absent even though 4 divides 12.

## An infinite family with varying geometric source

Take
\[
 C_t:\quad X^4+3Y^4+63Z^4+27t X^2Y^2=0,\qquad t\in\mathbb Z .
\]
Every member is smooth since \((27t)^2\ne12\). At 3 the pure valuations are \(0,1,2\) and the competing unit ratio is a nonsquare. All members therefore have source index four, minimum rational genus-one degree two, and minimum rational elliptic degree four. At \(t=1\) the traces of the three quotient Jacobians at 5 are \(2,-4,0\). The full comparison uses the following prime data in axis order \(X,Y,Z\).

| Prime | 5 | 11 | 17 | 19 | 23 | 29 |
|---|---:|---:|---:|---:|---:|---:|
| \(a_p(E_X)\) | 2 | 0 | −2 | 0 | 0 | −4 |
| \(a_p(E_Y)\) | −4 | 0 | 8 | 0 | 0 | 10 |
| \(a_p(E_Z)\) | 0 | −4 | 0 | −2 | 0 | 4 |

Both degree lists hold for every \(t\equiv1\bmod M\) with \(M=5\cdot7\cdot11\cdot17\cdot19\cdot23\cdot29\). The listed good prime reductions and all odd-order exclusions persist by congruence. Write \(m=27t\), \(\delta=m^2-12>0\) on this progression. On the first two axes the binary quartics have \(B=0\), \(A=\delta\) and \(C=-252,-756\), respectively. Thus \(AC<0\) and \(-4AC\) has odd valuation at 7. On the third axis \(AC\) has square class 3 and \(B^2-4AC\) has square class \(\delta\), strictly between \((|m|-1)^2\) and \(m^2\). This checks the rational 2-subgroups uniformly. The degree-two target on the third axis has
\[
j=16\,\frac{(m^2+36)^3}{3(m^2-12)^2}.
\]
This nonconstant function takes infinitely many values on the progression. A fixed genus-three source has finite automorphism group, hence finitely many involution quotients. The progression consequently contains infinitely many geometric isomorphism classes of source curves.

## Sixteen conjugate cubic directions on a diagonal twist

On \(C_0:X^4+3Y^4+63Z^4=0\) the same local and rational-isogeny tests separately give the two complete degree lists. Geometrically it becomes Fermat on rescaling \(Y\) by \(\sqrt[4]3\) and \(Z\) by \(\sqrt[4]{63}\). The monograph’s Theorem 128.1 determines the sixteen geometric cubic differential lines of Fermat. On the present quartic they are
\[
[1:i^j\sqrt[4]3:i^k\sqrt[4]{63}],\qquad 0\le j,k<4.
\]
Over \(\mathbb Q(i)\) the valuation rows for \(3,63\) at primes above 3 and 7 are \((1,2),(0,1)\). These primes are unramified in \(\mathbb Q(i)\). The determinant is odd, so the classes generate a copy of \((\mathbb Z/4)^2\) modulo fourth powers. Kummer theory gives all sixteen sign-and-phase choices transitively. Each cubic direction therefore has field of definition of degree exactly sixteen over \(\mathbb Q\). One is defined over \(\mathbb Q(\sqrt[4]3,\sqrt[4]{63})\) by transport from Fermat. The other directions are its conjugates. The [original classification](../../monograph-publication/README.md) supplies the geometric count; the local index obstruction and orbit calculation explain the arithmetic non-descent on this twist.

## Replay and input scope

From this folder run

~~~sh
python descent.py diagonal.json
python descent.py deformation.json
python descent.py mixed.json --recover-frame
python verify.py
~~~

Inputs use the same rational quartic-monomial format as the atlas. The verifier checks smoothness over the algebraic closure. The local obstruction also applies directly to smooth inputs with non-even monomials in the displayed coordinates. For the stronger isogeny classification, the optional frame recovery factors the Hessian of a diagonal twist or recovers rational sign centers using the [involution atlas](../involution-atlas/README.md). Each accepted substitution is checked by exact polynomial arithmetic. A limited prime search or an unavailable sign frame can leave the degree sets undecided. The output includes only conclusions established for that input.

