# Recovering the 3-torsion gluing from Howe's genus-two curves

Starting with Howe's genus-two equation \(C_{70}\), exact elimination recovers both rational degree-three elliptic maps. Pulling back the nonzero three-torsion points recovers the full Frey–Kani graph over \(\mathbb Q\), including its signs. The search uses the source equation alone. The same algorithm recovers both covers and the graph for all five of Howe's order-28 examples and one of his order-27 examples. [Howe's Theorems 2.1 and 2.4 and Table 2](https://doi.org/10.1112/blms/bdu107)

## The order-70 curve and its twists

Complete the square in Howe's equation to obtain
\[
v^2=(x+4)(4x^2-36x+89)(x^3+2x^2-15x+36).
\]
Put \(D=x^3+2x^2-15x+36\) and \(R=(x+4)(4x^2-36x+89)\). The two recovered maps are \(X=N(x)/S(x)\), \(Y=H(x)v/S(x)^2\) with the following coefficients.

| Target | \(S\) | \(N\) | \(H\) |
|---|---|---|---|
| \(Y^2=X^3+149013X+25726950\) | \(D\) | \(-3(23x^3+694x^2-2505x-16452)\) | \(1944(x-13)(x^2+8x+24)\) |
| \(Y^2=X^3-7483623723X+249446508217254\) | \(R\) | \(-3(20924x^3-104620x^2+1829311x-6605828)\) | \(76212576(2x^3-17x^2+40x-34)\) |

They give the factors \(E_1=66\mathrm{c}2\) and \(E_2=858\mathrm{k}1\) in this order. At the good prime 17 the traces are \(-2\) and \(4\). The reductions are ordinary, with distinct endomorphism fields \(\mathbb Q(i)\) and \(\mathbb Q(\sqrt{-13})\). An isogeny over a finite extension would preserve the geometric endomorphism algebra. Specialization therefore rules out an isogeny over \(\overline{\mathbb Q}\). Each geometric elliptic subvariety of the split Jacobian is one of its two factors. The entire optimal quotient spectrum is thus **two classes of degree three**. The independent local degree certificates give degree three at 17 and 23 from the sextic.

Let \(\psi_3(t)=t^4+298026t^2+102907800t-7401624723\). The kernel of \(E_1\times E_2\to J(C_{70})\) is the graph of
\[
(t,w)\longmapsto
\left(
\frac{7t^3-2427t^2+2768805t+178611399}{7776},
\frac{t^3-265t^2+384267t-10028907}{31104}\,w
\right),\qquad \psi_3(t)=0.
\]
The identity maps to the identity. The signs use the displayed choices of \(H\).

Howe's rational divisor has pushforwards of exact orders 10 and 7 under these maps. The dual gluing isogeny has kernel killed by three, so the divisor is killed by 210. Exact good-reduction Jacobian orders \(280,280,630\) at \(17,19,23\) bound rational torsion by 70. The pushforward orders force 70, giving \(\mathbb Z/70\mathbb Z\).

For every squarefree \(d\equiv5\bmod391\), the twist \(v^2=dDR\) has the same geometric spectrum and has two rational degree-three maps to the corresponding twisted targets. For a short target \(Y^2=X^3+aX+b\), evaluate the same map formulas on the twisted source, multiply both outputs by \(d\), and use target coefficients \(d^2a,d^3b\). The reduction orders are now \(352\) at 17 and \(486\) at 23. They bound the torsion by two, and the rational two-torsion of \(66\mathrm{c}2\) survives twisting and the degree-nine quotient. Thus these twists have group \(\mathbb Z/2\mathbb Z\). This gives an infinite family with the same geometric degree data as \(C_{70}\) and different rational torsion.

## Exact reconstruction

For a rational degree-three map, the hyperelliptic involution of the source induces a translate of negation on the target. On the quotient lines, the map has degree three. Each of the four fibres over target branch points contains an odd number of source Weierstrass points, by the ramification parity in the double covers. All six source Weierstrass points lie over those fibres, so their partition is \(3+1+1+1\). The unique triple is Galois stable. Move its image to the elliptic origin and write the map on the quotient lines as \(P/D\), with \(D\) a cubic factor of the sextic. An affine change in the target \(X\) makes \(P\) a projective polynomial of degree at most two.

At each root \(a\) of the complementary cubic, the residual quadratic
\[
\frac{D(a)P(x)-P(a)D(x)}{x-a}
\]
has discriminant zero. Reducing this discriminant modulo the complementary cubic gives three homogeneous quadrics on the coefficients of \(P\). The solver enumerates every rational cubic factor and all rational solutions of the projective charts. A resultant gives the three other target branch values. A square-factor test constructs the \(Y\) coordinate, and direct substitution certifies the map and its degree.

The projective charts are finite. In fact, if \(P,D\) shared a root \(u\), every residual quadratic would have a double root at \(u\). Since \(D\) has simple roots, this would give \(P(a)/D(a)=P'(u)/D'(u)\) at all three complementary roots \(a\). The cubic \(P-P'(u)D/D'(u)\) would vanish at those three roots and at \(u\), contradicting \(P\ne0\) and \(\deg P\le2\). With \(P,D\) coprime, two complementary roots cannot map to the same value: each residual quadratic would have to be the double of the other point, which would give incompatible multiplicities in the same cubic fibre. The resulting three distinct values give the target, and the even residual multiplicities give the square factor over \(\mathbb C\). Thus every projective solution gives a geometric degree-three cover. Such covers are finite modulo translation and isomorphism of target, since their classes are primitive representations of 9 by the positive-definite refined Humbert form. [Kani's correspondence](https://doi.org/10.1007/s40316-018-0105-6)

For the graph, use each irreducible factor of the first division polynomial as the field of an \(X\)-coordinate. Factor the second division polynomial there to enumerate possible partners. Compute the reduced Mumford representation of each pulled-back three-torsion point. Moving a source Weierstrass point to infinity allows ordinary Cantor reduction. The linear part of the pair has the form \(wV\), where \(w^2\) is the elliptic cubic at the first coordinate. Thus the arithmetic only needs \(w^2\). Equality of the reduced pairs with opposite signs gives the signed multiplier \(B\) in \((X,Y)\mapsto(A(X),B(X)Y)\). If a temporary extension is needed to choose a source Weierstrass point, exact linear algebra expresses \(B\) back in the first coordinate field.

The two prime-degree maps induce elliptic embeddings \(i_1,i_2\) and the recovered pairs give nine points of \(\ker(i_1+i_2)\). Distinct targets give distinct elliptic subvarieties in the replay. The degree of \(i_1+i_2\) is \(9-\deg(i_1^\vee i_2)\) by the adjoint matrix, with \(\deg(0)=0\). Hence the kernel has exactly nine points and the off-diagonal homomorphism is zero. The pulled-back polarization is three times the product polarization. This proves the graph is the full kernel and its map is a Weil anti-isometry. [Bröker–Howe–Lauter–Stevenhagen, §5 and Appendix A](https://arxiv.org/abs/1403.6911)

The order-27 example tests the prime where gluing matters. On the completed-square model its divisor \([(-1,6)]-[(-1,-6)]\) pushes forward to points of orders three and nine, so the divisor is killed by 27. Reduction at 7 shows its ninefold multiple is nonzero. The Jacobian orders 27 and 108 at 5 and 7 give an exact upper bound of 27. Here the Jacobian has a rational point of greater order than every rational point of the product of the two factors. The other two controls have no rational cubic branch factor and hence no rational map of degree three.

## Forward search and the orientation test

Enumerate Kubert–Tate order-7 and order-10 normal forms at reduced rational parameters of height at most 100. Every nonsingular specialization is checked by computing the order of the marked point. Modulo \(\mathbb Q\)-isomorphism, the families give 6,087 and 9,130 elliptic curves. The **55,574,310** possible pairs leave **183,020** after the first twelve prime tests and **two** after all primes from 5 through 251. A test uses a prime only when both curves have good reduction there. A Galois-equivariant gluing implies congruence of the two traces modulo three. The two surviving pairs use \(858\mathrm{k}1\) and the order-10 models \([1,0,0,115,561]\) and \([1,0,0,-45,81]\). [Howe's database and specialization search](https://arxiv.org/pdf/1407.2654) · [Kubert–Tate order-ten parametrization](https://alozano.clas.uconn.edu/wp-content/uploads/sites/490/2014/01/gonzalez-lozano-torsion_quartic_v10b-REV1-E2.pdf)

The two order-10 curves are rationally 2-isogenous. Their traces therefore coincide at every common good prime. The 2-isogeny multiplies the exponent in the Weil pairing on three-torsion by \(2=-1\bmod3\). Also, \(66\mathrm{c}2\)'s four projective three-torsion lines have full \(S_4\) Galois action. Its centralizer in \(S_4\) is trivial, so every equivariant automorphism of its three-torsion is scalar \(\pm I\) and preserves the pairing. Composing the recovered anti-isometry with the 2-isogeny gives an isometry to the second target, and every equivariant isomorphism to that target has this orientation. Only Howe's pair can glue anti-isometrically in the specified height box.

## The only point order above 70 in a rational split 3-gluing

Suppose \(A=E_1\times E_2\), and a rational 3-torsion graph gives dual isogenies
\(\chi:A\to J\) and \(\widehat\chi:J\to A\), with both compositions equal to multiplication by 3.
As Galois modules, \(J[3]\) has a two-step filtration with each term isomorphic to \(E_1[3]\): the submodule is the image of \(A[3]\) under \(\chi\), and the quotient is the graph under \(\widehat\chi\).

If the elliptic factors have no rational point of order 3, then \(J[3](\mathbb Q)=0\). The prime-to-three part of rational torsion is isomorphic to that of the product, so Mazur's list bounds every point order by 70. If the elliptic factors have rational 3-torsion, their group exponents are in \(\{3,6,9,12\}\). Since \(\chi\widehat\chi=[3]\), the order of a point on \(J\) divides three times the least common multiple of the two elliptic exponents. The only multiple candidate above 70 is **108**, from exponents 9 and 12. [Mazur's theorem, as stated by Darmon](https://www.claymath.org/wp-content/uploads/2022/03/cmip08c.pdf)

There is a local test for excluding 108. Let \(p\ne3\) be good for both factors. Suppose the order-9 factor has \(3\)-primary group of size nine at \(p\), and the order-12 factor has cyclic \(3\)-primary group of size at least nine. A rational point \(P\) of order 27 on \(J\) would have dual image \((V,W)\) with \(V\) of order nine and \(3W=0\). The stipulated groups make \(W\) divisible by three locally and \(V\) indivisible by three locally. To see the contradiction, choose \(R=(R_1,R_2)\) above \(P\) under \(\chi\), so \(3R=(V,W)\). Subtract a graph point from \(R\) to make \(R_2\) the rational local third of \(W\). The Frobenius difference of \(R\) then lies in the graph with zero second coordinate, so it vanishes. This makes \(R_1\) a rational local third of \(V\), which cannot exist.

At height at most 100 the marked order-9 and order-12 families give 6,087 and 7,608 curves, or **46,309,896 pairs**. An isomorphism on three-torsion must identify the finite algebra of the four projective lines. Exact number-field isomorphism tests leave four pairs, each of which fails the local lift test.

| Order-9 parameter | Order-12 parameter | Obstructing prime | Order of the first group | Invariant factors of the second group |
|---:|---:|---:|---:|---|
| \(67/6\) | \(-9\) | 71 | 72 | \(36,2\) |
| \(-1\) | \(1/3\) | 37 | 36 | \(36\) |
| \(5\) | \(1/3\) | 37 | 36 | \(36\) |
| \(88/5\) | \(1/3\) | 41 | 36 | \(36\) |

Thus no rational split 3-gluing in these order-9 and order-12 height boxes can reach 108. The unrestricted order-108 case remains open in this study. [Kubert–Tate order-nine parametrization](https://deepblue.lib.umich.edu/bitstream/handle/2027.42/100026/zscherr_1.pdf)

## Replay

Use Python with SymPy, cypari2 and python-flint. Run from the repository root.

    cd research/elliptic-quotients/torsion-gluing-reconstruction
    python -B reconstruct.py '[12816,-7320,817,690,-155,-12,4]' --graph
    python -B verify.py --local --receipt /tmp/howe-replay.json
    python -B sieve.py --height 100 --prime-bound 251
    python -B sieve.py --height 100 --prime-bound 251 --order-108

The [receipt](receipt.json) has both maps and every signed graph component for all seven split examples. It uses ascending coefficients. This study recovers known curves and gives an exact inverse implementation for rational degree-three maps on the supported sextic charts. An equation-only degree-four gluing for the earlier Richelot neighbour remains open.
