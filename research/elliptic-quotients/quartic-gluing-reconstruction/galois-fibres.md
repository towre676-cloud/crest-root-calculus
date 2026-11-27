# Rational degree-four quotients with Galois-permuted fibres

Every rational optimal degree-four elliptic map from a genus-two quintic or sextic can be recovered by exact enumeration of the Galois-stable Weierstrass pairings and a two-variable incidence. The incidence has full numerator rank at every genuine solution, even when the elliptic 2-torsion is irreducible over \(\mathbb Q\). The implementation computes both maps and their signed 4-torsion graph on an irreducible sextic. [Code](etale_fibres.py) · [Replay](verify_galois.py) · [Exact receipt](galois-receipt.json)

## The three fibres as a cubic algebra

For an optimal map the unique empty branch fibre is rational, and the three other fibres partition the six source Weierstrass points into pairs. Their partition is Galois-stable. Work first over the exact splitting field and enumerate the 15 pairings. For each pairing interpolate two functions in \(A=\mathbb Q[x]/(f)\) using the six Lagrange idempotents. Their values on each pair are the sum and product of its coordinates. Both functions are rational precisely when the partition is Galois-stable, because the two values together recover the monic quadratic of each pair. A rational linear combination \(\lambda\) of the functions separates the three pairs. It generates a rank-three étale subalgebra \(B=\mathbb Q[\lambda]\subset A\). The three quadratics are one equation \(q=x^2-s(\lambda)x+p(\lambda)\) over \(B\).

Put \(D=S^2\) for the empty branch fibre, with \(S\) a normalized binary quadratic dehomogenized to degree zero, one or two. The abscissa \(N/D\) takes the three values of some \(b\in B\) at the Weierstrass points. Translation removes the constant term of \(b\). For \(i=1,2\) let \(N_i\) be the remainder of \(D\lambda^i\) modulo \(f\) and \(c_i\) its \(x^5\) coefficient. A degree-four numerator has zero \(x^5\) coefficient. Thus \(N=c_2N_1-c_1N_2\) and \(b=c_2\lambda-c_1\lambda^2\), up to common scale, whenever \((c_1,c_2)\ne(0,0)\).

The solver takes the quadratic discriminant of \((N-bD)/q\) in \(B\). Its three rational coefficients must vanish. The solver discards coincidences among the three values of \(b\), common factors of \(N\) and \(S\), poles at a branch point, and degree loss at infinity. It extracts rational roots of the remaining equations by factored resultants, with exact saturated elimination available for common components. The incidence after these exclusions is finite because every point determines an optimal degree-four map and the refined Humbert form has finitely many vectors of norm \(16\). The normalized pole and branch-value data are unique for that quotient direction and pairing. The three distinct values of \(b\) give the target cubic by field norm, and the square-factor test verifies both coordinates exactly.

## Why rank loss cannot occur

Work over the algebraic closure, where the fibre quadratics split. Evaluation of binary quartics on the six Weierstrass points is injective. If the numerator equations lose rank, the space of allowed quartics has dimension three, so each of the three individual fibre values can be prescribed independently. A quartic supported on the value for pair \(i\) vanishes on the other four points, so it is proportional to \(q_jq_k\). Consequently rank loss requires
\[
  q_jq_k \bmod q_i \quad\text{to be proportional to}\quad S^2\bmod q_i
  \qquad(\{i,j,k\}=\{1,2,3\}).
\]
This is a projective condition. There are two cases for the two zeros of \(S\).

**A double zero.** Move it to infinity, so the abscissa is a quartic polynomial. Scale and translate to give it leading coefficient one and zero cubic coefficient. The three occupied fibres have residual double roots \(a,b,c\). They are distinct, since their branch values are distinct. They are the three finite critical roots, so \(a+b+c=0\). If the quadratic coefficient of the abscissa is \(A\), then \(q_a=x^2+2ax+(A+3a^2)\). The linear coefficient of \(q_bq_c\bmod q_a\) is \(-2a(a-b)(a-c)\). The rank-loss condition would force \(a=0\). Applying it to the other two fibres would force \(b=c=0\) as well.

**Two simple zeros.** Move the zeros to \(0,\infty\), so the abscissa is \(N/x^2\) with \(N\) monic quartic and nonzero constant coefficient. The three occupied fibres still have distinct finite nonzero residual double roots \(a,b,c\). Let \(d\ne0\) be the fourth root, with multiplicity, of \(xN'-2N\). This quartic has no \(x^2\) term, so \(ab+ac+ad+bc+bd+cd=0\). The pair quadratic is
\[
 q_a=x^2-2(b+c+d)x-abcd/a^2.
\]
In the basis \(1,1/x\) modulo \(q_a\), the coefficient of \(1/x\) in \(q_bq_c/x^2\) is
\[
 -2d(a-b)(a-c)(a+d)/a.
\]
The rank-loss condition would force \(d=-a\), and the other two conditions would force \(d=-b=-c\), again impossible. Both coefficient identities are independently checked symbolically in the replay. The proof also covers residual roots that coincide with one of the two marked Weierstrass points in their own fibre.

## An irreducible test curve

For the input
\[
 y^2=x^6-6x^5-11x^4+14x^3+12x^2-12x+1
\]
the solver finds one Galois-stable pairing. Its label algebra is the field defined by \(u^3-6u^2-16u-8\), so the three marked pairs are cyclically permuted. The two covers have half-poles \(x^2-2\) and \(x^2-1\) and target cubics \(X^3-3024X-12096\) and \(X^3-9072X+326592\). Both cubics are irreducible. At the good prime 5 the two ordinary elliptic reductions have traces \(-3,1\) and distinct geometric endomorphism fields \(\mathbb Q(\sqrt{-11})\) and \(\mathbb Q(\sqrt{-19})\). This excludes a geometric isogeny and proves that the two optimal maps exhaust the geometric spectrum of this example.

On nonzero elliptic 2-torsion the graph sends the first abscissa \(t\) to \(-t^2/36+5t/3+56\). Identical labels mean the pullbacks agree in the Jacobian. For exact order four, the signed Mumford comparison works over the abscissa field, which here has degree six. The graph gives every point in the kernel, including the ordinates, and the replay checks that its order-four part doubles to the stated 2-torsion map. The squarefree source comes from a Richelot pairing of the even sextic defined by the roots of \(z^3+z^2-2z-1\). Only the resulting source equation enters the inverse algorithm. For a separate exclusion, the replay checks that the sextic \(x^6-x-1\) with Galois group \(S_6\) has no stable pairing.

The reconstruction concerns rational maps. Maps defined only over number fields and algorithms avoiding explicit splitting fields are separate problems. The known forward classification of split Jacobians supplies context for the inverse computation. [Bruin and Doerksen](https://doi.org/10.4153/CJM-2011-039-3)

    python -B etale_fibres.py '[1,-12,12,14,-11,-6,1]' --graph
    python -B verify_galois.py
