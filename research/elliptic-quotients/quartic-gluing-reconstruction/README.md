# Degree-four elliptic maps and the signed gluing from one genus-two equation

The degree-four Richelot control \(C'\) now has an exact reconstruction from its equation. The solver finds both elliptic curves, both maps and the full graph identifying their four-torsion. It receives no target coefficients or Richelot correspondence. An independent replay verifies the local degree at 37 and 61 and compares one recovered map with the archived Richelot map.

The same reconstruction works on all **36 degree-four neighbours** in the five-source sweep. It gives two rational optimal maps and the full signed 16-point graph on every equation. On 32 curves the factors are geometrically distinct, as verified at a good prime, so the two recovered maps give the complete geometric optimal spectrum. On the four repeated-factor curves the two recovered rational maps give a splitting. Their additional geometric degree-four classes are described in the separate finite-kernel spectrum. [Earlier sweep](../../torsion-denominators/README.md#family-sweep) · [Repeated-factor spectrum](../exact-completeness/README.md)

## The control \(C'\)

Write
\[
f(x)=11(x^2-10x+13)(2x^2+2x-7)(7x^2+18x+3),\qquad C':y^2=f(x).
\]
The two maps have coordinates \((N(x)/S(x)^2,H(x)y/S(x)^3)\).

| \(S\) | Target \(Y^2=g(X)\) | \(11N\) | \(121H\) |
|---|---|---|---|
| \(x^2-3\) | \(X^3-24003X+1296702\) | \(-324(14x^2+30x-79)-1866S^2\) | \(486(x+5)(2x-1)(7x-9)\) |
| \(x^2-5\) | \(X^3-21168X-494208\) | \(-144(30x^3+139x^2-393)-552S^2\) | \(432(x+7)(3x-1)(5x+13)\) |

Both identities \(S^6g(N/S^2)=fH^2\) hold in \(\mathbb Q[x]\). The graph below maps the first elliptic factor to the second, using the signs of these maps. The identities map to one another. The three nonzero two-torsion abscissas are paired as \((-177,-132),(66,-24),(111,156)\); their ordinates are zero.

| Abscissa \(t\) on the first factor satisfies | Abscissa on the second factor | Multiplier on the ordinate |
|---|---|---|
| \(t^2-222t-639=0\) | \(2t-66\) | \(2t/243+58/81\) |
| \(t^2-132t+15291=0\) | \(-4t/3+64\) | \(t/81-31/27\) |
| \(t^2+354t-38655=0\) | \(2t/3-14\) | \(-2t/405-74/27\) |

The replay verifies both maps and the graph again on the reciprocal chart \(x=2+1/u\), checks both targets and matches the abscissa functions by target projectivities in the original coordinates. The second map agrees with the archived Richelot direction by exact projective and differential comparisons. [Archived map](../../torsion-denominators/map4.json) · [Exact receipt](receipt.json)

## Elimination in the paired branch chart

Suppose each of three fibres over elliptic two-torsion images contains two source Weierstrass points, and the fourth contains none. Choose the empty branch value as the target origin. The map on the quotient lines has degree four. Its empty fibre has even multiplicity everywhere, so its denominator is a square \(D=S^2\) of a binary quadratic. On the three other fibres the numerator satisfies
\[
N-r_i D = q_i L_i^2 c_i
\]
where \(q_i\) is a quadratic on the two marked source Weierstrass points, \(L_i\) is a binary linear form, and \(c_i\ne0\). This follows from the ramification parity in the two double covers. When the elliptic two-torsion is rational, each pair is stable over the rationals. The branch factorizations in this study give a unique rational pairing.

Normalize \(S\) and use translation on the abscissa to delete the coefficient of \(x^{2\deg S}\) from \(N\). Its four remaining coefficients satisfy one linear equation per \(q_i\): the remainder of \(N\) must be proportional to the remainder of \(D\) modulo \(q_i\). The signed \(3\times3\) minors give the projective coefficients of \(N\) wherever the rows have rank three. For each \(q_i\), remove it from \(N-r_iD\) and require discriminant zero in the residual quadratic. Clearing the abscissa denominators gives three polynomial constraints on the at most two coefficients of \(S\). The solver intersects rational factors by resultants, then uses the square-factor test and the full equation identity to check each candidate. The charts include degree one and degree zero for \(S\), since the marked fibre may contain source infinity. A quintic input is converted to an equivalent sextic chart.

This implementation covers the rank-three paired branch charts with isolated solutions obtained from the nonzero resultants. Nonrational labels, additional rational branch pairings and degenerate numerator or resultant charts require further work. The degree-four splitting classification and Richelot factorization are due to Bruin and Doerksen. Here the computational contribution is the exact inversion of the paired-fibre equations from the source alone and the signed four-torsion readout. [Bruin and Doerksen](https://doi.org/10.4153/CJM-2011-039-3)

## Optimality and the entire kernel

Let \(H_0\) be a hyperelliptic fibre. The pullback of the elliptic origin is equivalent to \(2H_0\), because all its poles are over unmarked branch values on the source quotient line and are paired by the hyperelliptic involution. Subtracting it from a fibre over any nonzero two-torsion point leaves the two distinct Weierstrass points on its \(q_i\), minus \(H_0\). These three classes are nonzero and distinct in the Jacobian. Thus pullback is injective on elliptic two-torsion. Since its kernel is killed by four, it is injective on the entire elliptic curve. This proves the maps are optimal within this chart, including the repeated-factor controls.

Conversely, for any separable degree-four elliptic map in characteristic zero, the number of source Weierstrass points in each of the four target branch fibres is even by ramification parity. There are six in total. The only distributions are \(2+2+2+0\) and \(4+2+0+0\). In the latter distribution both empty fibres are equivalent to \(2H_0\). Their target points differ by nonzero two-torsion, which is therefore killed by pullback. Thus **a degree-four elliptic map is optimal exactly when its Weierstrass distribution is \(2+2+2+0\)**. The same argument works over the algebraic closure, so the remaining general reconstruction problem is to cover Galois-permuted pairs and degenerate elimination charts.

For order-four points choose a source Weierstrass point as infinity over a temporary field. The elliptic origin then pulls back to \(4\infty\) in the divisor class. If \(P=(t,w)\) has order four, its pullback class is represented by the quartic fibre \(N-tS^2\) and the ordinate polynomial \(w S^3/H\) modulo the fibre polynomial. The algorithm reduces these Mumford pairs exactly using \(w^2=g(t)\). Matching opposite reduced pairs for both covers supplies the signed map on the points of exact order four. The verifier also doubles each block in the targets and matches its two-torsion label. Temporary extensions disappear by rational linear algebra.

Let \(i_1,i_2\) be the two pullback embeddings. The different empty-fibre divisors distinguish the two covers and hence the two elliptic subvarieties. The map \(\chi=i_1+i_2\) is an isogeny, of degree \(16-\deg(i_1^\vee i_2)\), by the adjoint polarization matrix. The sixteen exhibited kernel points force degree sixteen and a zero off-diagonal term. Therefore the graph is the full kernel and the induced torsion isomorphism reverses the Weil pairing. Geometric non-isogeny makes the spectrum complete for the 32 distinct-factor cases. A squarefree good-prime reduction with distinct geometric elliptic Frobenius blocks certifies that condition separately on each model.

## Replay

From this directory, with SymPy, cypari2 and python-flint installed:

    python -B reconstruct_four.py '[-3003,-14850,11968,6600,-3575,-990,154]'
    python -B verify_four.py
    python -B verify_four.py --sweep

The last command regenerates the checks for all 36 frozen equations and compares them with the stored receipt. The optional reciprocal chart changes coordinates by \(x=a+1/u\), \(v=u^3y\), with \(a\) stored in the receipt. Maps are written in that chart when it is used.
