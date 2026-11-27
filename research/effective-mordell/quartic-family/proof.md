# Every rational point on the quartic family through Sutherland row 36518

For every integer \(k\), the Crest–Root family through row 36518 has exactly two or four rational points. Its best uniform upper bound is four. The proof below starts from the equation and gives a terminating algorithm for the complete set of points.

## Equation and theorem

Use projective coordinates \([u:v:t]\). Write
\[
A=4u^2-8uv+2v^2,\qquad
B=-2u^4+2u^3v-4u^2v^2+8uv^3-3v^4,\qquad r=t^2+5kuv.
\]
Let \(C_k\) be the plane quartic \(r^2+Ar+B=0\). The archived row 36518 corresponds to \(k=0\), with catalogue coordinates \([x:y:z]=[u-2v:u:v+t]\). The original quartic and the family equation are in the Crest–Root family-sextics archive under the row-36518 directory. The underlying catalogue is Andrew Sutherland's database of smooth plane quartics. [Catalogue](https://math.mit.edu/~drew/genus3curves.html)

**Theorem.** The rational points on \(C_k\) are \([0:1:1]\) and \([0:1:-1]\), together with \([1:1:m]\) and \([1:1:-m]\) exactly when \(1-5k=m^2\) for a positive integer \(m\). Both outcomes occur infinitely often. Every positive \(k\) gives two points. Taking \(m=5j+1\) and \(k=(1-m^2)/5\) for \(j\ge 1\) gives four.

## Proof

All members reduce to the same quartic modulo 5. Over the algebraic closure of the residue field its three partial derivatives have no common projective zero. The accompanying verifier computes the unit Gröbner basis in each of the three affine patches. The projective singular locus is closed over the base by properness, so smooth reduction implies smoothness in characteristic zero. A smooth plane quartic has genus three.

Put \(w=2r+A\). On \(C_k\) one has
\[
w^2=(2u)^3(u-v)+16(u-v)^4.
\]
Thus whenever \(u\ne v\), a rational source point gives a rational point on \(E:y^2=x^3+16\) by \(x=2u/(u-v)\), \(y=w/(u-v)^2\). The group \(E(\mathbb Q)\) consists of \(O,(0,4),(0,-4)\). An independent proof of that arithmetic input follows below. [LMFDB 27.a4](https://www.lmfdb.org/EllipticCurve/Q/27/a/4)

For \(u\ne v\) this forces \(u=0\), hence \(v\ne0\). Set \(v=1\). The equation becomes \((t^2-1)(t^2+3)=0\), with precisely the two rational roots \(t=\pm1\). For \(u=v\), the point \([0:0:1]\) fails the quartic equation, so scale to \(u=v=1\). The quartic becomes \((t^2+5k-1)^2=0\). An integer is a rational square exactly when it is an integer square, and \(1-5k\) is nonzero for integral \(k\). This proves the classification, including the fibre at infinity. At \(k=0\) the four points in the original catalogue coordinates are \([1:0:0],[-1:0:1],[-1:1:0],[-1:1:2]\).

## Independent proof of the elliptic input

Suppose \(x,y\) are rational with \(y^2=x^3+16\) and \(x\ne0\). Define
\[
X=\frac{x^3+64}{x^2},\qquad Y=\frac{y(x^3-128)}{x^3}.
\]
Expansion gives \(Y^2=X^3-432\), so \(X\ne0\). The rational numbers \(p=(36+Y)/(6X)\) and \(q=(36-Y)/(6X)\) satisfy \(p^3+q^3=1\). Fermat's theorem for exponent three, after clearing denominators, gives \(pq=0\). Consequently \(Y=\pm36\) and \(X=12\). The original \(x\) must now be a rational root of the monic polynomial \(x^3-12x^2+64\), hence an integer. Its values at the five residue classes modulo 5 are \(4,3,4,3,1\), a contradiction. The remaining affine possibility \(x=0\) gives \(y=\pm4\), and there is one point at infinity. [Fermat's theorem for exponent three in Mathlib](https://leanprover-community.github.io/mathlib4_docs/Mathlib/NumberTheory/FLT/Three.html)

## Relation to the open problems

ProofAtlas #87 asks for a terminating method for every smooth projective curve of genus at least two over every number field. Problem #151 asks for a uniform bound for all curves of a fixed genus over a fixed number field. The theorem gives the complete point set and the sharp bound four for this family over \(\mathbb Q\). Pulling back points from a rank-zero elliptic quotient is classical. The role of Crest–Root here is the explicit quotient and the family on which it persists. No priority claim is made for the family point count. [ProofAtlas ranking](https://www.proofatlas.ai/open-problems/) · [Poonen, Computational aspects of curves of genus at least 2](https://math.mit.edu/~poonen/papers/ants2.pdf)

## Replay

With Python 3.11 and SymPy installed, run

    python verify_points.py 0
    python verify_points.py 1
    python verify_points.py -3

These give the catalogue row, a two-point member, and a four-point member, respectively. The script checks the identities in symbolic \(k\), the common smooth reduction, the arithmetic reduction to the Fermat cubic, and the listed points at the requested parameter. Completeness follows from the proof above.
