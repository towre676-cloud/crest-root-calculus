# Opposite secants and elliptic-quotient spectra

Distinct genus-two curves can have identical counts of optimal elliptic quotients in **every degree prime to 3**. An explicit infinite family of such pairs has an unbounded first distinguishing degree, even within a single CM field as the conductor varies. The same construction explains all 31 coincidences in the level-three census.

## The space observed by good degrees

Write $r_L(n)$ for the number of vectors of norm $n$ in a ternary lattice $L$, including both signs. Write $u_L(m)$ for the number of primitive vectors of norm $m^2$. On the refined Humbert lattice of a genus-two Jacobian the primitive vectors correspond to the geometric optimal elliptic subcovers. On a principally polarized product the two norm-one vectors correspond to the two elliptic factors. [Kani, The isogeny defect](https://doi.org/10.1007/s40316-018-0105-6)

Fix one full ternary genus. Let $A_p$ sum a function over the $p$-neighbours, with multiplicities. For a good prime $p\nmid 2D\ell$, put $\chi_D(p)=(D/p)$ and define
$$
 P_{p,-1}=0,\qquad P_{p,0}=1,\qquad P_{p,1}=X-\chi_D(p),\qquad
 P_{p,k+1}=XP_{p,k}-pP_{p,k-1}\quad(k\ge1).
$$
Let $R$ be the set of product classes. For every $m$ prime to $2D\ell$,
$$
 (r_L(m^2))_L=2\prod_{p^k\parallel m} P_{p,k}(A_p)\mathbf1_R,\qquad
 (u_L(m))_L=2\prod_{p^k\parallel m}
          (P_{p,k}-P_{p,k-1})(A_p)\mathbf1_R.
$$
Apply the [good-prime neighbour recurrence](../good-degree-spectrum/README.md) at $n=p^{2j}m'^2$ with $p\nmid m'$. Its middle symbol $(Dn/p)$ is $\chi_D(p)$ for $j=0$ and zero for $j>0$. Iterate in $j$ and then over the prime factors of $m$. The initial vector is $2\mathbf1_R$. Indeed, all forms in the genus take only $0,1\bmod4$. Two independent norm-one vectors would have integral inner product of absolute value less than one, hence their sum would have norm two. The second formula follows from $u_L(m)=\sum_{d\mid m}\mu(d)r_L((m/d)^2)$. These are signed polynomial weights on walks in the ternary neighbour graph.

Both families of count vectors span exactly $O=\mathbf H_{\mathrm{good}}\mathbf1_R$. Each prime-power polynomial is monic of its indicated degree, so products of powers of the Hecke operators expand in products of the prime-power polynomials. This gives the reverse inclusion. The joint orthogonal projector $E_\lambda$ for a good Hecke eigensystem $\lambda$ contributes the one-dimensional span of $E_\lambda\mathbf1_R$ when that vector is nonzero, even when the eigensystem has multiplicity.

For a basis matrix $B$ of $O$ and $W=\operatorname{diag}(1/|\operatorname{Aut}(L)|)$ the observed fraction of a pair difference is
$$
 \operatorname{Vis}(i,j)=
 \frac{(B^{\mathsf T}(e_i-e_j))^{\mathsf T}
       (B^{\mathsf T}WB)^{-1}
       (B^{\mathsf T}(e_i-e_j))}
      {W_{ii}^{-1}+W_{jj}^{-1}} .
$$
Here $W^{-1}(e_i-e_j)$ represents the difference of evaluation functionals in the mass inner product, and the orthogonal projector onto $O$ is $B(B^{\mathsf T}WB)^{-1}B^{\mathsf T}W$. A basis computed from a subset of good primes gives a rigorous lower bound. The verifier checks the two polynomial identities against direct enumeration for 1,120 ordinary and primitive entries across three full lattice levels. [Exact certificate](certificate.json)

## The level-three classification

Let $G$ be an ambient ternary lattice in the genus of $\mathbb Z\oplus O_D$ with norm $u^2+4N(c)$, and suppose 3 is inert. Its isotropic conic $X$ over $\mathbb F_3$ has four points. A split plane is a secant indexed by a two-element subset $e$ of $X$. Let $L_e$ be its integral inverse image. Isometries between such lattices extend to their unique integral index-three ambient overlattices. In a fixed $G$, the classes are therefore the secant orbits of $O(G)$. [Ambient lattices and forced incidence relations](../forced-square-kernel/README.md)

An internal point is a projective point whose orthogonal plane misses the conic. It is exactly a point of nonzero square norm in this local model. For any internal point and conic point there is a unique secant through both. The three internal points therefore give the three matchings of $X$. Distinct internal points cannot give the same matching because two secants intersect in just one point. The only internal point on $e$ is consequently the only internal point on the opposite secant $e^*=X\setminus e$.

**Theorem 1.** The ordinary and primitive norm-$m^2$ counts of $L_e$ and $L_{e^*}$ agree for every $3\nmid m$. Within each ambient block the forced split-genus kernel has a basis consisting of the differences of the distinct secant orbits exchanged by $e\mapsto e^*$.

Every ambient vector of norm $m^2$ with $3\nmid m$ reduces to an internal point. It lies in $L_e$ exactly when it lies in $L_{e^*}$. This proves ordinary equality, and divisor inversion gives primitive equality. The involution on secants commutes with the ambient automorphisms. Over each matching orbit there are either two exchanged secant orbits or a single fixed one. Matching orbits are the internal-point orbits, hence the nonsplit plane orbits by polarity. The number of exchanged pairs is $h_{\rm split}-h_{\rm nonsplit}$ in this block. Their differences have disjoint supports and this number is the dimension of the forced space, which proves the basis assertion.

The verifier enumerates full ambient genera at the six class-number-one levels where 3 is inert. Good square counts through degree 128 have ranks meeting the incidence bounds. Thus the nonsplit kernel vanishes, the split kernel consists exactly of the opposite-secant relations, and the same finite prefix separates classes across the genera. Exact isometries identify all 31 supplied pairs.

| CM discriminant | Split classes | Nonsplit classes | Opposite-secant pairs |
| ---: | ---: | ---: | ---: |
| $-4$ | 2 | 2 | 0 |
| $-7$ | 3 | 3 | 0 |
| $-19$ | 7 | 6 | 1 |
| $-43$ | 13 | 9 | 4 |
| $-67$ | 19 | 12 | 7 |
| $-163$ | 43 | 24 | 19 |

Each pair lacks norm-one vectors and represents a square with root prime to the lattice discriminant. Its forms take only $0,1\bmod4$. Kani's realization theorem gives principally polarized CM abelian surfaces with these refined Humbert invariants. The absence of norm one makes their polarizations indecomposable, so the surfaces are genus-two Jacobians. The first differing degrees are 3, 6, 9, 12 and 21, occurring in 10, 10, 7, 3 and 1 pairs respectively. Ordinary and primitive counts have the same first difference by divisor inversion. [Kani, realization theorem for primitive CM refined Humbert forms](https://mast.queensu.ca/~kani/papers/refHum5.pdf)

The source ambient block $\mathbb Z\oplus O_D$ contains all nonscalar self-gluings of the given order. It has no exchanged secants at these six levels. The same-genus self-gluing pair differences are orthogonal to the full good-degree kernel, so their observed fraction is **one**. Different-genus gluing pairs separate by the checked prefix.

## An infinite family with an unbounded first difference

Take $f\equiv181\pmod{300}$ and $f>157$. Put $k=(1+4f)/25$ and let $G_f$ be the lattice
$$
 Q_f(x,y,z)=25x^2+4y^2+2xz+kz^2 .
$$
Define $L_0=\{x\equiv0\bmod3\}$ and $L_1=\{x\equiv z\bmod3\}$ inside $G_f$. The parent has determinant $4f$ and belongs to the genus of $\langle1,4,4f\rangle$. In fact
$$
 Q_f=(5x+z/5)^2+4y^2+4f(z/5)^2
$$
is a local change of basis of determinant one at every prime other than 5. At 5 both lattices are unimodular of the same rank and determinant. The sublattices have determinant $144f$. Modulo 3 the form is $y^2+(x+z)^2+z^2$, so the displayed planes are the two split planes through $[0:1:0]$. Directly, on either plane the value of $(x+z)^2+z^2$ is $2z^2$ modulo 3. If a vector has a nonzero square norm modulo 3 then $z=0$ and hence $x=0$ modulo 3. The two lattices therefore contain exactly the same vectors of square norm prime to 3.

The two classes are distinct. The completed square gives $Q_f\ge4fz^2/25>25$ for $z\ne0$. The only norm-four and norm-25 vectors are respectively $\pm(0,1,0)$ and $\pm(1,0,0)$. An automorphism sends each of these pairs to itself. The image of $(0,0,1)$ has third coordinate $\pm1$ and zero second coordinate. Pairing with $(1,0,0)$ forces its first coordinate to vanish and its sign to agree with the sign on $(1,0,0)$. Thus the projective automorphism group is generated by changing the sign of $y$ and fixes both planes individually. The canonical ambient extension excludes an isometry between the two sublattices. The absence of norm one in the parent places the pair outside the nonscalar self-gluing family of discriminant $-4f$, whose parent has norm one.

Both forms are $0,1\bmod4$ by their local structure at 2. Choose $r\equiv16\pmod{30}$ and $r\equiv0\pmod f$, which is possible since $f$ is prime to 30. The vector $\big((r^2-1)/5,r,0\big)$ belongs to both lattices and has norm $(r^2+1)^2$, with $r^2+1$ odd and prime to $3f$. Kani's realization theorem thus gives principally polarized abelian surfaces over $\overline{\mathbb Q}$. Its binary genus description gives discriminant $-36f$, so the CM field is $\mathbb Q(\sqrt{-f})$. The absence of norm one gives Jacobians. Distinct refined Humbert forms imply distinct curves. [Kani, CM refined Humbert forms, Theorems 1, 2 and 20](https://mast.queensu.ca/~kani/papers/refHum5.pdf)

**Theorem 2.** For $M=24t+9$ with $t\ge0$ and $f=(25M^2-101)/4$, the first distinguishing degree is exactly $M$. At that degree $L_0$ has four more signed primitive vectors of norm $M^2$ than $L_1$.

We have $f\equiv181\pmod{300}$ and $k=M^2-4$. The layer $z=0$ contributes equally at every norm. For $|z|\ge2$, the bound $4fz^2/25$ exceeds $M^2$. For $|z|=1$ and $x\ne0$ one has $25x^2+2xz\ge23$, again giving a norm above $M^2$. The remaining norms at most $M^2$ are $M^2-4$ at $y=0$ and $M^2$ at $y=\pm1$. The former lies strictly between $(M-1)^2$ and $M^2$. The four vectors $(0,\pm1,\pm1)$ of norm $M^2$ are primitive and belong only to $L_0$.

There are infinitely many fundamental CM discriminants in Theorem 2. The quadratic polynomial $f=3600t^2+2700t+481$ has nonsquare discriminant $3600\cdot101$. It has no fixed square divisor because its value at zero is squarefree. The classical squarefree-value theorem for quadratics gives infinitely many squarefree positive values. Since $f\equiv1\pmod4$, each gives a fundamental discriminant $-4f$. [Squarefree values of irreducible quadratics](https://research-explorer.ista.ac.at/download/15312/18794/2024_JourNumberTheory_Browning.pdf)

The first difference also grows without bound in the **fixed field** $\mathbb Q(\sqrt{-481})$. The element $U=(964140+43961\sqrt{481})^2$ has norm one and is congruent to 1 modulo 120 in $\mathbb Z[\sqrt{481}]$. Write
$$
 (45+2\sqrt{481})U^j=5M_j+2s_j\sqrt{481},\qquad j\ge0.
$$
Then $M_j\equiv9\pmod{24}$, $s_j\equiv1\pmod{60}$, and $25M_j^2-1924s_j^2=101$. So $f_j=481s_j^2$ satisfies Theorem 2 and $-4f_j$ is the discriminant of the order of conductor $s_j$ in that field. The sequence $M_j$ grows without bound. The two forms now have [uniform bielliptic sextic equations](bielliptic.md) over algebraic number fields. The same construction applies to the entire fixed-field subsequence.

## Uniform exclusion of forced two-class collisions above level three

For the nine class-number-one fields and every odd inert $\ell>3$, let $F$ be the forced kernel in inverse-mass coefficient coordinates and $c=e_i-e_j$ any same-genus pair difference. Then
$$
 \frac{\|P_Fc\|^2}{\|c\|^2}\le\frac23.
$$
In the nonsplit genus $F=0$. In the split genus the [forced diagonal bound](../forced-square-kernel/README.md#concentration-and-the-analytic-question) gives $\rho_F(i)\le2g_D/(\ell+1)$, where $g_D$ is the largest $|O(G)/\{\pm1\}|$ in the ambient genus. Cauchy–Schwarz bounds the pair projection by $\rho_F(i)+\rho_F(j)$. This is at most $1/2$ when $\ell+1\ge8g_D$. Below the threshold the verifier computes the exact rational diagonals at 59 field-prime pairs. The group bound is $g_D\le12$ in every field, and twice the maximum diagonal over the checks is at most $2/3$. Combined with Theorem 1, the 31 pairs exhaust the forced two-class collisions over all odd inert primes in those nine fields.

For the full kernel $V=F\oplus R'$, the observed fraction is
$$
 \operatorname{Vis}(i,j)=1-\frac{\|P_Fc\|^2}{\|c\|^2}
                          -\frac{\|P_{R'}c\|^2}{\|c\|^2}.
$$
The first subtraction is now bounded uniformly. A uniform lower bound for the full observation on self-gluing pairs still requires control of the residual term. The product-generated projector gives an exact expression for the mass norm that survives.

## Reproduction and prior work

~~~sh
python -m pip install cypari2 numpy python-flint
study=research/elliptic-quotients/opposite-secant-spectra
python -B "$study/verify.py" --verify "$study/certificate.json"
~~~

To match the supplied walk census, pass `--archive` with `collision_counts.json` extracted from commit `c754b46fa6ef037248a9fdd63fcca6540392789f`. The all-degree equality is proved by incidence. The first differing degrees, finite rank certificates and exact isometry matching are replayed independently.

Kim and Oh previously gave genera of ternary forms indistinguishable by squares, as well as examples with equality restricted to odd square roots. The results here identify the opposite-secant obstruction at 3 for these Humbert genera and give a genus-two family with an explicitly delayed first difference. [Kim and Oh, square representation numbers](https://arxiv.org/abs/1509.09111)
