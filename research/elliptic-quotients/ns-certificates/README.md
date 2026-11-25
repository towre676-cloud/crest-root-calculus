# Exact Néron–Severi certificates: the complete elliptic spectrum of any split genus-two Jacobian

**25 November 2026.** This study removes the remaining restrictions on certified completeness for presence. The [exact-completeness study](../exact-completeness/README.md) handled two isomorphic or non-isogenous targets through one finite kernel. Here the object certified is the Néron–Severi lattice $\mathrm{NS}(J_{\overline{\mathbb Q}})$ itself, with its intersection form. The replay covers non-isogenous factors, distinct targets linked by a degree-three isogeny, and CM factors with $j=1728$. By Kani's description, every elliptic curve in $J$ is a primitive class $D$ with $D^2=0$ and $D\cdot\theta>0$, and its optimal map has degree $D\cdot\theta$. An exact lattice therefore determines every degree and multiplicity.

## The certificate

**Classes from maps.** Let $f:C\to E$ be a certified map of degree $d$. Put $N_f=f^*f_*\in\operatorname{End}(J)$. The symmetric endomorphism $\alpha_f=d-N_f$ corresponds to a class $D_f\in\mathrm{NS}$ with $D_f^2=2\det\alpha_f=0$ and $D_f\cdot\theta=\operatorname{tr}\alpha_f=d$, using analytic traces and determinants. For two maps,
$$
D_f\cdot D_g=\operatorname{tr}\alpha_f\operatorname{tr}\alpha_g-\operatorname{tr}(\alpha_f\alpha_g)=d_fd_g-\deg(f_*g^*),
$$
and $\deg\mu$ for $\mu=f_*g^*:E_g\to E_f$ is the scalar by which $\hat\mu\mu=g_*f^*f_*g^*$ acts on points. Nothing here assumes that $f$ is optimal; a non-optimal $D_f$ is divided automatically by saturation. Pullback $f^*$ and pushforward $f_*$ are computed exactly on the reduction at a good prime. They come from the explicit formulas of Proposition 6.1, through Cantor arithmetic on the odd model, with conjugate point pairs handled in $\mathbb F_q[t]/(u)$. They are self-tested by additivity and $f_*f^*=[d]$. Sums $g=f_1+\hat\rho\circ f_2$ with a global isogeny $\rho$ are also providers: $g^*=f_1^*+f_2^*\rho$ and $g_*=f_{1*}+\hat\rho f_{2*}$. Here $\rho$ and $\hat\rho$ are the exact reductions of PARI's Vélu isogeny and its dual, self-tested by $\hat\rho\rho=[k]$. The Gram matrix is computed at two primes and must agree.

**Rank.** If $D_f\cdot D_g>0$, then $J\sim E_f\times E_g$ and $\rho(J)=2+\operatorname{rank}\operatorname{Hom}(E_f,E_g)$. The witnesses are:
- for $\operatorname{Hom}=0$: unequal traces of $\pi^{120}$ at a good prime;
- for isogeny: $\deg\mu>0$ or an explicit global isogeny;
- for non-CM: two ordinary primes with different Frobenius fields;
- for CM: $j$ in the list of 13 rational CM $j$-invariants.

The Gram rank must equal this upper bound.

**Saturation.** $D\in\ell\,\mathrm{NS}$ iff $\alpha_D\in\ell\operatorname{End}(J)$, because $\theta$ is principal, and that holds iff $\alpha_D$ kills $J[\ell]$, because $[\ell]$ is separable. For a basis vector with denominators, multiply by the common denominator $M$ with $\ell$-part $m$. Then the test is whether $M\alpha_D$ kills $J[\ell m]$. Only primes with $\ell^2\mid\mathrm{disc}$ can divide the index. Each is tested on a Smith basis of $J[\ell^s]$ over a field where it is rational, and found classes are adjoined until none remain.

**Spectrum.** The exact enumerator chooses unimodular coordinates in which $D\cdot\theta=g h$, with $g$ the gcd of the polarization row and $h$ an integer. For each $h>0$ with $gh\le X$, write the remaining coordinates as $z$. The orthogonal complement is negative definite, so after completing the square the equation $D^2=0$ is an affine ellipsoid in $z$. Rational LDL decomposition bounds all but one coordinate by exact integer intervals. The final coordinate is obtained by an exact square-root and divisibility check. A gcd test selects the primitive classes. This computes the complete bounded spectrum without a hard-coded cap.

## A caught error

The first saturation test skipped every combination whose generator coefficients had a denominator $\ell$. On the identity-type $k=3$ gluing this returned the transposition spectrum: two degree-two classes instead of six. The period enumeration and the bielliptic theorem both disagreed with it. The Gram matrix before saturation is the same for every gluing type, so the whole gluing lives in saturation, and the skipped combinations were exactly the missing classes. The corrected test works at level $\ell m$ and recovers discriminant 6 and all six classes.

## Results

| Curve | Factors | $\rho$ and its witness | Saturation | $\mathrm{disc}\,\mathrm{NS}$ | Spectrum (exact) | Independent agreement |
|---|---|---|---|---:|---|---|
| $C$ | non-isogenous | 2; traces $2,-6$ at $13$ | none | $-4$ | $2\,(2)$ | torsion certificate |
| $C'$ | non-isogenous | 2; traces $-6,2$ at $13$ | none | $-16$ | $4\,(2)$ | torsion certificate |
| $T_1$, $T_2$ | isomorphic, non-CM | 3; Frobenius fields $\mathbb Q(i)$, $\mathbb Q(\sqrt{-7})$ | $(D_1+D_3)/2$ on $J[2]$ over $\mathbb F_{17^2}$ | $32$ | $4(4),20(8),52(8),68(8),100(8),116(8),148(8)$ | Theorem 4.1, $q=0$ |
| $k=3$, identity gluing | 3-isogenous, non-CM | 3; explicit Vélu isogeny | two steps, levels 1 and 2 | $6$ | $2(6),14(12),26(12),38(12),\dots$ | bielliptic law and periods |
| $k=3$, three transpositions | 3-isogenous | 3 | one step | $24$ | $2(2),4(4),14(4),26(4),28(8),\dots$ | bielliptic law and periods |
| $k=3$, two 3-cycles | 3-isogenous | 3 | one step | $24$ | $2(2),8(4),14(4),26(4),38(4),\dots$ | bielliptic law |
| CM, $j=1728$, $\beta=(1,0,2)$ and $(1,2,0)$ | isomorphic, CM by $\mathbb Z[i]$ | 4; $j=1728$ in the CM list | three steps at $\ell=2$ | $-16$ | $2(4),3(8),4(2),6(8),7(16),10(16),11(24),\dots$ | periods, all 25 degrees through 40 |

Every $k=3$ spectrum agrees with the bielliptic closed law through degree $150$. For both CM gluings the two composite providers are formed with the target isomorphism and the same isomorphism composed with the CM automorphism. They give Gram rank four. Both curves have eight optimal maps of degree three in their certified spectra. The exact degree-slice calculation gives 11,574 classes in 94 degrees through 150 for each curve. The independent period comparison covers the 25 degrees through 40.

## Replay

Install SymPy, cypari2 and Python-FLINT. Run from this directory.

~~~sh
python -B verify.py --full
python -B verify.py --cm-bound 150 cm_102 cm_120
~~~

The full check starts from the explicit period-map records or from the defined bielliptic equations. It checks the map identities in characteristic zero, recomputes the Gram matrices at two primes, computes a rank bound from checked good reductions or the explicit CM input, and replays torsion saturation. It compares the resulting bases and spectra with the stored receipts. Running without the **--full** option checks the spectra from the stored lattices. The examples exercise rational and Gaussian coefficient fields in the full Néron–Severi replay. The current finite-field engine uses a rational branch point and bounded torsion searches.
