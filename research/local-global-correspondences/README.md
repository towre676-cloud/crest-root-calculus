# From local Frobenius to global genus-two correspondences

Two integral binary sextics suffice to recover every rational polarized Richelot isogeny between their genus-two Jacobians, up to sign. The solver returns a projective change of coordinates, a cubic algebra describing the branch pairing and two equations for a correspondence. It permits Galois to permute all three pairs. On the two benchmarks below it also proves that the displayed map generates the *entire* geometric Hom group. A separate all-prime argument determines every rational isogeny from each source to any genus-two Jacobian. A third source has a four-vertex rational isogeny class with explicit primitive maps of degree sixteen. [Discovery](discover.py) · [Reproducible certificates](receipt.json)

Here a class of isogenies between principally polarized Jacobians is taken up to isomorphism of the target. An isogeny is **primitive** if it is not an integer multiple by an absolute value greater than one in the geometric Hom lattice. The degree of an isogeny is its finite kernel order. The correspondences below induce isogenies of degree four between two-dimensional Jacobians. Their curve graphs have bidegree $(2,2)$.

| Source | Geometric Hom to the displayed target | All rational isogeny degrees to any genus-two Jacobian | Number of rational classes |
|---|---|---|---|
| $y^2=x^6+4x^4-2x^3+7x^2-4x+5$ | $\mathbb Z\phi$, with $\deg\phi=4$ | $n^4$ and $4n^4$, $n\ge1$ | One in each occurring degree |
| LMFDB 88343.a.88343.1, $y^2=5x^6+4x^5+10x^4+2x^3+5x^2-2x+1$ | $\mathbb Z\phi$, with $\deg\phi=4$ | $n^4$ and $4n^4$, $n\ge1$ | One in each occurring degree |

For the first source the rational Richelot target in its reconstructed coordinates is $v^2=z^6-8z^4-2z^3+13z^2-10z-1$. The [input target](example.json) instead has all seven coefficients, a change $z=(3w+2)/(w+1)$ and a different vertical scale. The program recovers them. A quadratic twist of that target gives the same geometric Hom lattice, with generator defined over $\mathbb Q(\sqrt2)$ and zero rational Hom. A different test adds $11x^2$ to the *source* sextic. The two resulting curves agree over every extension at 11. At 13 their prime-field counts still agree, but the exponent-two counts are 152 and 162. A cross-root resultant proves that their geometric Jacobian Hom group vanishes. [Exact reductions and graphs](receipt.json)

The public source is an independent catalog input. Completing the square in its recorded model gives the sextic above. From its coefficients the solver constructs the rational target $v^2=z^6+16z^5+46z^4+34z^3-19z^2-2z+1$. Reductions at 17 and 19 give different real quadratic fields in the endomorphism test. The all-prime proof below uses point counts on the equations. The catalog independently records that the sole nonmaximal mod-$\ell$ image occurs at 2. [Catalog source](https://www.lmfdb.org/Genus2Curve/Q/88343/a/88343/1) · [Public input](public.json) · [Banwait et al., computed nonsurjective primes](https://doi.org/10.1090/conm/796/16000)

## Four rational targets and degree-16 correspondences

Starting from $A_0=J(y^2=x^6+67x^4-2x^3+1588x^2+140x+13249)$, the same algorithms find precisely three rational Richelot neighbors $A_1,A_2,A_3$. Each neighbor has a single edge back to $A_0$. The four branch sextics are irreducible over $\mathbb Q$. The source's geometric endomorphism ring is $\mathbb Z$ by reductions at 13 and 29, and the finite all-prime test excludes every odd prime from every primitive rational multiplier.

| Primitive rational isogeny degree | $A_0$ | $A_1$ | $A_2$ | $A_3$ |
|---|---:|---:|---:|---:|
| From $A_0$ | 1 | 4 | 4 | 4 |
| From $A_1$ | 4 | 1 | 16 | 16 |
| From $A_2$ | 4 | 16 | 1 | 16 |
| From $A_3$ | 4 | 16 | 16 | 1 |

Multiply an entry by $n^4$, for any positive integer $n$, to obtain every rational isogeny degree for that ordered pair. From any leaf, degree sixteen has three target classes. Two come from composing the adjoint edge to the center with an edge to a different leaf. The third comes from multiplication by two on the starting Jacobian. The receipt gives the four integral equations and three edge correspondences, including the coordinate transport and graph identities. Those graphs and their transposes specify both degree-sixteen compositions on each leaf. [Input](four-targets.json) · [Component search](rational_component.py) · [Models and correspondence equations](receipt.json) · [Distance and primitivity proof](arithmetic.md#rational-two-power-kernels-and-rational-components)

## All geometric quotient fields and all degrees

The equation-only atlas reconstructs **all fifteen geometric Richelot quotient models**. In each benchmark the Galois orbits of kernels have sizes **1, 6 and 8**. The certificate gives the model over each of the three fields. Its conjugates give all fifteen. On the first curve the two nontrivial fields have defining polynomials

$$
u^6+3u^5-35u^4-219u^3+106u^2-48u+8,\qquad
u^8+10u^7+10u^6+56u^5-555u^4+1662u^3+2888u^2+1048u+64.
$$

Two different kernels give nonisomorphic geometric targets here. Indeed the source has geometric endomorphism ring $\mathbb Z$. Any isogeny identifies its rational Hom to a given target with a one-dimensional space. Two maps to isomorphic targets of the same degree differ only by sign and therefore have the same kernel. The kernel stabilizer is consequently also the target's moduli stabilizer. The number-field models realize this minimal field of definition. [Orbit and model algorithm](kernel_fields.py) · [Independent number-field checker](check_orbits.py)

There is a count for **every** geometric isogeny degree, with the target allowed to vary. Write $m$ for the polarization multiplier. The isogeny degree is $m^2$. For $m=\ell^r$ the number of *primitive* targets is
$$
P_{\ell,1}=(\ell+1)(\ell^2+1),\qquad
P_{\ell,r}=(\ell+1)(\ell^2+1)
\left(\ell^{3r-3}+\ell^{2r-3}\frac{\ell^{r-1}-1}{\ell-1}\right)\quad(r\ge2).
$$
The total at multiplier $\ell^r$ is $N_{\ell,r}=P_{\ell,r}+N_{\ell,r-2}$ with $N_{\ell,0}=1$, $N_{\ell,1}=P_{\ell,1}$. Counts multiply across prime powers. For example there are 150 primitive targets and 151 total in degree 16, then 1,320 primitive and 1,335 total in degree 64. These statements count targets up to geometric isomorphism. The degree-four atlas also gives their equations. [Counting proof](arithmetic.md#counting-all-geometric-targets) · [Calculator](spectrum.py) · [Independent Hermite-lattice enumeration](lattice_check.py)

## Infinite family

Let $u^3-u-1=0$ and
$$
C_a:\quad y^2=\operatorname{Norm}_{\mathbb Q(u)/\mathbb Q}(x^2+ux+u^2+a).
$$
Every integer $a\equiv1\pmod{493}$ gives a smooth geometrically simple Jacobian with endomorphism ring $\mathbb Z$, hence the all-degree geometric counting formula. At $a\equiv1\pmod{1479}$ its only primitive rational two-power isogenies to genus-two Jacobians have degrees 1 and 4. At $a\equiv1\pmod{29087}$ the full geometric two-torsion pairing orbits have sizes 1, 6 and 8. On the progression
$$
a\equiv1\pmod{30\,557\,478\,886\,935}
$$
the **complete rational spectrum** consists of the degrees $n^4$ and $4n^4$, with one class in each. All four conclusions hold on the last progression. It contains infinitely many geometric isomorphism classes. Every rational quadratic twist of any source with the all-prime conclusion has the same degrees, with the corresponding twists as its two rational targets. [Family, semistability and all-prime proof](arithmetic.md) · [Symbolic family verification](family.py)

## Reproduction

Install SymPy and NumPy. From this directory run

```sh
python discover.py example.json
python kernel_fields.py public.json --models
python all_primes.py example.json
python rational_component.py four-targets.json
python spectrum.py 2 3 4 6 8
python verify.py
```

A new input has integral coefficient arrays named `source` and `target` in increasing powers of the variable. It accepts smooth quintics or sextics. The sextic atlas also works with the `source` entry alone. Its certificates cover product quotients by detecting the degenerate determinant, and give equations for the genus-two quotients. The optional `python verify.py --extended` separately checks every odd prime below ten million using Frobenius irreducibility. The finite all-prime proof already runs in the default replay. Worked controls include a fully split sextic with fifteen rational kernels, a quintic with orbit degree fifteen, a symmetric sextic where the generic covariant frame degenerates, and the quadratic twist and local collision above.

The genus-two Richelot recognition theorem is a fixed-polarization result. The complete Hom and all-degree conclusions use additional checked hypotheses. On other inputs the programs distinguish an inconclusive upper bound from an actual exclusion. The computations specialize the general isogeny-class problem studied by van Bommel, Chidambaram, Costa and Kieffer to exact certificates for these equations and progression. General correspondence discovery between arbitrary motives remains open. [Research context](https://arxiv.org/abs/2301.10118) [Construction, descent and completeness](construction.md) · [Arithmetic certificates](arithmetic.md)
