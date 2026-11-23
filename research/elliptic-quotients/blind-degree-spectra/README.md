# Every elliptic map degree on the twelve blind quartics

Ten of the twelve preselected plane quartics admit **no elliptic map of any positive degree**, even over an algebraic closure. The remaining two have one optimal elliptic quotient each, of degree two. Their rational degree spectra differ. The original list was selected using row numbers before examining outcomes. The separately supplied Klein positive control is outside these twelve. [Source catalogue](https://math.mit.edu/~drew/genus3curves.html) · [Frozen input](../../../applications/atlas/example.json) · [Original bounded decision](../../../applications/atlas/original/catalogue.json)

In the table $n$ ranges over the positive integers and an elliptic target of the indicated field may vary. An empty spectrum proves that no nonconstant morphism to any elliptic curve exists over that field.

| Sutherland row | Degrees over $\overline{\mathbb Q}$ | Degrees over $\mathbb Q$ | Good prime used |
|---:|---|---|---:|
| 42243 | Empty | Empty | 2 |
| 80491 | Empty | Empty | 5 |
| 55578 | Empty | Empty | 2 |
| 16968 | Empty | Empty | 3 |
| 78269 | Empty | Empty | 2 |
| 36518 | $2n$ | $2n^2$ or $6n^2$ | 5 |
| 54964 | Empty | Empty | 7 |
| 23904 | Empty | Empty | 7 |
| 72418 | Empty | Empty | 3 |
| 52173 | Empty | Empty | 7 |
| 32703 | Empty | Empty | 2 |
| 23640 | $2n$ | $2n^2$ or $4n^2$ | 23 |

Seven absolute exclusions already appeared in the original portable atlas. The three additional rows 54964, 23904 and 52173 had geometric degree-three exclusions in the archived Part XXX scan. Moving to the good prime 7 settles all degrees on each. The two positive rows have different elliptic quotients and different rational isogeny obstructions. [Row 36518 proof](../row-36518/README.md) · [Row 23640 proof](../row-23640/README.md) · [Historical portable results](../../../applications/atlas/receipt.json)

## The three new absolute exclusions

Let $p=7$ and write $N_i=\#C(\mathbb F_{p^i})$. Write the sextic characteristic polynomial of Jacobian Frobenius as
$$
 h(T)=T^6+aT^5+bT^4+cT^3+pbT^2+p^2aT+p^3 .
$$
Exact enumeration and Newton identities give

| Row | $(N_1,N_2,N_3)$ | $(a,b,c)$ |
|---:|---|---|
| 54964 | $(12,54,360)$ | $(4,10,24)$ |
| 23904 | $(11,61,332)$ | $(3,10,17)$ |
| 52173 | $(10,58,424)$ | $(2,6,36)$ |

Each $h$ is irreducible over $\mathbb Q$. To test geometric simplicity, compute $R(U)=\operatorname{Res}_T(h(T),h(UT))$. Its roots are the ratios of Frobenius roots. In each case the primitive factorization of $R$ has the shape $(U-1)^6g_6(U)g_{12}(U)^2$ up to a nonzero rational scalar. The subscript denotes degree. Both $g_6$ and $g_{12}$ are primitive irreducibles in $\mathbb Z[U]$ with leading coefficient 343. A root of unity has a monic integral minimal polynomial, so none of these nontrivial ratios is a root of unity. The [machine-readable certificates](receipt.json) give the full coefficients.

Hence for every $m\ge1$ the six $m$-th powers of the roots of $h$ are distinct and form a single Galois orbit. The Frobenius characteristic polynomial stays irreducible after every finite residue-field extension. The reduced Jacobian is geometrically simple, and a geometric elliptic factor in characteristic zero would specialize after finite extension at this good prime. Any elliptic map would induce such a factor. This proves the three empty spectra. The seven historical absolute exclusions also admit this resultant certificate. [Howe and Zhu on absolute simplicity](https://arxiv.org/abs/math/0002205)

## An exact local coincidence with no global homomorphism

Rows 78269 and 32703 have identical Frobenius polynomials at 2 and hence identical point counts over **every** extension of $\mathbb F_2$. Both curves also have five points over $\mathbb F_3$. Nevertheless their Jacobians have no nonzero homomorphism over $\overline{\mathbb Q}$. This gives F03 and F13 an explicit pair drawn from the fixed blind dataset.

The source counts over the first three extensions of $\mathbb F_3$ are $(5,7,23)$ for row 78269 and $(5,19,35)$ for row 32703. Let $h_1,h_2$ be the corresponding Frobenius sextics. The primitive part of $\operatorname{Res}_T(h_1(T),h_2(UT))$ is irreducible of degree 36 and has leading coefficient $19{,}683$. No root of this cross-resultant can be a root of unity. Thus the two sets of Frobenius roots have disjoint $m$-th powers for every $m\ge1$. A homomorphism after any finite residue-field extension would induce an intertwiner on the rational Tate modules. Their disjoint eigenvalues force that intertwiner to vanish. Specialization at good reduction proves the same Hom vanishing over $\overline{\mathbb Q}$. The receipt includes all 37 primitive cross-resultant coefficients and the common polynomial at 2.

## Replay

From this folder, with Python 3.11 or later and SymPy installed, run

```sh
python verify.py
```

The [verifier](verify.py) repeats the row-number selection from all 82,240 indices, checks the twelve equations against the original frozen inputs, verifies smoothness at every chosen prime, and recounts projective points over all three extension fields for each of the ten absolute exclusions. It factors the exact sextics and their root-ratio resultants. It then reruns and matches both independent bielliptic case receipts. The twelve determinations establish this benchmark's complete degree spectra. A general inference about prevalence on Sutherland’s full dataset would need a different statistical study.
