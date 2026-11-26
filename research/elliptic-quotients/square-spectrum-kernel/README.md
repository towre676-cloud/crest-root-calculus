# Linear relations among elliptic quotient spectra

The elliptic quotient spectrum can distinguish every Humbert isometry class in a census while still losing information about sums of those classes. For the 463 refined Humbert forms at the [17 original inert levels](../inverse-spectrum/README.md), the rational space of all-degree relations has total dimension **202**, summed over the levels. Relations supported within a single genus account for six dimensions. The full theta series are linearly independent at each of the 17 levels.

## A four-form example

Take $D=-3$, $\ell=5$, and $\omega=(1+\sqrt{-3})/2$. Denote a form by the scalar parameter $a$ in $O/5O$ which defines its [congruence lattice](../inert-gluing/README.md). If $S_a(m)$ counts the geometric optimal degree-$m$ classes, then
$$
S_{1+3\omega}(m)+2S_{1+2\omega}(m)
  =S_{0}(m)+2S_{1+4\omega}(m)
  \qquad(m\ge1).
$$
The four forms are pairwise inequivalent. The first two listed in the table belong to one genus, and the last two to the other.

| Scalar part | Coefficient in the relation | Ordinary representations of $12$ |
| --- | ---: | ---: |
| $1+3\omega$ | $1$ | $0$ |
| $0$ | $-1$ | $6$ |
| $1+4\omega$ | $-2$ | $2$ |
| $1+2\omega$ | $2$ | $0$ |

The signed full theta difference is $-10$ at index 12, its first nonzero term. The identity gives two weighted sums with identical quotient counts in every degree. Each genus contains two forms here, and its two individual spectra are distinguishable.

## Certification at each level

Let $r_a(n)=\#\{v\in L_a:q(v)=n\}$ count both signs. Every nonzero vector has a unique positive content, so $r_a(m^2)=\sum_{d\mid m}S_a(d)$. Möbius inversion gives the converse. Thus relations among the spectra are exactly relations among the ordinary square coefficients.

A relation must have zero coefficient sum in each genus. To see this, let $s_+,s_-$ be the two sums. At a fixed inert level the leading constant for $\sum_{m\le X}S_a(m)$ is positive and common to all $a$. Hence $s_++s_-=0$. Applying the relation only at degrees coprime to $\ell$ also gives $\rho_+s_++\rho_-s_-=0$, where
$$
\rho_\pm=\frac{\ell(\ell\pm1)}{\ell^2+2\ell-1}.
$$
The two fractions differ. The sums therefore vanish separately. The [inert genus detection theorem](../inert-gluing/README.md#the-unmarked-spectrum-determines-the-genus-at-fixed-level) supplies the fractions and the common leading constant.

For vectors satisfying those constraints, the theta combination is a sum of within-genus differences and is cuspidal. All forms have the same theta level $M=4|D|\ell^2$ and quadratic character. The $t=1$ Shimura lift is holomorphic of weight two with trivial character at level dividing $M$. Its coefficient at $m$ equals the square coefficient at $m^2$ plus a linear combination of coefficients at proper divisor squares. If the ordinary square coefficients vanish through $m=\lfloor[\mathrm{SL}_2(\mathbb Z):\Gamma_0(M)]/6\rfloor$, the weight-two Sturm bound makes the lift zero. Triangular inversion then gives the relation at every square. Together with the two genus equations this identifies the exact rational kernel. [Shimura mapping in weight $3/2$](https://doi.org/10.1017/S002776300001895X) · [Sturm bound](https://www.lmfdb.org/knowledge/show/cmf.sturm_bound)

The replay uses exact integer enumeration. For $c=x+y\omega$ and $T=\operatorname{Tr}(\omega)$, it counts
$$
 u^2+(2x+Ty)^2+|D|y^2=n,\qquad
 u+\operatorname{Tr}(a)x+\operatorname{Tr}(a\omega)y\equiv0\pmod{\ell}.
$$
It calculates $a$ from each archived anti-isometry, computes a rational nullspace and checks its basis vectors through the full Sturm bound. To test independence of the full theta series it then restricts the ordinary coefficient rows to that nullspace. Full rank there suffices, because any full-theta relation would already be in the square kernel.

| Discriminant | Prime and all-degree kernel dimension |
| --- | --- |
| $-3$ | $5:1,\ 11:6,\ 17:15$ |
| $-4$ | $3:0,\ 7:3,\ 11:8,\ 19:24,\ 23:35$ |
| $-7$ | $3:0,\ 5:2,\ 13:20,\ 17:36$ |
| $-8$ | $5:2,\ 7:5,\ 13:22$ |
| $-11$ | $7:4,\ 13:19$ |

## Recovering the information lost within a genus

The independent counter reproduces the six archived relations in five genera. Adding $r(5)$ and $r(8)$ makes the coefficient map injective on the span of forms in every individual genus in the census. The pair $r(5),r(20)$ also works, so the additional coefficients may both be taken from the single squarefree class $5$. These coefficients count lattice vectors of non-square norm. The recovery assertion concerns the ordinary theta series of combinations in a fixed genus. Across the two genera, the kernel dimensions can be much larger.

The relevant central-value question includes a quadratic twist. In the convention $q(v)=v^{\mathsf T}Gv$, the theta character is the Kronecker symbol $(2\det(2G)/\cdot)$. Here $\det G=4|D|\ell^2$, so the character on level units is $(|D|/\cdot)$. Waldspurger's twist in square class $t=1$ is the primitive character of $(|D|/\cdot)(-1/\cdot)=(D/\cdot)$. Identifying the kernel by central values requires studying the corresponding twisted values and the local coefficient factors at the bad primes. The rank computations establish the relations independently of such an identification. [Ternary theta character](https://www.impan.pl/shop/en/publication/transaction/download/product/92199) · [Waldspurger formula and local factors](https://math.berkeley.edu/~wodzicki/prace/Sza.pdf)

## Reproduction

The [individual-genus replay](replay.py) checks the archived ranks and [kernel basis](relations.json) from commit `f340eba16bdd9760739ef6d66a1db74dad045e0b`. The [joint replay](joint.py) uses the same [integer counter](count.cpp) and publishes a finite certificate for every kernel dimension. Both use the [archived representatives](../inverse-spectrum/imported-census.json).

~~~sh
python -m pip install cypari2 sympy
python -B research/elliptic-quotients/square-spectrum-kernel/replay.py
python -B research/elliptic-quotients/square-spectrum-kernel/joint.py
~~~

The verification outputs are [individual genera](within-receipt.json) and [fixed levels](joint-receipt.json).

The [neighbour and oldform study](../hecke-observation/README.md) gives a second separation method using only prime-power quotient degrees. It also computes the Shimura lifts on two eigensystems and determines which square coefficients distinguish copies of the same weight-two newform.
