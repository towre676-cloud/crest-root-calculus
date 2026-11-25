# First separating degrees in inert CM gluings

The exact census from the uploaded study separates all inequivalent forms at the 17 reported levels. The proposed bound $2\ell$ for the first separating degree fails within the maximal orders of class number one. With nonmaximal orders, first separating degrees are unbounded even at a fixed inert prime in a fixed CM elliptic isogeny class.

## An audited census at the original 17 levels

The uploaded lattice computation at commit `90728b818df457749e5d0ff3309737ac0a48f144` constructs the Néron–Severi group by requiring integral Riemann forms on the glued lattice. It obtains $q(D)=(D\cdot\theta)^2-2D^2$ by taking the quotient by the canonical polarization. [Kani's correspondence](https://doi.org/10.1007/s40316-018-0105-6) associates each primitive vector of norm $m^2$ with one geometric optimal degree-$m$ class.

The [compact imported census](imported-census.json) retains every representative gluing, its Gram matrix, its genus and every recorded first separating witness. The [audit](audit.py) reconstructs the form for every representative from the uploaded Riemann-form method. For each it verifies an exact integral isometry with the archived Gram matrix and another with the form from the [closed inert-level formula](../inert-gluing/README.md). It checks that the scalar parts cover exactly the predicted orbits, with the correct multiplicity for each. It verifies the archived genus partition by the formula for the Legendre symbol. Exact rational ellipsoid enumeration of the archived forms and direct integer enumeration of the closed forms agree on the first witness for every pair. That gives **463 forms and 5,368 verified first witnesses**. The Gaussian level-two benchmark also reproduces both existing Néron–Severi certificates.

The archived witness counts use one vector per sign. The audit checks that normalization in the receipts and reports both primitive vectors when stating geometric multiplicities. In the Humbert quotient, opposite vectors correspond to the two complementary elliptic subgroups.

## A maximal-order counterexample to the bound $2\ell$

Take $D=-163$, $\ell=5$, and $\omega=(1+\sqrt{-163})/2$. Take scalar parts $a=0$ and $a=\sqrt{-163}$ in $O/5O$. Both are admissible. They have the same genus since $1+N(a)$ is $1$ or $4$ modulo $5$, and their Humbert forms are inequivalent since their scalar orbit representatives differ.

For $c=x+y\omega$, the ambient form is
$$
 u^2+(2x+y)^2+163y^2.
$$
The two lattices impose $u\equiv0$ and $u+2y\equiv0$ modulo $5$, respectively. In the basis $1,(-163+\sqrt{-163})/2$ of the uploaded lattice checker, representative anti-isometries are
$$
 \begin{pmatrix}2&2\\1&3\end{pmatrix},
 \qquad
 \begin{pmatrix}2&3\\4&3\end{pmatrix}.
$$
The two exact enumerations give equal geometric multiplicities for every degree through $41$, then **eight versus zero** at degree $42$. The eight new vectors have $(u,2x+y,y)=(\pm40,\pm1,\pm1)$, with independently chosen signs. All eight satisfy the first congruence and are primitive. The second congruence excludes them.

This refutes the bound proposed from the smaller census while retaining a separating degree for the pair.

## Unbounded separation at fixed inert prime and fixed CM isogeny class

For every positive $f\equiv1\pmod7$, take $O_f=\mathbb Z+fi\mathbb Z$ and $E_f=\mathbb C/O_f$. These elliptic curves are all isogenous to $\mathbb C/\mathbb Z[i]$. The prime $7$ is inert and prime to $f$, so $O_f/7O_f$ is the quadratic field of order $49$. In the basis $1,fi$ of the residue field, use the two gluing matrices
$$
 \psi_0=\begin{pmatrix}1&0\\0&-1\end{pmatrix},
 \qquad
 \psi_1=\begin{pmatrix}1&0\\2&-1\end{pmatrix}.
$$
Both have determinant $-1$. They are nonscalar, with scalar parts $0$ and $fi$ and conjugate-linear parts $1$ and $1+fi$ respectively.

The Riemann-form proof of the closed formula works for this nonmaximal order. In fact $\operatorname{End}(E_f)=O_f$, the product's integral Hermitian forms have entries in $O_f$, and the residue graph still spans the full two-dimensional space over $O_f/7O_f$. The mixed pairings therefore force the Hermitian pullback matrix to be divisible by $7$. The restriction on the graph then gives exactly the trace congruence as for the maximal order. The resulting forms are
$$
 q_0(z,x,y)=(7z)^2+4x^2+4f^2y^2,\qquad
 q_1(z,x,y)=(7z+2y)^2+4x^2+4f^2y^2.
$$
Both polarizations give Jacobians since neither form represents $1$. They belong to the same genus. The ambient lattice is unimodular at $7$, the normal squared lengths reduce to $1$ and $2$, and $2$ is a square modulo $7$. At other primes both forms are locally isometric to the same ambient lattice.

If $d<2f$, any representation of $d^2$ has $y=0$. The two forms then give identical integer equations and identical primitivity tests. At $d=2f$, a vector with $y\ne0$ must have $y=\pm1$, $x=0$ and vanishing first square. This gives two primitive vectors for $q_0$. For $q_1$ the equation $7z+2y=0$ has no integer solution for $y=\pm1$. Consequently the **first separating degree is exactly $2f$**, with a multiplicity difference of two. The audit checks the lattice calculation at $f=1,8,15,36,71$.

The conclusion holds at a single gluing prime and inside a single elliptic isogeny class. A bound on the first separating degree in this setting must account for the conductor or other data that grows with it.

## Extended census

The [extended receipt](extended-census.json) covers every odd inert prime at most $199$ in all nine maximal imaginary quadratic orders of class number one. There are **210 field–prime combinations and 554,482 form classes counted across them**. Every pair of inequivalent forms at a fixed level has different spectra. For pairs in the same genus, the program partitions the representatives successively by their exact primitive counts and ends with singletons. Pairs in different genera have different limiting fractions of degrees coprime to $\ell$ by the [genus detection theorem](../inert-gluing/README.md#the-unmarked-spectrum-determines-the-genus-at-fixed-level).

In **24 levels** at least one pair first separates above $2\ell$. The largest first difference is degree **304**, at $D=-163,\ell=191$. One pair attaining it has scalar parts $53+85\omega$ and $3+185\omega$ and gives zero versus four geometric classes at that degree. The receipt records the cutoff and a pair attaining the largest first difference for each level.

## Replay

~~~sh
python -m pip install cypari2
python -B research/elliptic-quotients/inverse-spectrum/audit.py
python -B research/elliptic-quotients/inverse-spectrum/census.py --prime-limit 199
~~~

The [closed-form census](census.py) needs only the standard library. It enumerates scalar orbits for the nine maximal orders of class number one, sorts forms by genus and partitions each genus by exact degree multiplicities. If a pair remains tied, the reported result gives its cutoff. The routines treat a first difference as a separation certificate and use no modular-form computation for that decision. The [common theta level and Sturm comparison](../inert-gluing/README.md#uniform-level-for-the-square-comparison) apply to any remaining candidate at a fixed maximal-order inert level.

All-degree spectral injectivity for inequivalent forms in the maximal-order inert family is still a conjecture. The counterexample and conductor family concern the position of a difference.
