# Equal roots and different rational solutions

Two finite algebraic solution spaces can have the same point counts at **every good prime and every extension degree**, while one has no rational point and the other has two. Thus the complete root datum, including all extension counts and every statistic computed from them, cannot determine rational existence on arbitrary finite parameter spaces. This is an explicit application of the classical Klein-four Brauer relation. [de Smit, §2](https://pub.math.leidenuniv.nl/~smitbde/papers/desmit_factor-equivalence_1996.pdf) · [Bartel and Dokchitser, *Brauer relations in finite groups*](https://doi.org/10.4171/JEMS/563)

## Equations and the all-prime identity

Let

$$
S=\operatorname{Spec}\mathbb Q(\sqrt{13})\ \sqcup\
\operatorname{Spec}\mathbb Q(\sqrt{17})\ \sqcup\
\operatorname{Spec}\mathbb Q(\sqrt{221}),
$$

$$
T=\operatorname{Spec}\mathbb Q\ \sqcup\
\operatorname{Spec}\mathbb Q\ \sqcup\
\operatorname{Spec}\mathbb Q(\sqrt{13},\sqrt{17}).
$$

Here each `Spec` means the solution space attached to the indicated field, and the symbol $\sqcup$ keeps the components separate. Equations for the three components of $S$ are $x^2=13$, $x^2=17$ and $x^2=221$. Equations for $T$ are two separate copies of $x=0$ and

$$
x^4-60x^2+16=0.
$$

The quartic is the minimal polynomial of $\sqrt{13}+\sqrt{17}$. Its roots are the four independent sign choices. It has no rational root. Each quadratic is irreducible over $\mathbb Q$. Consequently $\#S(\mathbb Q)=0$ and $\#T(\mathbb Q)=2$. Both spaces have six geometric points. Both have three rational connected components, with respective field degrees $(2,2,2)$ and $(1,1,4)$.

For a prime $p\notin\{2,13,17\}$, put $a=(13/p)$ and $b=(17/p)$, the quadratic residue signs. Over $\mathbb F_{p^r}$ the counts are

$$
\begin{aligned}
\#S(\mathbb F_{p^r})
  &=(1+a^r)+(1+b^r)+(1+(ab)^r),\\
\#T(\mathbb F_{p^r})
  &=2+(1+a^r)(1+b^r).
\end{aligned}
$$

These expressions agree identically. If $r$ is even the common count is six. If $r$ is odd it is six when $a=b=1$, and two otherwise. This proves equality at every listed prime and every extension degree. It also proves equality of the entire local zeta functions. The models here are disjoint unions, so possible overlap between coordinate values in different quadratic equations never merges their components.

The Galois group is $G=C_2\times C_2$. For its three index-two subgroups $H_i$, the permutation-character identity is

$$
\mathbb Q[G/H_1]\oplus\mathbb Q[G/H_2]\oplus\mathbb Q[G/H_3]
\simeq \mathbb Q^2\oplus\mathbb Q[G].
$$

Both characters are $3\mathbf1+\chi_{13}+\chi_{17}+\chi_{221}$. The two $G$-sets are different, since their orbit sizes are different. Their coordinate algebras with Galois action therefore retain information that the associated linear representations have forgotten. The average fixed-point count is three in both cases. Burnside's counting lemma identifies this average with the number of Galois orbits.

The conclusion concerns rational points and algebra maps. The rational Artin motives of $S$ and $T$ are isomorphic, as the displayed representation identity shows. For abelian varieties over a number field, equality of all good-prime Frobenius polynomials implies isogeny, using Chebotarev and Faltings' theorem. The type of object and the type of map are essential parts of a reconstruction statement. [Kedlaya, §9.4](https://kskedlaya.org/weil-cohom/chapter-9.html)

## A local solution at every completion still leaves a global obstruction

The space $S$ even has a point over every completion of $\mathbb Q$. At an odd prime outside 13 and 17, the three residue signs $a,b,ab$ cannot all equal $-1$. A square unit residue lifts to a square in $\mathbb Q_p$. At 13 the integer 17 is the square residue 4. At 17 the integer 13 is the square residue $8^2$. At 2, the odd unit 17 is congruent to 1 modulo 8 and is a 2-adic square. At the real place all three radicands are positive. This is the classical zero-dimensional counterexample to the Hasse principle. [Colliot-Thélène, *Local-global principle for rational points and zero-cycles*, §3.1](https://swc-math.github.io/aws/2015/2015Colliot-TheleneNotes.pdf)

## Consequence for Crest–Root

The general reconstruction problem must specify what information accompanies the counts. For a finite map-parameter space, an explicit coordinate algebra and its defining equations allow exact rational reconstruction. Counts alone lose the multiplication and the individual Galois orbits. The [finite correspondence algorithm](../bounded-correspondences/README.md) retains that algebra, derives its reconstruction bounds, follows every allowed modular branch and checks the resulting equations over $\mathbb Q$.

This example establishes an obstruction for general finite solution spaces. The positive hyperelliptic map theorem in the companion work has its own geometric hypotheses and supplies an incidence space to which reconstruction applies. The classical Brauer relation and Hasse counterexample are credited above. The project contribution here is their explicit verification and their use to state the reconstruction interface precisely.

## Reproduce the calculation

From a checkout with SymPy installed, run

```sh
python research/local-pattern-obstruction/verify.py
```

The verifier checks the four group elements, rational factor degrees, the quartic minimal polynomial and every exceptional local square test. It independently computes finite-field root counts using $\deg\gcd(f,x^{p^r}-x)$ for 22 primes below 101 and extension degrees one through six. The all-prime and all-degree claims follow from the displayed character identity. [Computed receipt](receipt.json)
