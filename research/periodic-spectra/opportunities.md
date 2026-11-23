# Where Crest–Root can compete in materials and chemistry

## The question that the experiment answers

The strongest immediate target is an exact calculation on a structured
periodic model: recover the complete characteristic polynomial of a finite
Jacobi ring, retaining every coefficient and every algebraic multiplicity.
The [proof](proof.md) specifies that task. The measurements in this directory
determine how the implementation compares with PARI/GP on that task; they
do not establish superiority throughout crystallography or computational
chemistry.

The proposed strength is a precise one: recognize a repeated algebraic
structure, reduce the calculation to a small cell, assemble the global
answer, and verify the result against the original equations. Exact
arithmetic makes those identities checkable. It is not itself a distinction
from PARI, which already provides exact integer, rational, polynomial, and
number-field arithmetic [1].

## Which applications are a good fit?

| Application | Concrete mathematical task | Present assessment | Necessary comparator |
|---|---|---|---|
| Periodic structures | Exact spectrum of a finite repeated coupling operator | Best immediate fit; the Jacobi-ring calculation is a completed, narrow instance | PARI general characteristic polynomial and PARI with the same structural reduction |
| Quadratic energy landscapes | Characteristic polynomial and, in a further implementation, exact inertia of the Hessian | Direct interpretation for the synthetic spring models; inertia is a separate output task | Exact determinant/root-count methods and a structure-aware Hessian method |
| Reaction networks | Stability of a supplied positive steady state in a periodic diffusion network | Strong next application; multiple species require a block, generally nonsymmetric operator | Published Fourier-mode analysis, an equivalent GP implementation, and appropriate numerical eigenvalue methods under their own accuracy contract |
| Nonlinear polynomial constraints | Enumerate all positive steady states or all stationary points within a stated domain | Related to the finite-algebra reconstruction program, but not solved by this spectral experiment | Singular or other polynomial-system methods, with identical positivity, completeness, and multiplicity requirements |
| Crystal lattice reduction and symmetry | Reduced bases, lattice isometries, and complete automorphism groups | No comparative result established here | PARI's lattice routines; spglib for symmetry of supplied crystal coordinates |
| Discrete material ordering | Enumerate inequivalent atom or vacancy arrangements at fixed composition | Relevant symmetry mathematics, with mature specialist competitors | Hart–Forcade methods, enumlib, and its Enumlib.jl successor |
| Molecular enumeration | Enumerate constitutional isomers of a formula under valence constraints | No comparative result established here | Surge, MOLGEN, or MAYGEN; a handwritten GP enumerator alone is insufficient |
| General materials optimization | Find a physically meaningful optimal structure, including its property model | Requires an explicit objective, feasible set, and physical evaluation beyond exact spectral algebra | The appropriate optimization and materials-modeling workflow |

PARI already exposes lattice reduction through `qflll`, short-vector
enumeration through `qfminim`, and lattice isometry and automorphism
calculations through `qfisom` and `qfauto` [2]. Crystal symmetry from atomic
coordinates is a different problem: spglib supplies operations, space-group
identification, Wyckoff assignments, and primitive-cell searches [3]. A
benchmark of one problem cannot establish an advantage on the other.

For derivative structures, Hart and Forcade already use the quotient group
of a superlattice to enumerate inequivalent atomic configurations [4].
Their enumlib implementation and subsequent extensions are the substantive
competition. Its current README also identifies **Enumlib.jl** as a Julia
successor implementing the same algorithms, so a new comparison should
include that implementation where it supports the requested task. The
Fortran code retains displacement-direction and two-dimensional surface
enumeration; the successor also documents a symmetry-convention difference
for some overlapping sublattice label sets. Freeze those conventions and
compare equivalent structures, not just totals [5]. For molecular
structures, Surge uses canonical generation
and graph automorphisms to enumerate nonisomorphic constitutional isomers
[6]. Any future Crest claim in either area must face those algorithms on
the same output, including the difference between counting structures and
actually emitting them.

## What the spring model means physically

Let $u_i$ be the scalar displacement of site $i$ on a ring of $L$ identical
unit masses. Let $k_i>0$ be the stiffness of the spring joining sites $i$
and $i+1$, and let $p_i\geq0$ be a pinning stiffness. Indices are taken
modulo $L$. The harmonic energy is

$$
E(u)=\frac12\sum_{i=0}^{L-1}k_i(u_{i+1}-u_i)^2
      +\frac12\sum_{i=0}^{L-1}p_i u_i^2.
$$

Its Hessian has diagonal entries $k_{i-1}+k_i+p_i$ and neighboring entries
$-k_i$. It is therefore precisely a Jacobi ring. The equations
$\ddot u=-Au$ identify the eigenvalues of its Hessian $A$ with squared
normal-mode frequencies. Integer stiffnesses supply a mathematically
exact model for the present experiment. The examples are synthetic; no
measured crystal or molecule supplies their coefficients.

This connection has a substantial classical foundation. Aghaei, Dayal,
and Elliott reduce highly symmetric nanotube Hessians to block-circulant
form and then to Fourier blocks [7]. Sharma and Suryanarayana likewise
derive the periodic force-constant and mass-matrix description of phonons
within density-functional perturbation theory [8]. These papers motivate
the application. They do not turn a scalar spring-ring benchmark into a
first-principles phonon calculation.

A material calculation must obtain force constants from a specified
physical model, account for its units and uncertainty, and generally
handle several displacement coordinates per atom and multiple atomic
masses. With a positive diagonal mass matrix $M$, the relevant equation is
$Kv=\omega^2Mv$. Rational masses can be treated through the generalized
pencil $\det(K-xM)$ with explicit normalization; there is no need to
introduce irrational mass square roots merely to formulate an exact task.
That extension is outside this experiment's integral, unit-mass contract.

Exact arithmetic resolves roundoff within a supplied exact model. It does
not determine experimental force constants more accurately, validate a
potential-energy approximation, or make approximate atomic positions
exact. A future physical-data experiment must state how those inputs are
represented and what conclusions survive their uncertainty.

## What makes the PARI comparison meaningful

The general baseline is the default `charpoly` on the assembled integer
matrix. PARI selects a modular algorithm for this input type; replacing it
by an intentionally slow determinant expansion would weaken the test [2].
The other baselines implement the same periodic reduction in GP, with
linear recurrence and doubling variants. They test whether an observed
advantage persists after PARI is given the same structural information.

A separate [compiled libpari control](native-results.md) implements both
linear and doubling variants as native executables and reruns Crest on
the same 20 expanded inputs. Each timed solve includes validation,
period discovery, and the complete polynomial. **Native PARI doubling
wins the computation comparison on 19 of 20 inputs, including all 12
periodic spring cases.** This post-study control exposes the GP-language
cost in the original comparison; it supports pursuing structure recovery
and exact certificates, while providing no general arithmetic-speed
advantage over PARI. The report keeps computation and cold-process costs
separate.

All methods must return the complete expanded integer coefficient vector.
An unexpanded product, a few traces, approximate roots, or only a stability
decision is a different contract. The frozen protocol includes period
discovery and keeps solver time distinct from fresh-process overhead.
Broken-period cases reveal the effect of losing repetition. The underlying
transfer identity remains valid with one full-length cell, but a defect
must not be averaged away or silently treated as numerical noise.

The benchmark uses an optimized build of **PARI/GP 2.19.1** from the official
source distribution [9], with GMP 6.3.0 and a single threading engine. The
source archive is

<https://pari.math.u-bordeaux.fr/pub/pari/unix/pari-2.19.1.tar.gz>.

Its SHA-256, checked against the official download page, is

```text
419553b11ff8c3ecb8e9031b7d91d7cee003a5d721e685fcb499565c96fe2943
```

The upstream characteristic-polynomial and resultant tests passed for the
build used here. Full compiler and binary provenance belongs to the
benchmark receipt. The classical transfer/Floquet and Chebyshev identities
are credited in the proof; the tested contribution is their recognition,
implementation, verification, and measured behavior.

## A stronger algebraic target: discover the hidden spectral quotient

Write $D(x)$ for the cell trace polynomial and $B\ne0$ for the signed
product of its bonds, as in the proof. The normalized monodromy has
determinant one and trace $D/B$. If $z$ is a Floquet multiplier, then

$$
z^2-\frac{D(x)}{B}z+1=0.
$$

Eliminating $z$ gives the spectral curve

$$
Y^2=D(x)^2-4B^2,\qquad Y=B(z-z^{-1}).
$$

Periodic Jacobi discriminants and Floquet theory are established subjects
[10,11]. Their hyperelliptic description nevertheless makes the
connection to the existing Crest–Root map program exact. A useful next
problem is to discover a lower-degree spectral presentation from its
coefficients, without being handed the cell that produced it.

### Separate normalization from quotient recovery

Repeating a known cell creates Chebyshev composition and square factors.
The [normalization identity](spectral-curve.md) removes those factors and
gives a **degree-one birational identification of the same smooth curve**.
This follows from classical Chebyshev/Floquet theory [15]; it is neither
a new map-discovery theorem nor evidence for a map to a different,
lower-genus curve. Repeated-cell cases are useful normalization controls,
but they cannot alone establish success at quotient recovery.

For an expanded input $P\in\mathbb Q[x]$ having at least one root of odd
multiplicity, first return

$$
P=A^2F,\qquad A\in\mathbb Q[x]\text{ monic},\qquad
F\in\mathbb Q[x]\text{ squarefree}.
$$

Retain the rational leading unit in $F$. For each irreducible factor,
its multiplicity in $F$ is its multiplicity in $P$ modulo two. The
radical $P/\gcd(P,P')$ keeps every distinct factor and is generally the
wrong answer. Normalization needs no full irreducible factorization:
gcd-based squarefree decomposition supplies the required multiplicities.
Compare this stage with that algorithm in PARI and a dedicated routine
such as FLINT's integer squarefree factorization, accounting for rational
denominators and content [16]. Report its cost separately.

### Freeze a coefficient-only discovery contract

For a first implementable challenge, return every normalized polynomial
presentation

$$
F(x)=G(U(x)),\qquad U\in\mathbb Q[x]\text{ monic},\quad U(0)=0,
\quad d=\deg U\geq2,
$$

with squarefree $G$ and $1\leq g(G)<g(F)$, where
$g(H)=\lfloor(\deg H-1)/2\rfloor$ for these squarefree radicands.
Sources of genus below two cannot meet this target requirement.
It certifies the degree-$d$ map
$(x,y)\mapsto(U(x),y)$ from $y^2=F(x)$ to $v^2=G(u)$.
The inner normalization removes affine changes of the target coordinate;
a constant vertical scale can be absorbed into the newly discovered
$G$. A later prescribed-target task must instead track rational versus
geometric vertical multipliers and twists explicitly.

This contract concerns polynomial horizontal maps with constant vertical
multiplier over $\mathbb Q$. It does not certify the absence of rational
horizontal maps, nonconstant vertical multipliers, or maps defined only
over larger fields. In particular, the automatic presentation
$D^2-4B^2=(t^2-4B^2)\circ D$ has a genus-zero target and does not meet
the positive-genus target requirement. Nor does this task claim recovery
of a unique physical Jacobi cell from a spectrum.

The solver receives only the expanded coefficients of $P$ and the public
contract. Cell data, repeat count, bond product, target polynomial, map
degree, seed, and generator identity belong in a verifier-only sidecar.
A separate supplied-cell track may give metadata to every method, but
its results cannot stand in for coefficient-only discovery. Freeze the
algorithms before generating held-out seeds, degrees, coefficient sizes,
and an additional generating family. Neutral input identifiers must not
encode the answers.

Include actual Jacobi-derived polynomials with recorded physical-model
provenance, and label general algebraic composites separately. The panel
must distinguish repeated-cell normalization cases from squarefree
positive-genus quotient cases. Include square-factor inflation, rational
affine coordinate changes, nonminimal periods, and one-coefficient
perturbations. A perturbation is a prospective control, not automatically
a negative example: certify its actual outcome in the declared map class.

### Use the strongest matching decomposition baseline

Classical polynomial decomposition already addresses $F=G\circ U$;
Kozen and Landau give a polynomial-time algorithm in characteristic zero
[17]. The benchmark must therefore compare with the repository's
[corrected dense recurrence](https://github.com/towre676-cloud/crest-root-calculus/blob/main/research/solver-comparison/sympy-decomposition-regression.md)
and an independent implementation such as Singular's `decomp.lib` [18].
The existing [map adapter](https://github.com/towre676-cloud/crest-root-calculus/blob/main/research/solver-comparison/map_worker.py) receives a
target and infers an inner degree. For this challenge, its corrected
candidate recovery must instead scan **all admissible degree divisors**,
recover the outer polynomial by exact left division, restore the leading
unit, and verify each composition. No comparator receives the hidden
target or degree. A single complete decomposition chain is not a list of
all admissible quotient presentations; alternative degree divisors must
also be checked. Apply the same rule to Singular's fixed-degree routine.
For a claim against PARI, also run the same classical degree scan in GP
or libpari; comparison only with a generic elimination script is
insufficient.

Retain the documented stock SymPy regression as a correctness control.
An incorrect negative answer earns no speed advantage. Detection,
normalization, degree search, certificate construction, and independent
verification all need explicit cost accounting. Recovery within this
class is an implementation challenge against established mathematics;
novelty would need a separately demonstrated extension or advantage.

### Acceptance checklist for the next pass

- **Input isolation:** identical expanded rational coefficients for every
  method; no generator metadata in the discovery input; frozen code and
  independently generated held-out cases.
- **Exact normalization:** verify $P=A^2F$, $\gcd(F,F')=1$, and retention
  of the rational unit; compute genus from $F$, never the inflated degree.
- **Map certificate:** verify $F=G(U)$ coefficient by coefficient, the
  squarefreeness and genera of both models, and the degree of the induced
  function-field extension. Keep degree-one normalization separate.
- **Completeness within scope:** examine every admissible degree and
  record canonical answers or exact rejection certificates. Timeouts and
  bounded searches remain unresolved cases, not nonexistence proofs.
- **Fair controls and baselines:** include normalization-only positives,
  genuine quotient positives, certified negatives, and transformed and
  perturbed controls; compare corrected decomposition and dedicated
  squarefree routines under the same information and output contract.
- **Useful advantage:** report correctness and complete costs by cohort;
  require any claimed discovery advantage to survive the held-out panel
  and the best correct classical comparator.

## A concrete chemical challenge: periodic Brusselator stability

The next reaction-network task should start with a published model rather
than relabeling an arbitrary polynomial as chemistry. The Brusselator is a
canonical autocatalytic mass-action model originating with Prigogine and
Lefever [12], with archived executable equations available from the
Physiome project [13]. It is a model reaction system, not a fitted
pharmaceutical prediction.

To specify a reproducible next target, take the dimensionless kinetics

$$
f(u,v)=\alpha-(\beta+1)u+u^2v,
\qquad g(u,v)=\beta u-u^2v,
\qquad \alpha,\beta>0.
$$

The positive homogeneous equilibrium is
$(u_*,v_*)=(\alpha,\beta/\alpha)$. Its kinetic Jacobian is

$$
J=\begin{pmatrix}\beta-1&\alpha^2\\-\beta&-\alpha^2\end{pmatrix}.
$$

For a ring of compartments, let $L_C=2I-S-S^{-1}$ be the graph Laplacian
and let $D_m=\operatorname{diag}(d_u,d_v)$ contain positive diffusion
rates. The full linearization, using site-major ordering, is

$$
\mathcal J=I\otimes J-L_C\otimes D_m.
$$

For each Laplacian eigenvalue $\lambda\geq0$, the relevant block is
$J-\lambda D_m$. This is generally nonsymmetric, so its stability concerns
the real parts of eigenvalues, not the signs of the eigenvalues of a
symmetric Hessian. For a real two-by-two matrix, strict asymptotic
stability is equivalent to negative trace and positive determinant.

An exact rational fixture is

$$
\alpha=1,\quad\beta=\frac32,\quad d_u=1,\quad d_v=20.
$$

Direct substitution gives

$$
\operatorname{tr}(J-\lambda D_m)=-\frac12-21\lambda<0,
\qquad
\det(J-\lambda D_m)=20\lambda^2-9\lambda+1
                    =(1-4\lambda)(1-5\lambda).
$$

Thus the homogeneous kinetic system is stable, while a discrete diffusion
mode is unstable exactly when $1/5<\lambda<1/4$; an endpoint gives a zero
eigenvalue and must be classified separately. The graph-size question is
now precise: determine, with multiplicities, which eigenvalues
$2-2\cos(2\pi j/n)$ lie in this open interval. This reduction is an
algebraic specification of the next challenge, not a claim that its
implementation or benchmark is included in the present study.

The first acceptance gate is an independently verified reaction list and
Jacobian at a positive equilibrium, followed by exact stable, unstable,
and threshold fixtures. The second is agreement between the assembled
$2n$-dimensional operator and the mode calculation, including every
multiplicity. The third is a comparison with equivalent GP mode analysis
using exact algebraic root isolation. A boundary-mode decision must not be
made by a floating-point tolerance. Full-polynomial timings and
stability-only timings must be reported as separate tasks.

To establish relevance beyond a designed fixture, the next pass must also
include a published matrix/model instance with its original boundary
conditions and provenance. NIST's BRUSSEL collection supplies chemical
engineering Jacobians with **Dirichlet**, not periodic, boundary
conditions [14]. Those matrices are useful external controls only if the
boundary distinction is preserved and an appropriate method handles or
explicitly rejects them. They must not be converted into periodic inputs
and then described as unchanged published data.

## What would constitute a durable strength?

A useful result is an exact, complete answer obtained faster because the
implementation finds and exploits a reduction that matters for the
supplied task. A stronger result survives a comparator that knows the
same mathematics, parameter changes, defects, and externally sourced
instances. For real material or chemical decisions, the relevant model
and output also have to be validated. These are distinct acceptance
criteria: success on the exact periodic-ring calculation is valuable
without asserting that the later criteria have already been met.

## Primary sources

[1] PARI/GP project, scope and arithmetic software:
<https://pari.math.u-bordeaux.fr/>.

[2] PARI/GP reference manual, linear algebra and lattice routines:
<https://pari.math.u-bordeaux.fr/dochtml/html/Linear_algebra__Vectors__matrices_and_sets.html>.
Polynomial resultants and guaranteed exact-input root counting:
<https://pari.math.u-bordeaux.fr/dochtml/html/Polynomials_and_power_series.html>.

[3] Atsushi Togo, Kohei Shinohara, and Isao Tanaka, *Spglib: a software
library for crystal symmetry search* (2024); official documentation and
coordinate-tolerance definition:
<https://spglib.readthedocs.io/en/stable/>,
<https://spglib.readthedocs.io/en/stable/variable.html>.

[4] Gus L. W. Hart and Rodney W. Forcade, *Algorithm for generating
derivative structures*, Physical Review B **77**, 224115 (2008),
<https://doi.org/10.1103/PhysRevB.77.224115>;
author-hosted paper:
<https://bsg.byu.edu/docs/papers/GLWHart-enumeration.pdf>.

[5] Hart–Forcade and collaborators, enumlib source and related papers:
<https://github.com/msg-byu/enumlib>.
The current README identifies the Julia successor and its scope;
successor source and documented convention differences:
<https://github.com/glwhart/Enumlib.jl>. Consulted November 23, 2026;
no performance claim from either README is imported here.

[6] Brendan D. McKay, Mehmet Aziz Yirik, and Christoph Steinbeck, *Surge:
a fast open-source chemical graph generator*, Journal of Cheminformatics
**14**, 24 (2022); official project and source:
<https://structuregenerator.github.io/>,
<https://github.com/StructureGenerator/surge>.

[7] Amin Aghaei, Kaushik Dayal, and Ryan S. Elliott, *Symmetry-adapted phonon
analysis of nanotubes*, Journal of the Mechanics and Physics of Solids
**61**(2), 557–578 (2013),
<https://doi.org/10.1016/j.jmps.2012.09.008>;
university publication record:
<https://experts.umn.edu/en/publications/symmetry-adapted-phonon-analysis-of-nanotubes/>.

[8] Abhiraj Sharma and Phanish Suryanarayana, *Calculation of phonons in
real-space density functional theory*, Physical Review E **108**, 045302
(2023), <https://doi.org/10.1103/PhysRevE.108.045302>;
accepted manuscript:
<https://link.aps.org/accepted/10.1103/PhysRevE.108.045302>.

[9] PARI/GP official download page, stable source and checksum, consulted
November 23, 2026: <https://pari.math.u-bordeaux.fr/download.html>.

[10] Günter Stolz and Rudi Weikard, *Seminar on Jacobi Matrices* (2004),
Section 2.4:
<https://www.uab.edu/cas/mathematics/images/Documents/Weikard-Jacobi_matrices.pdf>.

[11] Vassilis G. Papanicolaou, *Periodic Jacobi Operators with Complex
Coefficients* (2019): <https://arxiv.org/abs/1909.09206>.

[12] Ilya Prigogine and René Lefever, *Symmetry Breaking Instabilities in
Dissipative Systems. II*, Journal of Chemical Physics **48**, 1695–1700
(1968), <https://doi.org/10.1063/1.1668896>.

[13] Physiome/Interagency Modeling and Analysis Group, Brusselator model
0329, including MML and SBML downloads:
<https://www.imagwiki.nibib.nih.gov/physiome/jsim/models/webmodel/NSR/brusselator>.

[14] NIST Matrix Market, BRUSSEL reaction–diffusion model and matrix
collection:
<https://math.nist.gov/MatrixMarket/data/NEP/brussel/brussel.html>.

[15] NIST Digital Library of Mathematical Functions, equations
18.5.1–18.5.2, classical Chebyshev representations:
<https://dlmf.nist.gov/18.5>.
Gerald Teschl, *Jacobi Operators and Completely Integrable Nonlinear
Lattices*, Mathematical Surveys and Monographs **72**, AMS (2000),
Chapter 7 and finite-gap spectral curves, author-hosted edition:
<https://www.mat.univie.ac.at/~gerald/ftp/book-jac/jacop.pdf>.

[16] FLINT project, squarefree factorization algorithms and
`fmpz_poly_factor_squarefree`:
<https://flintlib.org/doc/algorithms/polynomials/squarefree.html>,
<https://flintlib.org/doc/fmpz_poly_factor/index.html>.

[17] Dexter Kozen and Susan Landau, *Polynomial decomposition
algorithms*, Journal of Symbolic Computation **7**, 445–456 (1989),
author-hosted paper:
<https://www.cs.cornell.edu/kozen/Papers/poly.pdf>.

[18] Christian Gorzel / Singular project, `decomp.lib`, functional
polynomial decomposition, and the fixed outer-degree routine:
<https://www.singular.uni-kl.de/web/lib/decomp.html>,
<https://www.singular.uni-kl.de/web/cmd/decomp/decompunivmonic-23332686.html>.
