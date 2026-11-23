# Configuration curves and applied certificates

**Can a few checks settle a question about every possible case?** Turn a crank through 720 positions and you can still miss the place where it jams. Sometimes geometry tells you exactly where to look. Sometimes the best available test only gives you a safe answer on one side of a boundary. These companion Parts work out the difference, using the three released products.

| Read | Question | What you can reproduce |
| --- | --- | --- |
| [Part I · Configuration curves of planar mechanisms](mechanisms/README.md) | Will the input crank get all the way round? Can alternating two reflections ever return a wing to its starting pose? | A missed interval only about 0.26° wide, exact turning tests, and proofs of infinite folding order for the nominal and printable wings. |
| [Part I continued · Six-bar mechanisms](sixbar/README.md) | When do two coupled four-bars make a genus-three curve with a map of degree three? | Branch overlaps, a hyperelliptic example, a plane quartic with only even-degree elliptic maps, and a symmetric construction with a degree-three geometric quotient. |
| [Part II · Certificates for spectra and electricity networks](applied.md) | Which grain proportions fit a scan? What can a voltage bound promise a street of solar homes? | Exact wavelength exclusions, a weak-absorption fibre and a rational certificate for 18 modeled connections. |

## Worst-case margins

**Definition for these applications.** Let $\theta$ describe a design and $\omega$ a state in a specified set $\Omega$. Give the design a margin $M(\theta,\omega)$ with nonnegative values acceptable. Write

$$
M^*(\theta)=\inf_{\omega\in\Omega}M(\theta,\omega),\qquad
\mathcal W=\{\theta:M^*(\theta)=0\}.
$$

The application calls $M^*$ the crest and $\mathcal W$ its root wall. For a declared choice of units and coordinates, $\tau_\infty=\inf_{\theta'\in\mathcal W}\|\theta-\theta'\|_\infty$ measures distance to that wall. With a continuous margin and a separating boundary, this measures the clearance to first contact from the strict interior. Raw millimetres and changes of mineral proportion live in different parameter spaces. The historical monograph also uses roots and crests for arithmetic data and growth patterns. [Original terminology and examples](https://github.com/towre676-cloud/crest-root-calculus#where-do-crest-and-root-come-from)

For any finite sample $S\subseteq\Omega$, $M^*(\theta)\le\min_{\omega\in S}M(\theta,\omega)$. A negative sampled margin can refute a universal claim. To certify a pass, give a reason that the chosen tests cover every state: the linkage distance attains its extrema, the full affine spectral table reduces to extreme intervals, or a monotone voltage bound dominates every operating point. Sampling **can** settle the question when such a reason is available.

| Product | Margin and states | Scope of the verdict |
| --- | --- | --- |
| Linkwright | Triangle or transmission clearance at crank angle $\omega$ | Two endpoint signs exactly decide strict passage through the entire turn of the ideal rigid linkage. Algebraic invariants describe the curve of poses. |
| GeoAssay | Tolerance minus absolute residual at each measured wavelength, for a fixed candidate model | Intersecting affine-band intervals exactly decides feasibility for the table. In the nonlinear grain-size model, loss of derivative rank concerns parameter identifiability. That critical locus and the tolerance boundary answer different questions. |
| FlexExport | Squared-voltage headroom at every bus and every allowed operating state | The upper corner decides the chosen monotone voltage screen. Passing certifies upper voltage in the stated AC model. A rejected offer can still work. |

The analogy with an algebraic certificate is about the **direction of implication**. For instance, reducing generators of $I$ to zero against a candidate basis $G$ can show $I\subseteq\langle G\rangle$. Even $\langle x^2\rangle\subseteq\langle x\rangle$ passes that check. To establish equality one must also put every element of $G$ back in $I$, for example by explicit cofactor identities. To assert a Gröbner basis, its Gröbner property must also be established. Each application below states its own corresponding proof obligation.

## Replay the research

From this folder run the following command with Python 3.11 or later. It uses the standard library, works offline and compares independently recomputed results to all four stored receipts.

```sh
python verify.py
```

The six-bar checker counts points on the normalized covers and tests the quotient chart. The other checkers exercise the two wing designs, both measured powder witnesses, and the feeder's exact primal and dual bounds. Inside the combined ZIP the same command also rechecks the selected wavelength rows against the archived raw USGS files. The local notes distinguish general proofs, finite verification, conditional calculations and the obstruction to descending the exhibited cyclic quotient to the real physical form. The source for the three products is the [pinned expanded release](https://raw.githubusercontent.com/towre676-cloud/crest-root-calculus/4d311c055d98e8155d6d12d73588f747de034545/downloads/crest_root_three_builds_expanded.zip).
