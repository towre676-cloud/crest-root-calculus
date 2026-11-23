# Elliptic maps, projective graphs and the cubic tower

**23 November 2026.** The degree-four search recovers a map from a smooth quartic equation and proves that it stays primitive over an algebraic closure. One example has nineteen geometric degree-four classes. Fifteen factor through double covers, and four are primitive. In degree three the same source has four classes, all defined over $\mathbb Q(\sqrt{-2})$. On an infinite family, all primitive maps in every degree can be counted, and a single curve admits primitive maps of arbitrarily large degree. [Degree-four construction](elliptic-quotients/primitive-degree-four/README.md) · [The nineteen classes](elliptic-quotients/a3-spectrum/README.md) · [All-degree family theorem](elliptic-quotients/primitive-degree-spectrum/README.md)

## Degree four on an arbitrary smooth plane quartic

Three elimination charts now describe every primitive degree-four elliptic map over the algebraic closure, up to translation and isomorphism of the target. A second characterization uses four scalars and 57 Taylor coefficients at an unramified point. On the worked example the reconstruction returns explicit formulas for all four primitive classes, two over the rationals and two over $\mathbb Q(\sqrt{-2})$. The general incidence equations can be generated from the quartic alone. On a fresh positive input, the unrestricted positive chart reached its 90-second time budget. A two-variable elimination and a reduction certificate recovered its entire geometric solution in 7.037 seconds on the benchmark machine. [Equation-only elimination](elliptic-quotients/degree-four-incidence/README.md) · [Local reconstruction and exact formulas](elliptic-quotients/degree-four-jets/README.md)

Two fresh quartics now have complete coefficient-only answers. On the first there are no geometric involutions and the sole primitive elliptic quotient has degree four. An infinite family shares this geometric spectrum. In its smaller arithmetic progression the rational degrees are precisely $4n^2$ and $12n^2$. On the second quartic the solver found two primitive degree-four maps and a double cover, even after mixing all three plane coordinates. Optimal degrees 2, 4 and 4 persist on a second infinite family. [Positive elimination benchmark and first family](elliptic-quotients/flex-fiber-benchmark/README.md) · [Modular reconstruction and second family](elliptic-quotients/cartier-degree-four/README.md)

A separate computation rules out elliptic maps of **every degree** for all integral quartic equations in one explicit residue class modulo 35. Each of its reductions at 5 and 7 has an elliptic factor. At 7 the associated optimal map has degree 41. At 5 the prime-to-5 part of a lifted optimal degree divides 6. These restrictions cannot both hold for an elliptic quotient. [Family, proof and verifier](elliptic-quotients/polarized-reduction/README.md)

## The twelve blind quartics

| Inputs | Possible degrees over an algebraic closure | Possible degrees over the rationals |
|---|---|---|
| Ten Sutherland rows | None | None |
| Row 36518 | Every positive even integer | $2n^2$ or $6n^2$ |
| Row 23640 | Every positive even integer | $2n^2$ or $4n^2$ |

Here $n$ is any positive integer and the elliptic target may vary. Thus degree six occurs rationally on row 36518 and only after field extension on row 23640. A separate elimination computes every geometric degree-two quotient of an arbitrary smooth rational plane quartic. It recovers the ten zero cases, both rational positives, all fifteen involutions of the Fermat quartic and all twenty-one of the Klein quartic. [Complete degree spectra](elliptic-quotients/blind-degree-spectra/README.md) · [Equation-only involution atlas](elliptic-quotients/involution-atlas/README.md) · [Row 36518](elliptic-quotients/row-36518/README.md) · [Row 23640](elliptic-quotients/row-23640/README.md)

Rows 78269 and 32703 have the same point counts over **every** extension of $\mathbb F_2$. At another good prime, a root-ratio certificate proves that there is no nonzero homomorphism between their Jacobians over $\overline{\mathbb Q}$. This is a concrete obstruction to inferring a global correspondence from a complete local zeta function. [Equations and cross-resultant](elliptic-quotients/blind-degree-spectra/README.md#an-exact-local-coincidence-with-no-global-homomorphism)

## Graphs and tower factors

The projective graph program converts homogeneous formulas on an integral projective source into the saturated graph ideal. On the row 23640 family it resolves two base points, checks both projection degrees and computes an intersection of length eight with the graph of the negated map. [Projective graphs](projective-graphs/README.md)

The genus-13 tower stage has five elliptic factors and two absolutely simple fourfolds, all pairwise geometrically nonisogenous. Its complete zeta numerator at 11 is recorded. On the genus-393,661 stage, 21 elliptic curves and 90 absolutely simple fourfolds are individually certified and pairwise geometrically nonisogenous. Their dimension is 381, corresponding to 762 of the 787,322 degrees of the zeta numerator at 41. The subset and sign projectors apply at every finite stage and split the integral Tate modules away from 2 and 3. [Cubic-tower theorem and exact polynomials](cubic-tower-factors/README.md)

The stage-ten Jacobian contains **at least 2,047 distinct simple isogeny classes over the rationals**. The traces of its 2,047 signed subset factors, including the base, have full column rank modulo 127 across 2,055 good primes. This proves that the 1,936 higher-support factors each contain a constituent that can be chosen distinct from those assigned to all the other factors. Each signed factor also has an explicit geometric model up to isogeny for its square, using an unramified cyclic cover of degree three or six. Subsequent geometric certificates below determine many of the higher factors individually. The computations also give the full genus-55 zeta numerator at 11 modulo 3 and the 240 signed triple-support polynomials at 41 modulo 3. [Rank certificate, congruences and replay](cubic-tower-higher-supports/README.md)

The genus-55 stage now has a complete geometric decomposition. At 11 its entire degree-110 zeta numerator is known modulo 33, with its first eight nonconstant coefficients known over the integers. Among the new exact counts are 19,487,393 points over $\mathbb F_{11^7}$ and 214,404,341 over $\mathbb F_{11^8}$. Four central integer coefficients per signed top factor remain undetermined. [Differential matrices, exact counts and coefficients](cubic-tower-cartier/README.md)

The geometric decomposition is complete through stage seven. Its genus-10,207 Jacobian has 255 pairwise nonisogenous geometric simple factors, including two of dimension 448. At stage ten, **1,293 of the 2,047 displayed factors** are individually certified geometrically simple and pairwise nonisogenous. They account for 60,989 of the 393,661 dimensions. These include all 240 three-index, 420 four-index and 504 five-index factors, sixteen six-index factors and two seven-index factors. The geometric decomposition within the other 754 displayed factors remains open. All 420 six-index factors are pairwise geometrically nonisogenous as whole abelian varieties. A single reduction at 193 also distinguishes constituents inside partly classified factors. This gives **at least 1,599 distinct geometric simple classes** at stage ten, of combined dimension at least 102,691. [Seventh-stage decomposition](cubic-tower-stage-seven/README.md) · [Geometric separation and full recount](cubic-tower-geometric-separation/README.md) · [Five- and six-index computations](cubic-tower-size-five/README.md)

## Run the research

Install Python 3.11 or later and extract the packet, or clone the repository. From the directory containing `research`, run

```sh
python -m pip install numpy "sympy>=1.12,<2"
python research/replay.py
```

On Windows, `py` can replace `python`. The replay includes the original twelve-row solver with its historical limits, both positive and negative controls, new coefficient-only searches and exact agreement with the saved receipts. The default checks every dimension and separation certificate and rebuilds selected higher-support matrices. Add `--recount --workers 4` to rebuild every archived higher-support observation. The optional C++17 count over the field of 214,358,881 elements uses about 2.2 GB of table storage. The individual programs can also run separately.

| To use a research program | Command from the repository root |
|---|---|
| List all geometric double quotients and their equations | `python research/elliptic-quotients/involution-atlas/classify.py research/elliptic-quotients/row-36518/quartic.json` |
| Discover a primitive degree-four map | `python research/elliptic-quotients/primitive-degree-four/search.py research/elliptic-quotients/primitive-degree-four/sheared.json --all-rational-centers` |
| Generate the general primitive degree-four equations on one chart | `python research/elliptic-quotients/degree-four-incidence/incidence.py research/elliptic-quotients/primitive-degree-four/ciani.json 2 --singular chart2.sing` |
| Exhaust the elliptic maps on the new involution-free positive | `python research/elliptic-quotients/flex-fiber-benchmark/branch_fibres.py research/elliptic-quotients/flex-fiber-benchmark/quartic.json` |
| Discover both primitive quartic maps in mixed coordinates | `python research/elliptic-quotients/cartier-degree-four/cartier_search.py research/elliptic-quotients/cartier-degree-four/mixed.json` |
| Rebuild the higher-support certificate from its local matrices | `python research/cubic-tower-geometric-separation/verify.py --recount --workers 4` |
| Recover a map from its normalized local parameters | `python research/elliptic-quotients/degree-four-jets/reconstruct.py research/elliptic-quotients/degree-four-jets/example.json` |
| Test two good reductions for an obstruction in every degree | `python research/elliptic-quotients/polarized-reduction/screen.py research/elliptic-quotients/polarized-reduction/quartic.json --primes 5 7` |
| Rebuild all 504 dimension-80 certificates | `python research/cubic-tower-size-five/verify.py --recount --workers 4` |
| Recompute the higher-dimensional witnesses for stage seven | `python research/cubic-tower-stage-seven/verify.py --recount --workers 4` |
| Verify the stage-ten rational isogeny-class lower bound | `python research/cubic-tower-higher-supports/prime_rank.py` |
| Count the nineteen geometric degree-four classes | `python research/elliptic-quotients/a3-spectrum/map_classes.py 4 --directions` |
| Find the primitive degree-26 maps at parameter 4 | `python research/elliptic-quotients/primitive-degree-spectrum/degrees.py 4 26 --directions` |
| Compute the saturated graph of a map | `python research/projective-graphs/graph.py research/projective-graphs/row-23640.json` |

The general incidence and jet criteria apply to arbitrary smooth quartics. The two new positive algorithms each have a case-specific search and a geometric test of completeness. The full raw characteristic-zero elimination on an arbitrary positive input, the largest stage-ten factors and the integral higher zeta polynomials remain open. [The 57-item research register](backlog-2026-11-23.md)
