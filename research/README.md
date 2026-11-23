# Elliptic maps, projective graphs and the cubic tower

**23 November 2026.** The degree-four search recovers a map from a smooth quartic equation and proves that it stays primitive over an algebraic closure. One example has nineteen geometric degree-four classes. Fifteen factor through double covers, and four are primitive. In degree three the same source has four classes, all defined over $\mathbb Q(\sqrt{-2})$. On an infinite family, all primitive maps in every degree can be counted, and a single curve admits primitive maps of arbitrarily large degree. [Degree-four construction](elliptic-quotients/primitive-degree-four/README.md) · [The nineteen classes](elliptic-quotients/a3-spectrum/README.md) · [All-degree family theorem](elliptic-quotients/primitive-degree-spectrum/README.md)

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

## Run the research

Install Python 3.11 or later and extract the packet, or clone the repository. From the directory containing `research`, run

```sh
python -m pip install numpy "sympy>=1.12,<2"
python research/replay.py
```

On Windows, `py` can replace `python`. The replay includes the original twelve-row solver with its historical limits, both positive and negative controls, new coefficient-only searches and exact agreement with the saved receipts. The stage-ten test builds a table over a field of 2,825,761 elements. The individual programs can also run separately.

| To use a research program | Command from the repository root |
|---|---|
| List all geometric double quotients and their equations | `python research/elliptic-quotients/involution-atlas/classify.py research/elliptic-quotients/row-36518/quartic.json` |
| Discover a primitive degree-four map | `python research/elliptic-quotients/primitive-degree-four/search.py research/elliptic-quotients/primitive-degree-four/sheared.json --all-rational-centers` |
| Count the nineteen geometric degree-four classes | `python research/elliptic-quotients/a3-spectrum/map_classes.py 4 --directions` |
| Find the primitive degree-26 maps at parameter 4 | `python research/elliptic-quotients/primitive-degree-spectrum/degrees.py 4 26 --directions` |
| Compute the saturated graph of a map | `python research/projective-graphs/graph.py research/projective-graphs/row-23640.json` |

The first two programs take quartic monomial lists. Degree-four discovery requires the pair of commuting rational involutions and elliptic-quotient identifications described in the proof. Exhaustive primitive degree-four discovery on arbitrary plane quartics, the higher-support simple factors at stage ten, and general construction from local projector data remain open. [The 57-item research register](backlog-2026-11-23.md)
