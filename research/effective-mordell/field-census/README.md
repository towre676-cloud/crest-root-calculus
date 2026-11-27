# Effective Mordell on the curve y² = x⁶ − 1

The main PDF, crest-root-height-and-field-census.pdf, is the eight-page proof note. It gives explicit height bounds, an infinite prime family, exact point classifications, and a lower growth requirement on a universal quadratic-field bound. The ZIP contains the mathematical derivations, number-field coordinates, code and replay inputs. All point counts include the two points at infinity on the smooth projective model. A **partial** result records verified points and leaves completeness open.

## Results

| Input or hypothesis | Conclusion | Proof or executable |
|---|---|---|
| All signed squarefree quadratic radicands through 200 | 243 complete lists. Exactly 232 have four points, one has six and ten have eight | quadratic-results.md, verify_panel.py, final_three.m |
| All signed squarefree radicands through 1000 | 1,179 of 1,215 complete. The 36 outstanding fields occur in 18 reciprocal pairs. All 936 additional complete fields have four points | extended-quadratic.json, verify_extended_quadratic.py |
| An explicitly recognizable set of primes of natural density 927/1024 | Both signs of the prime yield exactly four points | quadratic-results.md, prime_density.py |
| Elliptic rank at most one over K(i), or at most two over K(ζ₁₂) | H(x) < 4.279, respectively H(x) < 4.965, for every finite abscissa of the fixed curve | height-and-census.tex, height_constants.py, cm_height.py |
| Q(ζ₂₄) | 46 points, Jacobian rank four, two independent exhaustions | zeta24.py, OCTIC-FIELDS.md |
| Octic fields Q(ζ₁₂, √d) for squarefree d from 2 through 200 prime to 3 | 49 complete of 91. Forty-eight have ten points; d = 2 has 46. The d = 7 and 11 lower bounds are 130 and 58 | octic_verify.py, octic-panel.json |
| Higher genus and other field classes | Smooth genus-55 example with a complete list over three fields, non-Galois towers of unbounded degree, class-group descent with 123 complete of 539 tested models | COVER-TRANSFER.md, tower_theorems.md, CLASSGROUP-DESCENT.md |
| Reconstructed degree-four maps on the positive sextic y² = x⁶ + 1 | Two of the maps independently give a finite height search under the stated rank hypothesis. The direct reciprocal pair gives the sharper bound | quartic_bridge.py, quartic-height.json |
| Real quadratic fields varying on the fixed negative sextic | The ratio of H(x) to the square root of the field discriminant is unbounded, even for rational abscissas. A logarithmic strengthening and explicit high-rank fields are proved | FIELD-HEIGHT-OBSTRUCTION.md, discriminant_obstruction.py, rank_growth.py |

## Run a field

Supply an irreducible monic polynomial as ascending integer coefficients, with the generator denoted by w in the answer. The field with polynomial w represents the rationals.

~~~sh
python solve_field.py '[-182,0,1]'
python solve_field.py '[1,0,0,0,-1,0,0,0,1]'
python cm_height.py '[0,1]' --enumerate
python cover_transfer.py plane '[-17,0,1]' 2 1/2 1/2
~~~

Each result lists affine coordinates. The two rational infinities are recorded separately. The height enumerator uses a finite coefficient box. In large field degree its raw enumeration can be extremely expensive. The practical cyclotomic classifier instead saturates an elliptic lattice and checks its finite set of possible fibres. The rank-dependent criterion uses rigorous upper bounds. A field can remain partial because its rank is too large or because the chosen test is inconclusive.

## Replay

The tested environment used Python 3.12, cypari2 2.2.4 with PARI 2.17.2, SymPy 1.14.0 and python-flint 0.9.0. Install these three Python packages in a PARI-capable environment, then run from this directory with assertions enabled.

~~~sh
python -m pip install -r requirements.txt
python verify_all.py
python verify_all.py --magma /path/to/magma
~~~

The second command runs the Python proofs and verifies archived outputs. The third also runs the Magma assertions for the elliptic Chabauty steps. The combined three-twist computation was executed in Magma V2.29-10; chabauty-run.json records the independent run. Python checks the hypotheses coming from the local sieve and the full lifts from the rational projections. Magma checks the ranks, subgroup index condition and projection sets. Periods, densities and printed decimal height bounds have exact identities, integer counts or interval comparisons behind them.

The classifier uses class-group and unit data after PARI certification succeeds. No uniform runtime or termination guarantee for rank or class-group computation on arbitrary inputs is claimed. The later [biquadratic proof](../biquadratic-38/crest-root-biquadratic-38-proof.pdf) establishes that all 38 points over Q(i, √7) are accounted for. The classifier includes that proof. The field Q(ζ₁₂, √7) of degree eight is a different input. Method-specific tables for Gaussian, octic and general Selmer runs can overlap. Their counts must not be summed.

## Mathematical context

Bogomolov–Tschinkel and Poonen give a reduction from general hyperelliptic and superelliptic input to the fixed negative sextic over varying number fields. The rank and field hypotheses here still exclude inputs. The public statement of ProofAtlas #87 also covers arbitrary smooth projective curves of genus at least two. No claim of priority is made for these classifications or their methods.

Poonen, [Unramified covers of Galois covers of low genus curves](https://math.mit.edu/~poonen/papers/etale.pdf). Garcia-Fritz and Pasten, [Effective Mordell for curves with enough automorphisms](https://arxiv.org/abs/2503.10443). [ProofAtlas open problems](https://www.proofatlas.ai/open-problems/).
