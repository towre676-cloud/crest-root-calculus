# Quotients of cyclic covers from one equation

The coefficient-only solver takes a rational polynomial $P$ and an integer $K$ defining $y^K=P(x)$. It removes $K$th powers with their exact rational normalization and requires the remaining degree to be coprime to $K$. The source can have mixed branch multiplicities. It finds the rational maps of degree at least two to unknown squarefree models $z^k=G(u)$ of positive genus and $\gcd(k,\deg G)=1$, with source infinity the only point above target infinity and with the specified cyclic actions intertwined. It searches all possible divisor exponents $k\mid K$ and all surjective characters unless filters are supplied.

Writing $G(U)=BV^k$ and $U'=dV^{k-1}W$ gives a parameter count of $\min(\deg V,\deg W)$. With a squarefree working branch this is also the elimination variable count. Mixed branches can use one extra variable to enforce target squarefreeness. A constant $W$ gives a direct coefficient calculation without Gröbner elimination. The report lists each admissible triple and its geometric support count, then gives the rational identities and original-coordinate formulas. [Theorem and scope](proof.md)

In some degrees, equivariance follows from the source equation. If $kD$ is less than the smallest pole order of any noninvariant function regular off source infinity, any degree-$D$ map in the stated model and infinity class has invariant horizontal coordinate. Its vertical coordinate then transforms by a character. The search covers every such map in that range. For squarefree source branches the pole threshold is the degree of the source radicand. [Pole-order proof](proof.md#when-equivariance-follows-from-the-pole-order)

## Reproduce the results

Use Python 3.11 or newer and SymPy 1.14.

~~~sh
python -B research/cyclic-cover-reconstruction/cyclic.py research/cyclic-cover-reconstruction/picard.json
python -B research/replay.py --quotients-only
~~~

The input coefficient list is in ascending powers. Use **--target-exponent**, **--character** and **--max-degree** to restrict the search. The degree cap is the total source-to-target degree. Save JSON output to a file and pass **--verify FILE** with the original input to replay all certified strata. The returned coordinates are $u=U(x)$ and $z=y^{\mathtt{vertical\_input\_power}}V(x)/D_{\rm input}(x)$, where the chart supplies $D_{\rm input}$ as its vertical denominator. Strata with an empty map list can still have geometric points if their coefficients have no rational solution.

| Source or exact control | Result |
|---|---|
| Smooth Picard quartic | Exactly two normalized maps across both characters, one of degree 2 and one of degree 5 |
| Same quartic with one coefficient changed | Both character searches have empty geometric support |
| Genus-six cubic-cover fixture | Degree 8 recovered with one parameter in place of three, and degrees 5 and 11 excluded in the weight-one chart |
| Saturated recovery | The zero-parameter route recovers exact maps for exponents 3, 4, 5, 7, 9 and 11 |
| Mixed branch multiplicities | Exact maps at exponents 3, 4, 5 and 7, plus a control for repeated target factors |
| Character and exponent changes | Both Picard maps recovered from the other vertical character, and an exponent-six source yields maps of total degrees 2 and 3 through cubic and quadratic subfields |
| Automatic equivariance | A genus-seven cubic cover has one degree-two quotient in the stated singleton-infinity target class at exponent 3. The pole bound covers the entire degree-two search. |

The verifier uses independent leading-coefficient elimination on the Picard positive and negative cases and a mixed exponent-four case. It also substitutes the exported maps in the original function-field equations, checks coverage and rejects a missing stratum. [Independent replay](verify.py)

## Measured against the same native elimination kernel

The frozen study gives each method the same source coefficients, exponent and test stratum. It includes 20 generated strata across six exponents, with positive saturated, nonsaturated and mixed controls and one-coefficient perturbations. The full study preserves 120 fresh-process observations with a four-second limit, two repetitions per input-method pair, and matching exact output checks. All routes use Singular 4.4.1 for the final nonempty multivariate incidence. The derivative route simplifies its preliminary derivative equations in SymPy and decides some strata directly. Solver CPU includes constructing the incidence, eliminating it, reconstructing the result and checking the maps.

| Route | Correct in both runs | Median solver CPU ratio over derivative route | Median fresh-process wall ratio |
|---|---:|---:|---:|
| Complementary derivative reduction | 20/20 | 1.00 | 1.00 |
| Formal-root elimination of $U,G$, without derivative reduction | 10/20 | 1.84 | 1.11 |
| Eliminate the original identity with all of $U,V,G$ unknown | 10/20 | 1.39 | 1.16 |

Each baseline ratio uses its ten jointly completed inputs. The derivative route won eight of ten paired solver-CPU comparisons against the formal-root baseline and seven against the raw baseline. All 80 completed observations agreed with the exact reference. Each baseline reached the time limit on four saturated fixtures, three nonsaturated fixtures and their three perturbations. This is a fixed-stratum benchmark on one generated panel and one host. [Frozen inputs](study/inputs.json) · [Provenance](study/provenance.json) · [Protocol](study/protocol.json) · [All measurements and executable fingerprints](study/measurements.json)

The stored benchmark audit runs without Singular. To rerun timing, use a fresh study directory with the frozen inputs and protocol and run **study.py run --directory** with the directory path. The timing runner requires Linux and Singular. Run **study.py audit** to recompute the mathematical answers and check the archived observations.

## Prior art and open problems

Classical polynomial decomposition covers the constant multiplier. Kozen, Landau and Zippel give a broader resultant framework for algebraic functions. Mendoza computes the Kummer semigroup that underlies the pole-order argument, and effective de Franchis–Severi already gives target enumeration in genus at least two. The complementary derivative recurrence and its explicit cyclic reconstruction scope are the precise candidate contribution. Priority remains open. [Detailed comparison and proof](proof.md#scope-and-antecedents)

The [rational-coordinate solver](branch-poles.md) also recovers quotients with several points over infinity when all finite horizontal poles lie on the source branch locus. It allows higher and mixed pole orders. General infinity fibres, target radicands with mixed multiplicities, and maps that break the cyclic symmetry beyond the automatic range remain open. The benchmark does not measure full unknown-degree discovery or imply a general speed advantage.
