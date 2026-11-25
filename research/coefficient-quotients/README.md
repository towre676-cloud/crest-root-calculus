# Recover a quotient from the source coefficients

Give the solver a single expanded rational polynomial $P$. It returns a squarefree normalization $P=A^2F$, then searches for positive-genus quotient curves of $y^2=F(x)$. The source coefficients supply the entire discovery input.

- **Constant vertical multiplier.** For every possible inner degree, find the normalized $F=G(U)$ or certify its exclusion. The solver accepts both odd and even source degrees.
- **Nonconstant vertical multiplier.** For odd-degree $F$, find every rational map of degree at least two to an unknown odd-degree hyperelliptic or elliptic model $v^2=G(u)$ whose only point over target infinity is source infinity. The certificate uses $G(U)=FV^2$ for the map $(u,v)=(U,yV)$. The output uses monic $U,V$, $U(0)=0$, and counts normalized maps up to pointed target isomorphism.
- **Derivative reduction.** For $d=\deg U$ and $e=\deg V$, write $U'=dVW$. The leading equations recover $V$ from $W$, so the incidence uses $\min(e,d-1-e)$ variables. A saturated stratum $e=d-1$ needs no multivariate elimination. [Proof and worked degree-five quotient](proof.md)

The arbitrary local-to-global problem remains open. This theorem uses the input equation and a bounded coefficient incidence. The [equal-count obstruction](../local-pattern-obstruction/README.md) explains why a general rational reconstruction cannot depend on point counts alone.

## Run and verify

From the repository root with Python 3.11+ and SymPy 1.14:

~~~sh
python -B research/coefficient-quotients/discover.py research/coefficient-quotients/example.json --nonconstant
python -B research/coefficient-quotients/verify.py
~~~

The input is a JSON list of integer or rational-string coefficients in ascending order. The example is the degree-nine branch polynomial from the proof. Its genus-four curve has a degree-five map to an elliptic curve, with $V$ cubic and $W$ linear. The output includes the second admissible stratum and its exact exclusion. For an inflated input $y_P^2=P$, the original vertical coordinate in the map is $y_P V/A$.

Each variable-multiplier stratum includes its geometric support count and its complete rational list. The verifier replays the rational coefficient algebra, checks each identity in the original equations, and checks four cubic strata by separate leading-coefficient elimination. It independently constructs all six Jacobi examples from full transfer products. SymPy is the only additional dependency for replaying the stored mathematical results.

## Frozen comparison

The held-out panel includes 36 constant-multiplier inputs and 20 nonconstant-multiplier inputs. Algebraic covers, rational affine changes, square-factor inflation, several admissible degrees, one-coefficient perturbations, a stock SymPy regression, and six Jacobi spring-model equations are distinguished in the provenance. The generator tests smoothness. The exact replay checks each perturbation's answer. The six Jacobi inputs have explicit stiffness and nonnegative pinning parameters, including a declared nonminimal cell. These are exact spring models.

Each constant method searches every admissible divisor and receives only the public coefficient list. The variable-multiplier comparison also supplies the common test stratum $(d,m,e)$ to its two solvers to isolate the effect of the derivative. The independent replay runs the full degree list for each variable input. The benchmark retains a separate Crest run with the modular sieve disabled, corrected dense decomposition, native Singular with its fixed-degree Newton routine, PARI/GP with exact series and division, and stock SymPy as a version-specific correctness control. The dense and Singular adapters share the host's exact squarefree preprocessing. GP uses its own monic-gcd normalization.

[Held-out inputs](study/inputs.json) · [Separate provenance](study/provenance.json) · [Freeze](study/freeze.json) · [Exact constant references](study/constant-reference.json) · [All observations](study/measurements.json) · [Verified summary](study/summary.json)

New timings require Linux, Singular, and PARI/GP:

~~~sh
python -B research/coefficient-quotients/benchmark.py --directory /path/to/copied-study
~~~

The runner pins one logical CPU and uses three fresh processes per input and method in seeded order. It measures both CPU and wall time, as well as normalization, search and verification phases where available. The GP millisecond native phases and its wrapper/launch-inclusive CPU are recorded separately. The archived environment uses Singular 4.4.1 and PARI/GP 2.17.2. CREST_QUOTIENT_SINGULAR, CREST_DECOMP_LIB and CREST_QUOTIENT_GP allow explicit executable and library paths. Only runs with three correct completions per method enter paired speed ratios.

A preliminary timing attempt was stopped after discovering that the Singular adapter could misparse unparenthesized rational coefficients. Its original freeze is retained in [preflight-freeze.json](study/preflight-freeze.json). The rendering was corrected, all four supported constant routes agreed on the separate training panel, and a fresh holdout was drawn after freezing the revised code. The preliminary attempt produced no completed measurement file.

### Observed results

The independent replay checks **123 constant-degree decisions across the full panel**, **34 nonconstant degree strata**, and **four raw cubic eliminations**. On the 36 constant-study inputs the exact methods recovered 29 normalized targets. Three of six Jacobi cases have a positive-genus degree-two quotient, including a repeated cell after normalization. In both commuting-composition controls there are three admissible quotient degrees. All ten constant one-coefficient perturbations have an empty quotient list in this scope.

| Route | Correct constant inputs | Solver-phase CPU wins against Crest | Median route CPU / Crest CPU |
|---|---:|---:|---:|
| Crest without modular sieve | 36/36 | 16/36 | 1.017 |
| Corrected classical dense | 36/36 | 32/36 | 0.773 |
| Singular adapter | 36/36 | 0/36 | 4.096 |
| PARI/GP adapter | 36/36 | 3/36 | 4.462 |
| Stock SymPy 1.14.0 | 23/36 | 19/23 | 0.882 |

The solver-phase ratios include normalization, quotient discovery and exact output checks. The native adapter figures also include subprocess launch, conversion and output parsing. The wrapper-inclusive comparisons therefore establish full-task measurements on this host, and a persistent native-kernel comparison would answer a different question. The corrected dense algorithm wins on 32 of 36 solver-phase measurements and is the strongest observed constant-route comparator. Its median advantage is about 1.29×. Fresh-process wall time for the two Python approaches is nearly equal. Stock SymPy misses at least one map on 13 inputs, each in all three repetitions.

The derivative solver completes all **20/20** supplied nonconstant strata within the four-second fresh-process limit. The ablation, which uses the same formal-root and finite-algebra engines while omitting the derivative restriction and complementary parametrization, completes **12/20**. All four quintic cover fixtures and four of their perturbations exhaust its budget in every repetition. On the 12 jointly completed inputs the derivative solver wins 10 and has a **1.80× median solver-phase CPU advantage**. The median fresh-process wall ratio is **1.06×**. The complete nonconstant replay additionally certifies every other admissible stratum for these sources and recovers ten rational quotient maps across the 20 inputs.

Across **768** attempts there were **705** correct completions, **39** stock-SymPy wrong answers and **24** ablation timeouts. The derivative comparison is an ablation of the implementation on a designated nonconstant class. It shows where the new coefficient reduction matters and leaves a native external comparison for that class open. [Summary by cohort and metric](study/summary.json)
