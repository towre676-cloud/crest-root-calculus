# Exact solver comparisons on three mathematical tasks

These experiments compare the project’s algorithms with native Singular, SymPy and PARI/GP on independently specified inputs. They cover marked degree-three elliptic pencils, polynomial hyperelliptic maps and complete local zeta numerators. A separately labelled project implementation of a corrected classical dense-decomposition recurrence supplies an additional comparator. Returned answers are checked exactly. The map study compares every method with a frozen, replayed reference and records incorrect outputs separately. The raw receipts retain every incorrect output, timeout, repetition and resource limit. [All measured results](results.md)

The panels contain eighteen marked quartics, eighteen polynomial-map problems and fifteen point-count problems: 312 planned measurements across the three published studies. Input files and executable adapters are hashed before their respective timed studies. Each method runs in a fresh process, serially on the same selected CPU, with two repetitions in a seeded order. The solver CPU timer includes normalization, mathematical preprocessing, exact elimination and result extraction. The separate process wall timer also includes Python startup. Timings refer to one Linux host. Two repetitions describe these deterministic runs and their observed variation. [Quartic protocol](quartic-results-protocol.json) · [Map protocol](map-results-protocol.json) · [Counting protocol](frobenius-results-protocol.json)

The marked-quartic study was interrupted after 55 recorded runs. The continuation preserved every recorded result and resumed the exact remaining seeded schedule on the same machine and software. The original input and adapter hashes were verified before continuation. The three imported quartic research modules were first hashed at resumption and checked again at completion; those dependency hashes were not part of the original protocol. The initial execution requested CPU eight, but its JSON protocol recorded only the available affinity set and omitted a separate selected-worker-CPU field. The resumption and both later protocols explicitly record CPU eight. The map and counting studies record their mathematical dependency hashes before timing and check them afterward. [Resumption record](quartic-results.json) · [Resumption runner](resume_quartic.py)

The CPU measurements use POSIX `getrusage`, summing user and system time for the worker and its reaped children. This includes native Singular and PARI work, with seconds represented at microsecond precision. The PARI adapter includes process launch and decoding; it does not use GP's millisecond timer. Per-input numbers describe this implementation, hardware and resource budget. Two repetitions support a reproducible comparison here and do not establish population-level statistical precision.

## Observed results

All 312 planned measurements are recorded. Across the published studies, 254 runs returned the verified answer, 38 reached the wall limit and 20 returned an incorrect answer. All 20 incorrect answers came from the stock SymPy 1.14.0 decomposition comparator, on the same ten inputs in both repetitions. There were no process errors. [Complete per-input tables](results.md) · [Machine-readable summary](summary.json)

The ratios below divide the comparator's median time by Crest–Root's median time on each input, then take the median across inputs that both methods completed correctly twice. A ratio above one means the comparator took longer. CPU and process-wall ratios use the same paired inputs.

| Task and comparator | Correct inputs completed twice: Crest / comparator | Paired inputs | Median CPU ratio | Median process-wall ratio |
|---|---:|---:|---:|---:|
| Marked quartics / Singular | 18/18 / 16/18 | 16 | 2.946 | 2.672 |
| Marked quartics / SymPy Buchberger | 18/18 / 4/18 | 4 | 4.448 | 2.092 |
| Polynomial maps / Singular | 18/18 / 15/18 | 15 | 2.084 | 1.285 |
| Polynomial maps / stock SymPy decomposition | 18/18 / 8/18 | 8 | 0.862 | 0.963 |
| Polynomial maps / corrected classical dense method | 18/18 / 18/18 | 18 | 0.883 | 0.935 |
| Small-prime zeta cohort / PARI | 12/12 / 12/12 | 12 | 1.844 | 1.081 |
| Larger-field zeta cohort / PARI | 3/3 / 3/3 | 3 | 0.183 | 0.536 |

The marked-quartic method has the strongest observed advantage: all eighteen cases finish correctly, and it uses less CPU than Singular on twelve of the sixteen jointly completed cases. Singular is faster on the other four. For polynomial maps, the corrected classical dense method is faster than the project implementation on fourteen of eighteen inputs. It computes the same complete answer and supplies a useful implementation target. The stock SymPy figures use only its eight correct inputs; its ten false-negative inputs remain explicit correctness failures.

Both counting methods return the same complete local polynomial on all fifteen inputs. The small-prime cohort uses primes five, seven and eleven; interpreter startup largely removes the solver-CPU difference in its process-wall comparison. On the three larger-field inputs, PARI is faster every time: the median direct-counter/PARI ratio is 5.463 on CPU and 1.864 on process wall. These are complete invocations, including PARI launch and decoding. The predeclared cohorts reveal a change in the faster method as the required finite fields grow. Further scaling and native-kernel measurements within a persistent process remain open.

## The same output for each method

### Complete transverse pencils on a marked quartic

The input is a smooth ternary quartic over the rationals together with a rational smooth point. The output is the reduced geometric support of every transverse degree-three elliptic pencil at that mark, with its target equation. The method receives the quartic coefficients and the point. The mark is specified as part of this task.

The six square-completion equations from the original chart are denoted by $N_{ij}$. On the open set $h\ne0$, two of them are equivalent to the osculation equation $q=0$ and the forced conic relation $Db-B=0$. Explicitly,

$$
N_{22}=16h^2(3w^2q),\qquad
N_{12}=16h^2\bigl(Db-B-bu(3w^2q)\bigr).
$$

The generic solvers receive those two lower-degree equations and the other four original equations. Powers of $h$ dividing an individual generator are removed. Adding $1-zh$ enforces the same open set. These changes preserve the solution scheme after inverting $h$ and give the generic solvers the elementary simplifications already justified by the chart identities. [Generator identities](../chart-one-cofactors/README.md)

Crest–Root uses its rank-four norm elimination. Singular uses its documented heuristic `groebner` routine in degree reverse lexicographic order, computes the radical, and converts with FGLM. SymPy uses its default Buchberger algorithm in graded reverse lexicographic order and converts with FGLM. The final exact fibre calculation removes nilpotents and reconstructs the maps. Complete outputs are compared as canonical rational lexicographic ideals in the remaining coordinates. [Singular Gröbner methods](https://www.singular.uni-kl.de/Manual/latest/sing_266.htm) · [Singular radical](https://www.singular.uni-kl.de/Manual/latest/sing_1541.htm) · [SymPy methods and FGLM](https://docs.sympy.org/latest/modules/polys/reference.html#sympy.polys.polytools.groebner)

Six inputs use different square constructions, eight are sampled from marked quartics with coefficient heights one through eight, and four are fixed classical curves. The generator rejects singular curves. No solver result enters that selection rule. The training positive used to configure Singular is absent from the frozen panel. [Input generator](panel.py) · [Frozen quartics](quartic-panel.json)

### Polynomial hyperelliptic maps

The input is a pair of squarefree branch polynomials defining $y^2=F(x)$ and $v^2=G(u)$. This comparison asks for all rational maps of the form

$$
u=U(x),\qquad v=c\,y,\qquad G(U)=c^2F,\qquad c\ne0.
$$

The degree of $U$ is forced by $\deg F/\deg G$. The output includes every rational coefficient list and the exact number of geometric maps. The panel varies target genus from one to four and map degree through one, three, five and seven. The largest branch polynomial has degree sixty-three. Each of six independently seeded constructions supplies a cover, a quadratic twist and a coefficient perturbation. Smoothness is checked before the panel is frozen. [Input generator](map_panel.py) · [Frozen map problems](map-panel.json)

The project extracts the uniquely normalized inner polynomial from the leading coefficients, then solves the remaining affine change of target and vertical scaling. The Singular baseline receives the full coefficient equations for an unknown polynomial $U$ and an unknown scalar $c$, with an inverse for the leading coefficient of $U$. Its result is converted to exact rational maps.

The stock SymPy comparator uses its dense fixed-degree functional-decomposition routines, underlying the public `decompose` function. On ten frozen inputs, SymPy 1.14.0 misses a decomposition that has an exact polynomial identity. The same failure reproduces through public `decompose`; it is not confined to the adapter. Those results are retained as `wrong_answer`, with their actual outputs and runtimes. They do not count as correct completions or enter speed ratios. [Minimal public reproduction and recurrence derivation](sympy-decomposition-regression.md)

A fourth comparator, `corrected-dense-decomposition`, uses the same dense operations with the corrected right-component recurrence. It is an explicitly labelled project-side implementation of the classical method, not unmodified SymPy. All three decomposition methods then use the same residual affine equations and rational-point extraction. Comparing their correct outputs measures the implementation of decomposition within this task. The primary classical reference is Kozen and Landau’s 1989 algorithm. [SymPy implementation](https://github.com/sympy/sympy/blob/1.14.0/sympy/polys/densetools.py) · [Kozen and Landau](https://www.cs.cornell.edu/kozen/Papers/poly.pdf) · [Map theorem and exact reconstruction](../bounded-correspondences/README.md)

The original map attempt stopped at the first observed disagreement. Its [11 recorded rows](map-aborted-results.json), [protocol](map-aborted-protocol.json) and [original adapter](map-worker-before-correction.py) are preserved. An audit also identifies an earlier false negative in that prefix. The final study keeps the same input panel, freezes [exact reference certificates](map-results-reference.json) before timing, adds the corrected comparator under a separate name and reruns every scheduled input-method pair. No input was removed or replaced because of the failure.

### Complete local zeta numerators

The input is a squarefree hyperelliptic equation and a prime of good reduction. The output is the entire local zeta numerator, with its first $g$ extension-field point counts. The project explicitly constructs the finite fields, counts quadratic characters over every required extension and applies Newton identities. PARI/GP computes `hyperellcharpoly`, and the adapter converts its characteristic polynomial to the same reciprocal-coefficient convention. Both methods include validation and all per-input initialization in their solver time.

The twelve original seeded examples cover genera two through four at primes five, seven and eleven. Three have even branch degree, so agreement also checks the full infinity fibre. A second seed adds three smooth inputs of genus two at prime 1009, genus three at prime 97 and genus four at prime 31. Their largest extension fields have 1,018,081, 912,673 and 923,521 elements. The small-prime and larger-field cohorts are fixed before timing and reported separately, as well as in the all-fifteen totals. The direct counter has a two-million-element field-table budget. Every method retains the same ten-second process and four-GiB memory limit. Each repetition starts a fresh process and an empty field cache. [Frozen counting inputs](frobenius-panel.json) · [Counter and adapter](frobenius_benchmark.py) · [PARI/GP documentation](https://pari.math.u-bordeaux.fr/dochtml/html-stable/Elliptic_curves.html#hyperellcharpoly)

The [original twelve-input panel](frobenius-preliminary-panel.json), [48 measured runs](frobenius-preliminary-results.json), [original protocol](frobenius-preliminary-protocol.json) and [original adapter](frobenius-preliminary-worker.py) are preserved as the preliminary study. The final panel keeps all twelve original inputs unchanged and appends the three declared larger-field inputs. Every final input was timed afresh; the published counting table uses only the final sixty measurements.

## Interpretation and reproduction

The completion tables and paired ratios answer different practical questions. Correct completion records whether an input returned the verified answer within the specified limit. Wrong answers and timeouts have separate columns. Paired CPU and process-wall ratios compare two finite runtimes on the same input. The latter includes interpreter startup and shows the practical cost of a fresh invocation. The report retains timed-out and incorrect-answer inputs in the outcome totals and displays them explicitly in the full per-input tables. Comparing medians from different sets of solved inputs can create selection bias, so the reported speed ratios use jointly completed correct inputs only.

This study addresses three defined tasks. The original M02 register also asks about other algorithms and problem families, including general-purpose symbolic systems, determinant methods and sparse constraint solvers. Those comparisons remain part of the wider program. The source code makes the present experiments reusable on additional inputs and machines. [Research register](../backlog-2026-11-23.md)

The runners use POSIX process groups, memory limits and optional CPU affinity. Use Python 3.11 or later; the published runs used Python 3.12.14. Install the pinned Python dependencies from `research/solver-comparison/requirements.txt`, Singular with its standard libraries and dynamic polynomial modules, and PARI/GP. The conventional decomposition adapter uses SymPy 1.14.0 private dense-polynomial routines, so keep that version for exact reproduction. Set `SINGULAR` and `PARI_GP` to executable paths when they are absent from `PATH`. From the repository root, run

```sh
python research/solver-comparison/verify.py
python research/solver-comparison/verify.py --recompute
python research/solver-comparison/report.py
python research/solver-comparison/verify_sympy_regression.py
```

The default verifier checks the regression’s exact identities without requiring newer SymPy versions to retain the old failure. The separate regression command checks the archived stock 1.14.0 implementation and behavior; use its `--observe` flag to inspect another installed version.

To time fresh studies while preserving the published measurements, supply fresh output paths.

```sh
python research/solver-comparison/benchmark.py --output quartic-new.json --cpu 0
python research/solver-comparison/map_benchmark.py --output maps-new.json --cpu 0
python research/solver-comparison/frobenius_benchmark.py --output counts-new.json --cpu 0
```

The quartic process limit is thirty seconds. The map and counting limits are ten seconds. Every run has a four-GiB address-space limit. The published studies used CPU eight, Python 3.12.14, SymPy 1.14.0, Singular 4.4.1 and PARI/GP 2.17.2. The protocols retain the exact system description and the raw receipts preserve solver CPU and process wall measurements. A fresh run should retain its own protocol and report its own observed costs.
