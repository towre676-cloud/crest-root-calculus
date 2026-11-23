# Exact periodic spectra: measured results

PARI/GP **2.19.1** (GMP kernel) and the new Crest C++/GMP module were run on the same CPU. The frozen panel has **20 inputs**, four implementations and **3 repetitions**, giving **240 planned rows**. Every returned integer coefficient was compared with an independent full-ring determinant recurrence. The complete schedule, calibration batches, software hashes and all outcomes are retained.

## What is being compared

The common input is a cyclic Jacobi operator supplied as integer diagonal and bond arrays. `crest-cpp` detects the least joint period and uses the cell transfer polynomial with a linear Chebyshev recurrence. `pari-dense` builds the full matrix and calls the recommended default `charpoly`, which uses modular arithmetic for integer matrices. `pari-periodic-linear` independently implements the same periodic reduction in GP. `pari-periodic-doubling` also uses binary powering identities, exploiting PARI's polynomial multiplication. These are separate implementation comparisons; the classical transfer/Floquet reduction is not a new theorem of this project.

Compute timing includes validation, period discovery, construction and complete polynomial computation on every solve; it excludes text parsing and serialization for every method. The dense baseline's matrix construction is included because it starts from the same arrays. GP's millisecond timer is used only for batches lasting at least 100 ms; the native solver uses its process CPU clock. Cold wall timing separately measures a single solve in a fresh process, including parsing and output. Compilation and independent verification are excluded. Three repetitions characterize this host and do not establish a population-level confidence interval.

## Completion

| Method | Correct rows | Timeouts | Incorrect outputs | Process errors |
|---|---:|---:|---:|---:|
| crest-cpp | 60 | 0 | 0 | 0 |
| pari-dense | 60 | 0 | 0 | 0 |
| pari-periodic-linear | 60 | 0 | 0 | 0 |
| pari-periodic-doubling | 60 | 0 | 0 | 0 |

## Paired comparisons

Ratios divide each comparator's median time by Crest's median time on the same input. **Above 1 means Crest is faster; below 1 means PARI is faster.** Only inputs correctly completed in every repetition by both implementations enter a finite ratio. The geometric mean is reported alongside the median. Broken-period cases still have cyclic nearest-neighbor structure; they are not generic dense matrices.

| Cohort | Comparator | Paired inputs | Median CPU ratio | Geometric mean CPU ratio | Median cold wall ratio | Crest CPU wins |
|---|---|---:|---:|---:|---:|---:|
| all | pari-dense | 20 | 31.893 | 37.686 | 6.556 | 20/20 |
| all | pari-periodic-linear | 20 | 2.295 | 2.370 | 4.430 | 20/20 |
| all | pari-periodic-doubling | 20 | 1.064 | 0.885 | 4.382 | 10/20 |
| periodic-springs | pari-dense | 12 | 28.506 | 35.988 | 5.722 | 12/12 |
| periodic-springs | pari-periodic-linear | 12 | 2.565 | 2.321 | 4.430 | 12/12 |
| periodic-springs | pari-periodic-doubling | 12 | 0.620 | 0.685 | 4.265 | 4/12 |
| structure-controls | pari-dense | 3 | 173.949 | 184.652 | 24.216 | 3/3 |
| structure-controls | pari-periodic-linear | 3 | 2.002 | 2.723 | 4.221 | 3/3 |
| structure-controls | pari-periodic-doubling | 3 | 0.390 | 0.579 | 4.775 | 1/3 |
| broken-period | pari-dense | 5 | 14.849 | 16.222 | 6.604 | 5/5 |
| broken-period | pari-periodic-linear | 5 | 2.159 | 2.293 | 4.597 | 5/5 |
| broken-period | pari-periodic-doubling | 5 | 2.263 | 2.116 | 4.443 | 5/5 |

## Every input

Times are median **milliseconds per complete polynomial** from the batched CPU observations. The full JSON also records the minimum and maximum of the three repetitions, and cold wall/CPU times. No case is omitted.

| Input | Sites | Found period | Maximum coefficient bits | Crest C++ | PARI dense | PARI periodic linear | PARI periodic doubling |
|---|---:|---:|---:|---:|---:|---:|---:|
| spring-m1-L24 | 24 | 1 | 31 | 0.0152 | 0.0492 | 0.0425 | 0.0235 |
| spring-m1-L96 | 96 | 1 | 130 | 0.2788 | 1.7167 | 0.4208 | 0.1000 |
| spring-m1-L384 | 384 | 1 | 529 | 6.1580 | 156.0000 | 7.6667 | 0.5167 |
| spring-m2-L24 | 24 | 2 | 58 | 0.0089 | 0.0604 | 0.0327 | 0.0271 |
| spring-m2-L96 | 96 | 2 | 237 | 0.1705 | 5.4000 | 0.2917 | 0.0975 |
| spring-m2-L384 | 384 | 2 | 954 | 4.1715 | 1084.0000 | 10.1429 | 0.7708 |
| spring-m3-L24 | 24 | 3 | 56 | 0.0089 | 0.0604 | 0.0319 | 0.0256 |
| spring-m3-L96 | 96 | 3 | 229 | 0.1396 | 5.0500 | 0.2396 | 0.0933 |
| spring-m3-L384 | 384 | 3 | 923 | 2.6559 | 1087.0000 | 7.1667 | 0.7500 |
| spring-m4-L24 | 24 | 4 | 63 | 0.0068 | 0.1108 | 0.0308 | 0.0281 |
| spring-m4-L96 | 96 | 4 | 258 | 0.1480 | 11.9000 | 0.2025 | 0.1233 |
| spring-m4-L384 | 384 | 4 | 1040 | 2.3048 | 2906.0000 | 6.7778 | 0.8625 |
| hidden-period-L192 | 192 | 3 | 517 | 0.6159 | 130.0000 | 1.2333 | 0.2400 |
| large-coefficients-L96 | 96 | 3 | 2058 | 0.2026 | 35.2500 | 1.3833 | 0.2625 |
| prime-repeat-L183 | 183 | 3 | 323 | 0.5948 | 102.0000 | 0.8783 | 0.2292 |
| single-defect-L24 | 24 | 24 | 51 | 0.0305 | 0.0883 | 0.0850 | 0.0842 |
| single-defect-L96 | 96 | 96 | 206 | 0.6229 | 9.2500 | 1.0333 | 1.0567 |
| single-defect-L384 | 384 | 384 | 829 | 11.5776 | 2107.0000 | 25.0000 | 26.2000 |
| aperiodic-L24 | 24 | 24 | 61 | 0.0310 | 0.1383 | 0.1117 | 0.0750 |
| aperiodic-L96 | 96 | 96 | 221 | 0.7008 | 22.5000 | 1.2333 | 1.1583 |

## Interpretation and limits

This panel measures synthetic one-dimensional harmonic/signed cyclic operators and exact characteristic polynomials. It does not benchmark full crystallographic space-group recognition, LLL, DFT, measured phonon data, molecular-isomer enumeration or drug discovery. The native implementation's comparison with GP also includes language/dispatch and data-representation differences; it is not a comparison with every possible hand-written libpari or gp2c implementation. The structured GP comparators are essential: a gain over a general matrix routine alone does not establish a superior periodic algorithm. Coefficient equality establishes the same finite spectral polynomial; the physical model's accuracy remains a separate question.

[Frozen panel](panel.json) · [Protocol](results-protocol.json) · [Raw observations](results.json) · [Exact reference polynomials](results-reference.json) · [Summary JSON](summary.json) · [Proof](proof.md) · [Field comparison and next challenge](opportunities.md)
