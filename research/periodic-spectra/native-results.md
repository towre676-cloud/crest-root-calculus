# Compiled PARI control

After observing the original GP comparison, we tested the concrete language/dispatch confound: both the project module and the PARI comparators now run as optimized native executables. The PARI programs link the same official PARI 2.19.1 / GMP build and implement the same validated period reduction, with linear or binary trace powering. All 20 original inputs are retained. A contemporaneous Crest rerun is the anchor. The separate protocol was frozen before its first timing; this is an explicitly post-study control, not a new independently selected confirmation panel.

## Main finding

**Compiled PARI with binary trace powering is faster on 19 of 20 paired inputs**, including 12 of 12 periodic spring inputs. Its median PARI/Crest CPU ratio is 0.206, and its geometric mean ratio is 0.263. These are ratios in the stated PARI/Crest direction; no reciprocal-median equivalence is assumed. Crest's compute wins against that comparator are: `single-defect-L384`. The periodic reduction therefore provides useful structure to exploit, but this native control does not support a general arithmetic-speed advantage over PARI.

The median cold-process wall ratio is 1.259. Cold execution includes initialization and output, so it answers a different question from repeated polynomial computation. The tables retain both costs.

All three implementations receive identical expanded integer arrays. Every timed solve repeats validation, period discovery and complete characteristic-polynomial computation. Process CPU batches last at least 100 ms; parsing, library initialization and output are included separately in cold process wall time. The native PARI source and compile command are included. No original observation has been replaced.

## Completion

| Method | Correct rows | Timeouts | Incorrect outputs | Errors |
|---|---:|---:|---:|---:|
| crest-cpp | 60 | 0 | 0 | 0 |
| pari-native-linear | 60 | 0 | 0 | 0 |
| pari-native-doubling | 60 | 0 | 0 | 0 |

## All paired comparisons

Ratios are PARI time divided by Crest time. Values above 1 favor Crest; values below 1 favor PARI. Every input must complete correctly in all three repetitions to be paired.

| Cohort | Compiled comparator | Inputs | Median CPU ratio | Geometric mean CPU ratio | Median cold wall ratio | Crest wins | PARI wins |
|---|---|---:|---:|---:|---:|---:|---:|
| all | pari-native-linear | 20 | 0.800 | 1.109 | 1.403 | 7/20 | 13/20 |
| all | pari-native-doubling | 20 | 0.206 | 0.263 | 1.259 | 1/20 | 19/20 |
| periodic-springs | pari-native-linear | 12 | 0.762 | 1.003 | 1.403 | 4/12 | 8/12 |
| periodic-springs | pari-native-doubling | 12 | 0.158 | 0.153 | 1.200 | 0/12 | 12/12 |
| structure-controls | pari-native-linear | 3 | 1.328 | 1.899 | 1.429 | 2/3 | 1/3 |
| structure-controls | pari-native-doubling | 3 | 0.156 | 0.259 | 1.133 | 0/3 | 3/3 |
| broken-period | pari-native-linear | 5 | 0.854 | 1.020 | 1.350 | 1/5 | 4/5 |
| broken-period | pari-native-doubling | 5 | 0.799 | 0.978 | 1.344 | 1/5 | 4/5 |

## Every input

Median milliseconds per full characteristic polynomial. Per-input ranges and all calibration batches remain in the JSON.

| Input | Crest C++ | Native PARI linear | Native PARI doubling |
|---|---:|---:|---:|
| spring-m1-L24 | 0.0158 | 0.0116 | 0.0026 |
| spring-m1-L96 | 0.2793 | 0.1785 | 0.0124 |
| spring-m1-L384 | 5.1754 | 6.0852 | 0.2134 |
| spring-m2-L24 | 0.0089 | 0.0062 | 0.0027 |
| spring-m2-L96 | 0.1698 | 0.1245 | 0.0160 |
| spring-m2-L384 | 3.5101 | 8.4539 | 0.4886 |
| spring-m3-L24 | 0.0076 | 0.0045 | 0.0026 |
| spring-m3-L96 | 0.1307 | 0.1056 | 0.0163 |
| spring-m3-L384 | 2.6657 | 6.3709 | 0.3947 |
| spring-m4-L24 | 0.0072 | 0.0048 | 0.0037 |
| spring-m4-L96 | 0.1137 | 0.0900 | 0.0228 |
| spring-m4-L384 | 2.3845 | 6.0625 | 0.5063 |
| hidden-period-L192 | 0.6005 | 0.7975 | 0.0935 |
| large-coefficients-L96 | 0.1955 | 1.2804 | 0.1572 |
| prime-repeat-L183 | 0.5802 | 0.4569 | 0.0810 |
| single-defect-L24 | 0.0317 | 0.0303 | 0.0253 |
| single-defect-L96 | 0.6306 | 0.5384 | 0.5256 |
| single-defect-L384 | 11.6267 | 26.5267 | 26.6334 |
| aperiodic-L24 | 0.0334 | 0.0264 | 0.0255 |
| aperiodic-L96 | 0.6962 | 0.5220 | 0.5354 |

This compiled comparison removes the main GP-language/native-language mismatch. It still compares these explicit implementations, polynomial representations, algorithms and this machine. It does not prove that either implementation is optimal or establish a result for general crystallography or chemistry.

[Protocol](native-results-protocol.json) · [Raw results](native-results.json) · [Summary](native-summary.json) · [Native source](pari_native.c) · [Build provenance](pari-native-build.json) · [Original GP study](results.md)
