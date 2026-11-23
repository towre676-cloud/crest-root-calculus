# Measured solver comparisons

The three frozen panels contain 51 inputs and 312 planned runs in the published studies. Returned answers are checked exactly, and the map study uses a frozen verified reference. Wrong answers and timeouts remain in the recorded schedule. CPU seconds include each worker and its reaped native children; the adapter descriptions identify the work included.

The quartic study resumed after 55 recorded runs without replacing any earlier measurement. Its input and adapter were frozen before the original start; imported research-module hashes were captured at resumption and verified at completion. The initial protocol stored the available affinity set, but not the selected worker CPU; the resumption records CPU eight. The later studies froze their mathematical dependencies before timing. The original 12-input counting panel and its 48 runs are preserved separately; the final counting panel retains those inputs and adds a predeclared three-input larger-field cohort, with every case retimed. A first map attempt stopped on an exact disagreement; its 11 recorded rows and original protocol are preserved separately. The published map study reran the unchanged input panel with verified references and a separately labelled corrected dense comparator. See the [protocol and reproduction notes](README.md).

## Marked quartic pencils

Each input was run twice per method with a 30-second process limit. Entries below give median solver CPU seconds when both runs completed correctly. T records a timeout. W marks an incorrect answer and retains its CPU time in parentheses. E would mark a process error. A mixed entry preserves both outcomes in repetition order.

| Input | Exact outcome | crest | singular | sympy-buchberger |
|---|---|---:|---:|---:|
| square-1 | 1 geometric pencil | 2.8069 | 20.0805 | T / T |
| square-2 | 1 geometric pencil | 2.5491 | 11.6464 | T / T |
| square-3 | 1 geometric pencil | 2.2064 | 10.5215 | T / T |
| square-4 | 1 geometric pencil | 1.7292 | 3.9044 | T / T |
| square-5 | 1 geometric pencil | 2.3837 | 14.5896 | T / T |
| square-6 | 1 geometric pencil | 2.5242 | 17.1404 | T / T |
| random-1 | 0 geometric pencils | 1.7782 | 6.4624 | T / T |
| random-2 | 0 geometric pencils | 1.3275 | 2.2902 | T / T |
| random-3 | 0 geometric pencils | 1.5690 | 2.4175 | T / T |
| random-4 | 0 geometric pencils | 1.4848 | 2.6494 | T / T |
| random-5 | 0 geometric pencils | 3.0267 | 28.9110 | T / T |
| random-6 | 0 geometric pencils | 3.2206 | T / T | T / T |
| random-7 | 0 geometric pencils | 4.4080 | T / T | T / T |
| random-8 | 0 geometric pencils | 3.4175 | 25.5363 | T / T |
| fermat | 0 geometric pencils | 0.1180 | 0.0942 | 0.5264 |
| klein | 0 geometric pencils | 0.0930 | 0.0812 | 0.4126 |
| cyclic | 1 geometric pencil | 0.3576 | 0.1797 | 2.6311 |
| diagonal-2 | 0 geometric pencils | 0.1269 | 0.0938 | 0.5377 |

| Method | Correct inputs completed twice | Timeout runs | Wrong-answer runs | Process-error runs |
|---|---:|---:|---:|---:|
| crest | 18/18 | 0 | 0 | 0 |
| singular | 16/18 | 4 | 0 | 0 |
| sympy-buchberger | 4/18 | 28 | 0 | 0 |

Paired ratios use only inputs that both methods completed correctly in both repetitions. Each ratio divides the comparator’s median time by the project method’s median time on the same input. CPU ratios use solver CPU; process-wall ratios also include interpreter startup. Both use the same jointly completed inputs. A ratio above one means the comparator took more time. Timeouts and incorrect outputs contribute to their respective table columns and remain outside these finite ratios.

| Comparator | Paired inputs | Median CPU ratio | Geometric mean CPU ratio | Median process-wall ratio | Project faster on CPU |
|---|---:|---:|---:|---:|---:|
| singular | 16 | 2.946 | 2.611 | 2.672 | 12/16 |
| sympy-buchberger | 4 | 4.448 | 4.984 | 2.092 | 4/4 |

## Constant-vertical hyperelliptic maps

Each input was run twice per method with a 10-second process limit. Entries below give median solver CPU seconds when both runs completed correctly. T records a timeout. W marks an incorrect answer and retains its CPU time in parentheses. E would mark a process error. A mixed entry preserves both outcomes in repetition order.

| Input | Exact outcome | crest | singular | sympy-decomposition | corrected-dense-decomposition |
|---|---|---:|---:|---:|---:|
| g1-d1-cover | 2 geometric, 2 rational | 0.0629 | 0.0986 | 0.0568 | 0.0670 |
| g1-d1-twist | 2 geometric, 0 rational | 0.0649 | 0.0846 | 0.0578 | 0.0594 |
| g1-d1-perturbed | 0 geometric, 0 rational | 0.0487 | 0.0705 | 0.0547 | 0.0427 |
| g1-d3-cover | 2 geometric, 2 rational | 0.0761 | 0.1441 | W (0.0550) / W (0.0342) | 0.0666 |
| g1-d3-twist | 2 geometric, 0 rational | 0.0617 | 0.1085 | W (0.0342) / W (0.0649) | 0.0734 |
| g1-d3-perturbed | 0 geometric, 0 rational | 0.0501 | 0.0944 | 0.0480 | 0.0387 |
| g1-d5-cover | 4 geometric, 2 rational | 0.1014 | 0.2113 | W (0.0608) / W (0.0419) | 0.0901 |
| g1-d5-twist | 4 geometric, 0 rational | 0.0921 | 0.2312 | W (0.0455) / W (0.0422) | 0.0855 |
| g1-d5-perturbed | 0 geometric, 0 rational | 0.0582 | 0.1257 | 0.0466 | 0.0421 |
| g2-d3-cover | 2 geometric, 2 rational | 0.0707 | 0.2196 | W (0.0493) / W (0.0340) | 0.0951 |
| g2-d3-twist | 2 geometric, 0 rational | 0.0954 | 0.1500 | W (0.0523) / W (0.0336) | 0.0911 |
| g2-d3-perturbed | 0 geometric, 0 rational | 0.0497 | 0.1477 | 0.0414 | 0.0602 |
| g3-d5-cover | 2 geometric, 2 rational | 0.3029 | 2.4201 | W (0.0558) / W (0.0449) | 0.2053 |
| g3-d5-twist | 2 geometric, 0 rational | 0.1768 | 1.3879 | W (0.0557) / W (0.0655) | 0.1335 |
| g3-d5-perturbed | 0 geometric, 0 rational | 0.1848 | 1.2792 | 0.0510 | 0.0600 |
| g4-d7-cover | 2 geometric, 2 rational | 0.9296 | T / T | W (0.0880) / W (0.0579) | 0.8720 |
| g4-d7-twist | 2 geometric, 0 rational | 0.4745 | T / T | W (0.0759) / W (0.0602) | 0.3417 |
| g4-d7-perturbed | 0 geometric, 0 rational | 0.4375 | T / T | 0.0842 | 0.0691 |

| Method | Correct inputs completed twice | Timeout runs | Wrong-answer runs | Process-error runs |
|---|---:|---:|---:|---:|
| crest | 18/18 | 0 | 0 | 0 |
| singular | 15/18 | 6 | 0 | 0 |
| sympy-decomposition | 8/18 | 0 | 20 | 0 |
| corrected-dense-decomposition | 18/18 | 0 | 0 | 0 |

Paired ratios use only inputs that both methods completed correctly in both repetitions. Each ratio divides the comparator’s median time by the project method’s median time on the same input. CPU ratios use solver CPU; process-wall ratios also include interpreter startup. Both use the same jointly completed inputs. A ratio above one means the comparator took more time. Timeouts and incorrect outputs contribute to their respective table columns and remain outside these finite ratios.

| Comparator | Paired inputs | Median CPU ratio | Geometric mean CPU ratio | Median process-wall ratio | Project faster on CPU |
|---|---:|---:|---:|---:|---:|
| singular | 15 | 2.084 | 2.559 | 1.285 | 15/15 |
| sympy-decomposition | 8 | 0.862 | 0.647 | 0.963 | 1/8 |
| corrected-dense-decomposition | 18 | 0.883 | 0.780 | 0.935 | 4/18 |

## Hyperelliptic local zeta numerators

Each input was run twice per method with a 10-second process limit. Entries below give median solver CPU seconds when both runs completed correctly. T records a timeout. W marks an incorrect answer and retains its CPU time in parentheses. E would mark a process error. A mixed entry preserves both outcomes in repetition order.

| Input | Exact outcome | crest-direct | pari |
|---|---|---:|---:|
| g2-d5-p5 | degree 4; max field 25 | 0.0180 | 0.0405 |
| g2-d5-p7 | degree 4; max field 49 | 0.0181 | 0.0373 |
| g2-d5-p11 | degree 4; max field 121 | 0.0226 | 0.0429 |
| g3-d7-p5 | degree 6; max field 125 | 0.0236 | 0.0394 |
| g3-d7-p7 | degree 6; max field 343 | 0.0241 | 0.0403 |
| g3-d7-p11 | degree 6; max field 1331 | 0.0215 | 0.0446 |
| g4-d9-p5 | degree 8; max field 625 | 0.0460 | 0.0454 |
| g4-d9-p7 | degree 8; max field 2401 | 0.1129 | 0.0502 |
| g4-d9-p11 | degree 8; max field 14641 | 0.0291 | 0.0409 |
| g2-d6-p5 | degree 4; max field 25 | 0.0209 | 0.0538 |
| g3-d8-p7 | degree 6; max field 343 | 0.0211 | 0.0421 |
| g4-d10-p11 | degree 8; max field 14641 | 0.0301 | 0.0539 |
| g2-d5-p1009 | degree 4; max field 1018081 | 0.3632 | 0.0810 |
| g3-d7-p97 | degree 6; max field 912673 | 0.2996 | 0.0548 |
| g4-d9-p31 | degree 8; max field 923521 | 0.4067 | 0.0469 |

| Method | Correct inputs completed twice | Timeout runs | Wrong-answer runs | Process-error runs |
|---|---:|---:|---:|---:|
| crest-direct | 15/15 | 0 | 0 | 0 |
| pari | 15/15 | 0 | 0 | 0 |

Paired ratios use only inputs that both methods completed correctly in both repetitions. Each ratio divides the comparator’s median time by the project method’s median time on the same input. CPU ratios use solver CPU; process-wall ratios also include interpreter startup. Both use the same jointly completed inputs. A ratio above one means the comparator took more time. Timeouts and incorrect outputs contribute to their respective table columns and remain outside these finite ratios.

| Comparator | Paired inputs | Median CPU ratio | Geometric mean CPU ratio | Median process-wall ratio | Project faster on CPU |
|---|---:|---:|---:|---:|---:|
| pari | 15 | 1.675 | 1.021 | 1.011 | 10/15 |

### Predeclared counting cohorts

The small-prime cohort contains the twelve original inputs. The larger-field cohort contains three seeded inputs with maximum extension orders 1009² = 1,018,081, 97³ = 912,673 and 31⁴ = 923,521. The two-million-element table budget applies to the direct counter. Both cohorts use the same ten-second process and four-GiB memory limits. Cohort membership was fixed before the final timing study.

| Cohort | Direct correct twice | PARI correct twice | Paired inputs | Median PARI/direct CPU ratio | Median PARI/direct process-wall ratio |
|---|---:|---:|---:|---:|---:|
| small-prime | 12/12 | 12/12 | 12 | 1.8442 | 1.0807 |
| large-field | 3/3 | 3/3 | 3 | 0.1831 | 0.5364 |
