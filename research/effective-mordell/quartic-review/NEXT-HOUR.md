# One-hour research split

## Claude — reconstruction

Work from commit **093779c**. Keep the old receipt for comparison.

| Minutes | Work | Evidence to return |
|---|---|---|
| 0–25 | Fix the one-variable pole chart. Move the univariate branch ahead of the plane dimension test. Handle the zero ideal, unit ideal, genuine finite roots and roots supported on the degeneracy locus. | The original PERMUTED sextic and its rational transform with ascending coefficients \([1,0,-26,-70,-57,-14,-1]\) each give two maps. Transport the maps through \(x=1+1/u\) and compare. |
| 25–40 | Use a primitive integral multiple of the source to generate integral branch invariants. Preserve the original coefficients for the cover identity. | Show that \(f(x)=\prod_{j=0}^{5}(x-j/2)\) yields fifteen distinct rational pairings. Verify agreement under denominator clearing without changing the source twist. Make polynomial content nonnegative in rational_content and add the positive nonstable-pairing control from the review. |
| 40–60 | Make distinctness across pairing orbits exact at acceptance. Use a common separating linear form with pairwise coprime minimal polynomials of the right degrees, or evaluate the paired invariant exactly against its proposed resolvent factor. Rerun both verifiers and the 36-neighbour sweep. | Commit or patch, targeted regressions, exact-orbit certificate, versioned receipts and any unfinished runs. Document the actual stopping conditions. |

Use the supplied positive control on y² = x⁶ + 1. Check six quadratic map orbits with four over a genuinely nonstable pairing. The signed 16-point graph over extensions can follow these correctness fixes.

## Codex — arithmetic and external comparison

| Minutes | Work | Evidence to return |
|---|---|---|
| 0–15 | Starting from the completed Magma control and equation-only run on the \(d=26\) twist, extend the validation to the fixed set of positive squarefree twist parameters from 31 through 100. Retain every parameter. Give priority for direct Magma comparisons to cases with Jacobian rank at least two and an elliptic factor of rank zero. | Record both rank intervals, curve count, fields, rational point sets and available completeness flags. Keep timeouts distinct from false flags. |
| 15–40 | Start the positive-rank \(d=22\) case. Exploit the two elliptic factor classes as the extra rational divisor classes needed for quadratic Chabauty. Work at 5, where the curve has good reduction and the two elliptic Frobenius traces are \(-4\) and \(4\). | Rational points already include \((0,0),(1,\pm44),\infty\). Six residue disks at 5 require examination, including the disks containing known points. Record any rigorous exclusions and the outstanding local height or sieve calculations. |
| 40–60 | Assemble the family results and run the audited reconstruction patch on the frozen examples if available. Prioritize a held-out genus-two curve with a genuinely necessary nonrational quotient: a quotient available over a field \(L\), finite \(E(L)\), and no rational rank-zero quotient giving the same answer by the existing route. | Exact fibres, projective points, fields and rank bounds, with observed baseline outcomes separate from prior-work descriptions. A failed search or missing rank certificate is an open test. |

The quadratic Chabauty rank criterion supplies a reason to try \(d=22\). The two independent elliptic quotient maps give Jacobian rank two and Néron–Severi rank at least two. Computing the finite set of \(p\)-adic candidates and proving which are rational are additional steps. [Balakrishnan and Müller](https://math.bu.edu/people/jbala/2020BalakrishnanMuellerNotes.pdf)
