# The numbers behind the origin story

You can check the first two puzzles with pencil and paper. The later figures use a prime product, a measured family of curves, and two exact solver records. [Return to the illustrated story.](../../README.md#how-the-project-began)

## Two clocks, two curves

Take the starting pair $(a,b)=(1,1)$ modulo 5. Its 25 possible refinements modulo 25 are $(1+5i,1+5j)$ for $i,j$ from 0 to 4. Each little grid on the main page uses $i$ from left to right and $j$ from bottom to top.

Try one step right and two up. That gives $(x,y)=(6,11)$. For the parabola, $6^2=36$ leaves remainder 11 after division by 25. The dot is filled. For the cusp, $11^2=121$ leaves remainder 21 and $6^3=216$ leaves remainder 16. That dot stays open.

Now try the block above $(0,0)$. Here $x=5i$ and $y=5j$, so $y^2=25j^2$ and $x^3=125i^3$. Both sides of the cusp equation leave remainder 0 for every choice of $i$ and $j$. For the parabola, $y=x^2$ requires $y$ to leave remainder 0, which picks out one row.

| Curve | Starting pairs modulo 5, in the order drawn | Number of refinements per pair | Total modulo 25 |
|---|---|---|---:|
| $y=x^2$ | (0,0), (1,1), (2,4), (3,4), (4,1) | 5, 5, 5, 5, 5 | 25 |
| $y^2=x^3$ | (0,0), (1,1), (1,4), (4,2), (4,3) | 25, 5, 5, 5, 5 | 45 |

The two real-curve sketches show the ordinary shapes near zero. Their little grids count residue pairs. Both totals were checked by testing every one of the $25^2=625$ pairs against each equation. The starting experiment appears in *Root–Crest Calculus: exact prototype and research specification*, under “Foundations and exact obstruction.”

## The 80-slot chain

Squaring 0, 1, 2, 3, and 4 and keeping the remainder modulo 5 gives 0, 1, 4, 4, and 1. For neighboring squares to add to 1 modulo 5, the only allowed ordered pairs are **(0,1), (0,4), (1,0), and (4,0)**.

| Number of slots | Ends in 0 | Ends in 1 | Ends in 2 | Ends in 3 | Ends in 4 | Total |
|---:|---:|---:|---:|---:|---:|---:|
| 1 | 1 | 1 | 1 | 1 | 1 | 5 |
| 2 | 2 | 1 | 0 | 0 | 1 | 4 |
| 3 | 2 | 2 | 0 | 0 | 2 | 6 |
| 4 | 4 | 2 | 0 | 0 | 2 | 8 |

There’s also a direct proof for 80 slots. Pick whether the zeros occupy the odd slots or the even slots. That gives two choices. Each of the 40 remaining slots can independently hold either 1 or 4. The total is $2\times2^{40}=2{,}199{,}023{,}255{,}552$.

The original program tests $5\times5=25$ possible neighbor pairs at each of 79 joins. That accounts for the **1,975 pair checks**. It measures this implementation on a chain. Constraints that link far-apart slots can require keeping track of more information. The original research specification describes the counting algorithm under “Sparse high-dimensional computation.”

## Prime by prime

Let $p$ be a prime. A projective line has $p+1$ points over the field with $p$ elements. Divide by $p$ to get the factor $(p+1)/p$. A pair of projective lines has $(p+1)^2$ points. Dividing by $p^2$ gives the square of the line’s factor.

| Include primes through | New line factor | Product for one line | Product for two lines paired |
|---:|---:|---:|---:|
| 2 | $3/2$ | $3/2$ | $9/4$ |
| 3 | $4/3$ | $2$ | $4$ |
| 5 | $6/5$ | $12/5$ | $144/25$ |
| 7 | $8/7$ | $96/35$ | $9{,}216/1{,}225$ |

The graph draws these products for every prime up to 31. The step stays level between primes. At any cutoff $x$, the surface’s product is exactly the square of the line’s product. A separate theorem gives the latter’s long-run growth as a constant times $\log x$. Squaring gives the exponent 2 for the surface. [See the surface and the growth argument.](README.md#the-split-quadric)

## Checking the apparent law

Each dot is a related version of the elliptic curve 11a1, called a *quadratic twist*. The audit in §15 of the monograph uses the 201 rank-zero records with conductor $N>200$ and at least two listed positive zero heights. Both plots use the same records, horizontal positions, and axis scales. On the left, the vertical position is the whole measured score. On the right, it is the score after subtracting the nearest pair’s contribution. [Inspect all 201 plotted records (CSV).](data/lowest-zero-audit.csv)

To reproduce the axes, let $\gamma_1$ be the smallest listed positive zero height and put $L'=\log(\sqrt{N}/(2\pi))$. With $D_{\rm lo}$ the historical numerical lower estimate used in the audit, the quantities shown are

$$
x=\frac{2}{L'\gamma_1\sqrt{1+\gamma_1^2}},\qquad
y=\frac{D_{\rm lo}}{L'},\qquad
\text{remainder}=y-x.
$$

The highlighted record is twist $d=-2531$. Its whole score rounds to 10.73, the pair contributes about 9.79, and the remainder rounds to 0.94. The dashed line on the left is the fit to all 201 points. Its slope is about 0.997 and its $R^2$ is about 0.989. Because $x$ is already a piece of $y$, a close fit alone gives little evidence of a new prediction. What is still worth studying is the variation in the remainder.

The earlier claim in §12 used a different fit. This figure reproduces the later diagnostic using the archived measurements, including their original numerical precision. The archive inputs are `crest/twist_zeros_fine.json` and `crest/twist_resid.py`.

## When the bounds meet

In both number-line panels, an open circle marks an integer between the recorded lower and upper bounds. For instance, 16 candidates can allow totals from **0 through 16**, which gives **17** possible totals on the line. The last diamond marks the single remaining total.

| Input curve | Step in the decision record | Lower bound | Upper bound |
|---|---|---:|---:|
| C⋆ | First local precision at 5 | 0 | 1 |
| C⋆ | Second local precision and exact-map check | 1 | 1 |
| $F+13D=0$ | First local precision at 13 | 0 | 16 |
| $F+13D=0$ | Second local precision, modulo $13^2=169$ | 0 | 0 |

Here $F=X^4+Y^4+Z^4$ and $D=-5X^2Z^2-4Y^2Z^2$. The bounds count geometric degree-three quotient classes. For C⋆, the scan over the field with 25 elements accounts for the possible classes, and a verified exact formula supplies a class. For the perturbed quartic, every local candidate fails at the next precision. The field with 25 elements used in the C⋆ scan and the precision of a calculation modulo 25 play different mathematical roles.

The table gives the `lower` and `upper` entries from the two `levels` arrays in `crest/scan_completeness_xxx/decision/cstar_decision.json` and `nonlift_decision.json`. The earlier figures show the [C⋆ scan](README.md#the-c-star-scan) and the [Fermat lifting test](README.md#the-fermat-lifting-test).
