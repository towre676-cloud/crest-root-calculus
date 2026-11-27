# The mathematics behind the figures

The section numbers and receipt paths refer to the integrated Part XXX research archive. Open an illustration on its own to see the animation, or choose the still version. [Project overview](../../README.md) · [Autoplay exhibits](../visuals/README.md)

| Illustration | Still version | Equations and data |
| --- | --- | --- |
| [The cusp and the parabola](cusp-lifts.svg) | [Still](still/cusp-lifts.svg) | [Two clocks, two curves](origins.md#two-clocks-two-curves) |
| [Counting an 80-slot chain](chain-count.svg) | [Still](still/chain-count.svg) | [The neighbour rule and exact count](origins.md#the-80-slot-chain) |
| [The split quadric](quadric-rulings.svg) | [Still](still/quadric-rulings.svg) | [Point counts on the surface](#the-split-quadric) |
| [How the prime product grows](prime-steps.svg) | [Still](still/prime-steps.svg) | [Prime by prime](origins.md#prime-by-prime) |
| [Auditing the lowest-zero fit](lowest-zero-audit.svg) | [Still](still/lowest-zero-audit.svg) | [The 201 records and what the fit measured](origins.md#checking-the-apparent-law) |
| [A fibre of the Klein quartic](klein-fibre.svg) | [Still](still/klein-fibre.svg) | [The degree-three map](#the-klein-quartic) |
| [The C-star scan](cstar-modulo-five.svg) | [Still](still/cstar-modulo-five.svg) | [Residue counts and marked points](#the-c-star-scan) |
| [The Fermat lifting test](lifting-depth.svg) | [Still](still/lifting-depth.svg) | [The perturbation equations](#the-fermat-lifting-test) |
| [When the bounds meet](proof-bounds.svg) | [Still](still/proof-bounds.svg) | [The two decision records](origins.md#when-the-bounds-meet) |

## The Klein quartic

The plot uses the real curve $F=X^3Y+Y^3Z+Z^3X=0$. We choose coordinates with $X+Y+Z=1$, so specifying $X$ and $Y$ also specifies $Z=1-X-Y$. Cycling the coordinates from $(X,Y,Z)$ to $(Y,Z,X)$ preserves the equation. Applying the cycle twice gives $(Z,X,Y)$.

The quotient map uses two quantities that have the same values at all three points,

$$
u=XY+YZ+ZX,\qquad W=14XYZ+u^2-5u+1.
$$

| Input $(X,Y,Z)$ | $u$ | $XYZ$ | Output $(u,W)$ |
|---|---:|---:|---|
| $(1,0,0)$ | 0 | 0 | $(0,1)$ |
| $(0,1,0)$ | 0 | 0 | $(0,1)$ |
| $(0,0,1)$ | 0 | 0 | $(0,1)$ |

The blue curve has equation $W^2=Q(u)$, where $Q(u)=1-10u+27u^2-10u^3-27u^4$. You can verify what happens away from the three marked vertices too. In the displayed chart the polynomial identity is

$$
W^2-Q(u)=28\,(X^3Y+Y^3Z+Z^3X)(X^3Z+Z^3Y+Y^3X).
$$

The first factor vanishes at every point of the green quartic. Section 109 identifies the map as the quotient by the order-three cycle. The animation plots the three cyclic permutations of each sample and evaluates their common image at that same sample. The lines between consecutive screen positions approximate the motion.

**Archive source** · Part XXV, §109, in `crest/quadratic_fibre_xxv/RESEARCH_NOTE.md`.

## The C-star scan

To reproduce a cell, take these coefficients in arithmetic modulo 5. For each row, multiply the expression by the indicated power of $Z$, then add the five contributions.

| Power of $Z$ | Coefficient |
|---|---|
| $Z^0$ | $X^4+4X^3Y+3X^2Y^2+4XY^3+2Y^4$ |
| $Z^1$ | $2X^3+2X^2Y+3Y^3$ |
| $Z^2$ | $2X^2+2Y^2$ |
| $Z^3$ | $4X+Y$ |
| $Z^4$ | $1$ |

For the ring in the ordinary grid, substitute $(1,0,1)$. The five contributions are $1,2,2,4,1$. Their sum is 10, giving remainder 0. For the ring at infinity, substitute $(2,1,0)$. Only the first row contributes, and its value is $16+32+12+8+2=70$, again giving remainder 0.

Over the field with five elements, every projective position has exactly one of the following representatives. This explains why the drawing covers the entire plane.

| Coordinates tested | Number tested | Number on the quartic |
|---|---:|---:|
| $(x,y,1)$ | 25 | 6 |
| $(x,1,0)$ | 5 | 2 |
| $(1,0,0)$ | 1 | 0 |

The projection and contact equations mark two of the eight curve points. Part XXX proves that for this example the field with 25 elements suffices to see every geometric class in the reduction. The scan over that field checks 40 curve points, finds the same two marks, and groups them into one class. This gives an upper bound of one in characteristic zero. The program reconstructs a map and verifies it exactly, attaining the upper bound.

**Archive sources** · Part XXX, §§138–143. Under `crest/scan_completeness_xxx/`, see `scan/cstar.json` for the integer polynomial, `scan/cstar_25.json` for the larger scan and `decision/cstar_decision.json` for the field bound and the certified result.

## The Fermat lifting test

Let $F=X^4+Y^4+Z^4$ and $D=-5X^2Z^2-4Y^2Z^2$. Compare $F=0$ with $F+13^kD=0$ for $k=1,2,3$. Reducing the coefficients modulo $13^k$ erases the added term. Passing from $13^k$ to $13^{k+1}$ reveals one more digit.

| Quartic | Last displayed modulus with 16 candidates | First modulus with zero candidates | Exact geometric classes |
|---|---:|---:|---:|
| $F=0$ | 28,561 | — | 16 |
| $F+13D=0$ | 13 | 169 | 0 |
| $F+13^2D=0$ | 169 | 2,197 | 0 |
| $F+13^3D=0$ | 2,197 | 28,561 | 0 |

The calculation starts from a complete set of 16 candidate classes at 13, each with two marked points. It tests whether the equations for each class can still be satisfied at the next precision. For the three perturbations, the obstruction eliminates every class at the indicated step. For the unperturbed Fermat quartic, the 16 degree-three quotients constructed in characteristic zero supply the matching lower bound.

**Archive sources** · Part XXX, §§142–144. Under `crest/scan_completeness_xxx/`, see `adic/adic_receipt.json`, `hensel/deformation_receipt.json` and `decision/nonlift_decision.json`.

## The split quadric

The surface in the figure is the product of two projective lines. Its equation in four projective coordinates is $X_0X_3=X_1X_2$. On the patch $X_0=1$ we can set $X_1=u$, $X_2=v$, $X_3=uv$. The drawing projects a square piece of this real saddle onto the page. Holding $u$ or $v$ constant gives a straight line. Every point of the patch lies on one line from each family.

Over the field with $p$ elements, each projective line has $p$ finite positions and one position at infinity. The projective surface includes all pairs of positions, including the boundary missing from the saddle patch.

| Pair of positions | Number | At $p=5$ |
|---|---:|---:|
| Both finite | $p^2$ | 25 |
| Only the first is infinite | $p$ | 5 |
| Only the second is infinite | $p$ | 5 |
| Both infinite | $1$ | 1 |
| **Total** | $(p+1)^2$ | **36** |

To see the crest, divide each count by $p^2$ and multiply over the primes up to $x$. This gives

$$
\prod_{p\le x}\frac{(p+1)^2}{p^2}
=\left(\prod_{p\le x}(1+1/p)\right)^2
\sim C(\log x)^2
$$

for a positive constant $C$. The last step follows from the prime product formula of Mertens and the convergent product of $1-1/p^2$. The exponent of $\log x$ is 2.

**Archive source** · Part I, §§1–3 of `crest/index.html`, for the point-count and growth-exponent definitions. The drawing uses the explicit split model $X_0X_3=X_1X_2$.
