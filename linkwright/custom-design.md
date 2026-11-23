# Build a custom four-bar with Linkwright

Open [the workshop](index.html#designer). To fabricate the fixed hand-cranked pointer and follow its hardware, dimensions and bench criteria, use [the reference assembly guide](kit/assembly-guide.html).

## Design a motion

**Flap a wing.** Choose the swing, the longer-to-shorter stroke-time ratio, and the rocker length, then select **Design the mechanism**. Stroke times assume a constant crank speed. The page lists the crank angle for each direction of rocker travel. Switch **Drive direction** to move the fast stroke to the other direction. Which one moves an attached wing downward depends on the way you mount the arm. The pale wing silhouette shows a possible attachment.

**Trace a path.** Draw a closed loop or choose the built-in wingtip, set its horizontal width, and choose **Find a mechanism**. The purple line is the target. The coloured dot is the tracing point. The search repeats deterministically for a given drawing. Read the mean and largest sampled miss before committing material. Changing the target width starts a new design. Changing the overall scale scales both targets and parts, followed by fresh rounding and error measurement.

**Inspect a pose.** Move the crank-angle slider to pause at a chosen input angle. Zero points along the base from the fixed crank pivot O to the fixed rocker pivot D. The positive branch places the moving joint B on the left when looking from crank tip A toward D in the plan view. The dashed second assembly in the branch diagram moves differently.

## Set a fabrication tolerance

Measure pin spacings **from hole centre to hole centre**. The tolerance slider allows an independent change in each of the four lengths. If each centre can miss its mark by at most 0.25 mm in the plane, two errors pointing in opposite directions can change that centre spacing by 0.50 mm. Separate radial joint play, tilt and flexible material can introduce additional motion error.

The custom builder scales the ideal lengths and then rounds the four spacings and coupler tracing offsets to **0.1 mm**. All finalized poses, checks, drawings and design exports use those cut values. The base keeps its original direction and anchor, so world coordinates can contain additional decimals even though the fixed-hole spacing is on the tenth-millimetre grid.

**Worked boundary.** Select 60°, 1.80× and 70 mm, cardboard, and a ±5.00 mm spacing tolerance. At a 159% scale the cut lengths are 93.8, 46.6, 84.1 and 111.3 mm for ground, crank, coupler and rocker. Their folded margin is 20.0 mm. Four adverse 5 mm errors use all 20.0 mm and permit a toggle. **160% is the first whole-percent cut scale that passes** under these stock assumptions. The fixed 159% plywood reference uses a much tighter shop tolerance.

## Cut, assemble and check

1. Select stock and a pin size. Inspect the full-turn, stock and crank–rocker clearance results. **Review before cutting** flags a failed or unresolved choice and is carried onto the exported sheet.
2. Save **parts SVG** and **design JSON** together. Open the SVG in a print or vector application and print at 100%. Use printing across pages if necessary. Both scale rules must measure 50 mm. The coloured outside contour of the coupler plate includes the edge material around the holes. On a nearly collinear drawing the third tracing hole can lie close to a pivot; inspect actual spacing and material.
3. Lay out the ground, crank, coupler and rocker. The base template puts the fixed pivots horizontally for cutting. Rotate the board to the specified world bearing to reproduce the on-screen orientation. Use the assembly side printed on the sheet.
4. Use **Check measured pin spacings** to compare first-article geometry with the requested motion. Turn the unloaded linkage slowly in both directions. Check the complete assembly, including the tracing attachment or wing, fasteners and base. The current custom sheet sizes the planar parts. Designing custom bearings, axial stacks, loads and attachments remains part of building a different mechanism.

## What the reported numbers measure

For lengths $g,a,b,c$ in ground–crank–coupler–rocker order, full-turn closure requires both $|g-a|-|b-c|$ and $b+c-g-a$ to be strictly positive. Over a tolerance radius $t$, the two worst margins are $\max(0,|g-a|-2t)-(|b-c|+2t)$ and $b+c-g-a-4t$. A small numerical guard classifies a value effectively at zero as **at limit**. A stock must accommodate the maximum measured span, strip width and spacing tolerance.

At a rocker reversal the crank and coupler align. The distance from O to B is then $a+b$ or $|a-b|$. Intersect each corresponding circle with the rocker circle, and retain the two valid stops on the selected branch. The two input-angle intervals between them determine the stroke-time ratio. The minimum transmission angle follows from the extreme crank-tip distances $|g-a|$ and $g+a$. These angle results come from the fitted rigid geometry. Force, friction and material deformation require physical measurements.

The path score samples each direction between the target polyline and the simulated path. Each direction has 192 equal-arclength samples; the simulated curve has 720 poses. The largest sampled miss is a diagnostic. Timing, traversal order at self-intersections and error between samples need separate evaluation.

## Reuse a design export

Version 3 JSON carries millimetres and coordinate conventions, the actual cut spacings, signed assembly branch, along-coupler and perpendicular tracing offsets, and the current stock and tolerance settings. The top-level swing and timing give the realized rocker values when the selected branch has two distinct stops. The `requested` object preserves the intended motion or path. The `motion.stops` inputs and rocker angles are measured counterclockwise from O→D. In a path export the target coordinates and error distances are at the **current cut scale** in the same world frame as the fixed pivot coordinates.
