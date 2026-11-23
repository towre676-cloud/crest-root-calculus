# Part II · Certificates for spectra and electricity networks

**A certificate tells someone else how to check a particular answer.** For one mineral mixture, you can hand a reader just two wavelengths. For one solar offer, you can hand them the row of the voltage calculation where all the headroom is used up. What the reader learns is different in the two cases.

## 1. Two wavelengths that disagree

The first reading on a measured calcite–clay powder allows about **33–38%** calcite in a straight-line model of reflected light. Another wavelength requires **at least 97%**. Either band by itself looks possible. Their allowed intervals have no point in common, so that model fails at the stated tolerance. The lab weighed **80% calcite** into this powder. Those two light readings can challenge a model even when the ingredients are known. [Measured data, diagram and proof](spectra/README.md#1-a-measured-two-band-exclusion)

**Correction to the grain-size interpretation.** In a model for light scattered among grains, the weak-absorption signal identifies two products of grain-size scale and optical fraction. Different surface-area shares can produce those same products. Knowing the effective sizes and densities of the reference grains makes the mass ratio constant along that family of solutions. A scan with uncertain reference sizes still gives an uncertain mass fraction. The precision also depends on where the instrument noise enters the calculation. [Weak-absorption derivation and numerical experiment](spectra/README.md#2-grain-size-and-the-weak-absorption-fibre)

## 2. Eighteen houses sharing one voltage limit

For a balanced radial network, raising injections raises the linear upper bound on squared voltage. The argument evaluates the worst corner of an entire box of possible injections and accounts for the network’s shunts. On the released feeder, **2.71964344 kW per new home** satisfies all 42 upper-voltage rows for source magnitude at most 1.02 per unit and the stated active and reactive injection bounds. One row also proves the cap is the largest equal offer this screen accepts. A one-line counterexample shows an AC-feasible injection can fail the screen. [Assumptions, diagram, rational bounds and power-flow replay](power/README.md)

**Correction to the volt-var result.** At a supplied local gain of 1.256, an undamped discrete update multiplies a small error by **−1.256**. It crosses the target and grows. Moving halfway to the commanded setting multiplies that mode by **−0.128**. This proves a local stability distinction for the specified equations. The released data establish neither the supplied control-curve gain nor the draft’s zero-curtailment outcome in an installed inverter. [Control law, stability thresholds and provenance](power/README.md#3-when-a-voltage-controller-overshoots)
