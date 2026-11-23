# Which home upgrade earns more heat savings per unit of cost?

A home that loses 352 watts per degree could have leaky joints, poorly insulated surfaces, or both. Enter a measured range for the whole home and two retrofit quotes. The program calculates whether either quote buys more annual useful heat savings per unit of installation cost in every compatible home. When the answer depends on the source of heat loss, it ranks follow-up surveys by how much each can narrow the decision.

Run the example in this folder with Python 3.9 or later.

```sh
python app.py
python app.py example.json
python app.py surveyed.json
python field_replay.py
python -m unittest -v test_app.py
```

## Worked decision

DEEP measured home 07LT at **352 ±12 W/K** before renovation. The example takes that reported range as the accepted input. The example budget, climate, upgrade fractions and survey specifications are illustrative. The assumed fabric work saves one fifth of fabric heat loss at a cost of 8,000 units. The air-sealing quote saves 35% of unintended air loss for 4,000. A steady total alone leaves both ahead in different feasible homes. A targeted infiltration survey with ±8 W/K error lowers the worst possible width of the difference in value from **1.9656 to 0.1152** annual kWh per unit of installation cost. A second whole-house test with ±6 W/K still permits opposite splits at exactly the same total.

Try a fresh survey by adding `{"quantity":"infiltration","bounds_W_K":[32,48]}` to the `measured` array in a copy of the example. The choice becomes **Insulate fabric**. An interval of **[170,185]** favors **Seal unintended leaks**. Both follow-up readings are hypothetical. The file `surveyed.json` contains the first follow-up ready to run.

## Use your own survey and quotes

* `total_W_K` is the measured total of the modeled components. `heating_degree_hours_Kh` is annual heating degree-hours for your chosen base temperature. Multiply degree-days by 24.
* `measured` accepts any number of direct bounds on `fabric`, `infiltration` or `total`, in watts per kelvin. Each `tests` entry names the target, the survey's allowed absolute error in W/K, a unique ID, and the visit cost. A blower-door result at 50 Pa needs a site-appropriate conversion with uncertainty before use.
* Give two distinct upgrades with installed costs in the same currency and effective savings fractions in [0,1] on each whole component. Estimate those fractions for the specific quotes. Protect required ventilation and represent fixed or unaddressed losses explicitly by adjusting the inputs. A roof-only quote requires its incremental **whole-fabric** fraction, not the percentage improvement in the roof's U-value.
* `one_upgrade_dominates` identifies a strictly better quote across the full admissible region. `one_upgrade_never_worse` selects an upgrade that wins in some homes and can tie in others. `equal_value` indicates identical scores everywhere. `needs_measurement` means each upgrade wins somewhere in the feasible region. `model_mismatch` means the supplied bounds have no common solution. A missing next test means no available option guarantees a smaller decision interval or the decision already holds.
* Score is annual useful thermal kWh per unit of installation cost. Heat-system efficiency, tariffs, service life, moisture behavior, ventilation design and occupant comfort require separate inputs for a full works decision. If retrofit effectiveness or weather are uncertain, run distinct documented scenarios. The ranking is conditional on the entered values.

## Why the bounds cover the full range

Let F be fabric/bridge loss and A be modeled infiltration loss in W/K. The total measurement bounds F+A. Each survey clips the nonnegative (F,A) plane by two closed half-planes. Savings from a quote are linear in (F,A), so the extrema occur at polygon vertices, including point and line limits. The engine clips with exact rational arithmetic and converts final numbers for display.

For a proposed diagnostic measuring v·(F,A) with error e, two homes can share a reading exactly when their projections differ by at most 2e. The difference polygon consists of all differences between two feasible homes. Intersecting it with that strip and maximizing the score difference gives the worst remaining width after *any* reading. Tests are ranked by guaranteed width reduction per survey cost.

## Measured cross-check

The two DEEP homes underwent **room-in-roof insulation plus sealing together**. Table 3-5 gives baseline, post-work heat loss and reported uncertainty. The central changes are **42 W/K** in 07LT and **10 W/K** in 09LT; the study reports changes of **42 ±36** and **10 ±15 W/K**. Treating the separate before/after ranges as independent hard endpoints gives wider change ranges of **[−4,88]** and **[−10,30] W/K**. The code reproduces both central changes. These cases establish real measured total losses and the limits of attribution to individual upgrades. They do not validate the example quotes or predict counterfactual single-upgrade savings. Independent validation requires buildings with separate interventions, component measurements and predeclared tolerances. [DEEP Report 2.11, table 3-5 and section 3.3.1](https://assets.publishing.service.gov.uk/media/6717b9ace319b91ef09e3851/2.11_DEEP_07LT___09LT.pdf)

The public UK case measurements in `field.json` were transcribed from table 3-5, printed page 61, © Crown copyright, published under the Open Government Licence v3.0. The file records the source PDF hash. [Licence](https://www.nationalarchives.gov.uk/doc/open-government-licence/version/3/)

The command-line reader preserves decimal input text before exact comparisons. When integrating with an app that decodes numbers as floating point, pass precision-sensitive decimal values as JSON strings.
