# Battery capacity grader

**Does this used cell meet the capacity grade, or does it need a full test?** Enter voltages measured during an approved screening run. The grader keeps each reference cell compatible with every reading. It reports the lowest and highest reference capacities still possible after the stated allowances. A 4.80 Ah grade means a cell at exactly 4.80 Ah passes.

## Run a screen

Extract the kit, open this folder, and run these commands with Python 3.9 or newer. Use `py` in place of `python` on Windows if needed.

```sh
python app.py
python app.py example.json
python benchmark.py
python -m unittest -v
```

The example uses the first-minute voltage of a physical reference cell. Six profiles fit. Their capacity interval is **4.487–5.645 Ah**, spanning the 4.80 Ah cut. The next recommended measurement is at **300 seconds from the start**. At that checkpoint three pairs spanning the grade boundary have non-overlapping voltage bands. The candidate count is a count of reference traces and gives no probability for the incoming cell. The readout uses measured capacity in ampere-hours for the reference test protocol. Deliverable energy in watt-hours under a different load requires a separate labelled measurement.

A [grade-complete example](grade.json) and a [five-minute deferral](deferred.json) can be replayed with `python app.py grade.json` and `python app.py deferred.json`.

## Enter your readings

Copy [example.json](example.json). All voltage and time fields are absolute measurements on one consistent screening protocol. The five checkpoints are elapsed times of 60, 300, 600, 1200 and 1800 seconds.

| Field | Meaning |
| --- | --- |
| `min_capacity_ah` | Positive threshold for the grade. |
| `voltage_tolerance_v` | Absolute allowed difference between a query and each stored reference voltage. Include both instruments' error, rounding and the allowed mismatch between nominally similar cells. The example uses 0.020 V. |
| `capacity_error_ah` | Expand both ends of each reference capacity by this allowance. The example uses 0.050 Ah. |
| `observations` | Object mapping elapsed seconds to measured voltage. For example `{"60":3.207446,"300":3.393884}`. Empty is allowed for planning. |
| `policy` | Optional. `adaptive` chooses subsequent checkpoints. `fixed_300s` collects the 60 and 300 second readings, then defers unresolved cases immediately. |
| `reference` | Optional replacement library. Supply a list of objects with a unique `id`, measured `capacity_ah`, and `volts` mapping all five checkpoint strings to volts. |

The built-in references are 23 second-life LongTTech LFP cells with unknown degradation history. Each published label is the mean of three capacity tests. The reference charge trajectories were acquired at 25 °C with a 1C constant-current / constant-voltage protocol terminating at 3.65 V. These conditions identify the source data. Use qualified equipment, a compatible approved procedure and known starting conditions. Other chemistries, starting states, temperatures, and module sizes require their own reference library. Initial SOC varies in the archive. Only timestamps and terminal voltages enter the estimator. [Ji, Yue, Yan and Zhu, *Second-life LFP battery dataset*](https://doi.org/10.5281/zenodo.18630889)

| Result | Bench action |
| --- | --- |
| `meets_capacity_grade` or `below_capacity_grade` | All surviving capacity intervals agree on that grade **within the supplied reference model**. |
| `continue_screen` | Read the next elapsed-time checkpoint. Options rank the number of cross-grade reference pairs that the measurement must separate, divided by seconds spent waiting. If tied, prefer more separated pairs, then the earlier checkpoint. Even when two voltage bands overlap, a particular reading may separate them. |
| `reference_test` | Measure capacity directly. An empty survivor list indicates an unsupported reading. Remaining ambiguity at the last checkpoint or one surviving label straddling the threshold also calls for a reference test. |

The screening logic compares the decimal inputs including equality at band endpoints. Any guarantee depends on the actual cell lying within the specified allowances of a reference trace and capacity. An incomplete library can return a confident wrong grade. Safety, internal defects, thermal response, calendar life and vehicle serviceability require independent assessment.

## What happened on cells held out of the library?

The replay removes one entire physical cell, screens it against the other 22 and then reveals its published capacity. It repeats for all 23. Each fixed-time policy uses **all** available checkpoints up to its time limit. The adaptive policy starts at one minute and picks subsequent checkpoints from survivors. Every policy uses the same illustrative 4.80 Ah threshold and error allowances above.

| Screening policy | Correct early grades | False accepts | False rejects | Deferred | Total time on short screens |
| --- | ---: | ---: | ---: | ---: | ---: |
| Adaptive | 12 | 0 | 0 | 11 | 3 h 37 min |
| Fixed 5 minutes | 12 | 0 | 0 | 11 | 1 h 55 min |
| Fixed 30 minutes | 8 | 0 | 0 | 15 | 11 h 30 min |

**The five-minute screen graded the same twelve cells as the adaptive screen in less time.** In this archive, a held-out unit can be compatible with a reference early in the charge and incompatible later. That explains how the fixed 30-minute rule defers more units. Full-test time for deferred cells is additional and was not measured in this study. For these data the simpler operating policy is to collect 60 and 300 second readings together and send unresolved cells to a full capacity test. The per-cell outcomes are in [benchmark.json](benchmark.json). A larger prospective intake cohort with an agreed false-accept ceiling is needed to establish any shop-wide benefit.

## Data origin and rebuild

The small table in [data/measurements.csv](data/measurements.csv) contains actual timestamp-matched voltages rounded to six decimals, with the workbook filename and the original capacity label. The extractor checks the published file bytes by SHA-256 and takes each sample at its recorded elapsed second. The original 23 workbooks are distributed by Minghui Ji, Hong Yue, Yan Shao and Qunzhi Zhu under CC BY 4.0. This table is a reduced adaptation of their versioned dataset. [Original ZIP and metadata](https://zenodo.org/records/18630889) · [License](https://creativecommons.org/licenses/by/4.0/)

```sh
python prepare_data.py --out rebuilt.csv
```

If you already downloaded the upstream archive, use `python prepare_data.py SLBs_LFP_charging_data.zip --out rebuilt.csv`. The raw ZIP SHA-256 is `16afff956da104fcfb0acae5ee86a53f4c5c5e030c06f6fa8e7dfd566c715d07`. An additional EV-derived dataset covers six stationary duty cases with different pack protocols. It is a candidate for a separate library and external study. [Beatty et al., *Long-Term Sweat Testing Dataset for Second-Life Batteries*](https://doi.org/10.1038/s41597-025-05360-7)
