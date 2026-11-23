# Choose a telescope follow-up window

**Which observing slot would tell two possible orbits apart?** Give the scheduler periods that fit the evidence, ranges for their timing errors, and windows when the telescope can observe the target continuously. It finds a slot where catching a full transit would select one explanation and a sufficiently sensitive, uninterrupted observation with no full transit would select another.

## Run the historic example

```sh
python app.py
python app.py confirmed.json
python benchmark.py
python -m unittest -v
```

Run inside this folder with Python 3.9+. On Windows, `py` also works. [The editable example](example.json) contains a 55 Cancri e replay. Both proposed orbital periods are anchored at the same measured Spitzer transit. During the first illustrative 3.6-hour slot both models predict no complete transit. A 3.6-hour slot around the short-period return predicts a complete transit for the 0.736539-day orbit and none for the 2.81705-day alias. The longer-period return is another distinguishing slot. With equal cost and separating power the earlier slot wins. These are proposed windows on the historical time axis. Actual telescope availability at those windows is an assumption in this example.

[The confirmation input](confirmed.json) adds the independently reported MOST photometric transit epoch. Under the stated timing allowances the short-period model survives, and the daily alias is excluded. Both periods were given the same Spitzer phase anchor. This reconstructs a test between two published period hypotheses. The predicted visits were not actually carried out for this replay. Winn et al. obtained a folded-light-curve epoch from multiple transits; this single reported timing is the confirmation datum. The replay does not claim to recover the two periods from raw spectra.

## Make a new observing plan

Copy [example.json](example.json) and edit the following fields. Use quoted decimals for dates when another program may round a long numeric literal. All dates, periods and durations must be in the units specified below.

| Field | Meaning |
| --- | --- |
| `time_system` | Declare a common time standard for all timestamps, e.g. BJD-TDB. Convert epochs and the station's available windows to this standard upstream. |
| `hypotheses` | Each named model has `epoch_days` and `period_days` as closed lower/upper pairs. `duration_hours` bounds the complete first-to-fourth-contact transit duration. Give the same number twice to fix a parameter. |
| `events` | Each pair is the allowed **midpoint** interval of an observed transit from this same body. Empty is allowed. |
| `windows` | Each choice supplies a unique `name`, `days` as a start/end pair, and `available` as true or false. Supply contiguous usable coverage for the target; daylight, airmass, weather, overhead and instrument sensitivity must already have been assessed. Split visits with gaps. |
| `non_detections` | Optional list of uninterrupted start/end pairs in days. Use these only if the instrument would have detected a full transit of the depth and duration under consideration. A partially observed transit cannot count as a confirmed full event in this model. |

**Read the results.** A window gets one point for each pair of surviving hypotheses with disjoint possible yes/no outcomes. Rank points per observing hour, then total distinguished pairs, earlier start and name. `no_guaranteed_separation` means the offered slots cannot guarantee different outcomes for any pair. `model_mismatch` means no listed hypothesis fits. If you obtain a midpoint, append its measurement interval to `events` and rerun. If a qualified visit saw no complete transit, append its interval to `non_detections`. A surviving model is one of the supplied alternatives. An unlisted orbit may also fit.

## How the calculation works

At cycle number (n), the predicted transit midpoint is the reference epoch plus (n) times the period. One observed midpoint cuts a strip from the allowed rectangle in epoch/period space. For each candidate integer cycle the solver clips a convex polygon by that strip, keeping all compatible cycle histories across all measured midpoints. It uses rational arithmetic, inclusive boundaries and signed cycle numbers.

For a future window to contain a complete transit, the midpoint must be at least half a duration from either edge. The shortest allowed duration tests whether a detection is possible. The longest tests whether a missed complete event is possible, by looking for a gap between consecutive transit midpoints. The code checks linear inequalities on the surviving epoch/period polygons. Negative evidence removes each entire cycle history when a full transit was compulsory in that history. For a history partially consistent with an earlier non-detection the code retains the whole polygon. Future guarantees are therefore conservative in that case. Inputs with more than 512 cycle-history pieces or a single observation spanning more than 10,000 candidate cycles report a search-limit error. The ephemeris assumes a constant period and transit duration within its stated bounds. Stellar variability, transit timing variations and insufficient photometric depth sensitivity require additional modeling.

## Historical measurements and baseline

Dawson and Fabrycky's Tables 6 and 7 give periods 2.81705(5) and 0.736539(3) days. The demo uses three times each quoted last-digit fit error as an illustrative interval. Demory et al. measured the Spitzer transit near BJD 2455568.0265 (published +0.0015/−0.0010 d). We widen this to [2455568.022, 2455568.031]. Winn et al.'s Table 1 gives HJD 2455607.05562 ±0.00087 and a duration of 0.0658 ±0.0013 d. We use a symmetric ±0.00661 d allowance for the later timing, covering three quoted errors plus an illustrative 0.004 d timing-convention allowance. This broad window permits the historical comparison of reported BJD and HJD without claiming a precise time-standard conversion. The 1.40–1.75 h duration bracket is shared by both historical models and informed by the published transit; it is a retrospective design assumption. These are explicit model bounds, not confidence statements. [Dawson & Fabrycky](https://doi.org/10.1088/0004-637X/722/1/937) · [Demory et al., Spitzer measurements](https://doi.org/10.1051/0004-6361/201117178) · [Winn et al., MOST measurements](https://doi.org/10.1088/2041-8205/737/1/L18)

The benchmark checks the frequency difference between the candidate periods against **320 actual observation times** from Fischer et al. For each time, it forms a unit vector with angle (2\pi t(1/P_1-1/P_2)). The magnitude of their mean is **0.859588**, where 1 would mean perfect alignment at this frequency difference. The same observing cadence can therefore make the two periods difficult to distinguish. This samples the observing window only. The result uses the archived timestamps and does not fit the radial velocities or subtract signals from the four other planets. As a planning baseline, taking the first available slot in the worked case distinguishes zero pairs; the recommended slot distinguishes one for the same duration. The subsequent MOST epoch independently checks the surviving period in this example. [Recorded benchmark](benchmark.json) · [Fischer et al. observations](https://doi.org/10.1086/525512)

The NASA Exoplanet Archive source tables are included unmodified, with their source and instrument headers. Benchmark code checks SHA-256 before reading the 70 Keck and 250 Lick rows. To reproduce from the originals, save the two linked tables in `data/` under their published names and run `python benchmark.py`. [NASA archive, Keck](https://exoplanetarchive.ipac.caltech.edu/data/ExoData/0043/0043587/data/UID_0043587_RVC_001.tbl) · [NASA archive, Lick](https://exoplanetarchive.ipac.caltech.edu/data/ExoData/0043/0043587/data/UID_0043587_RVC_002.tbl)
