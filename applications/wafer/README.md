# Chip inspection

**Can these optical readings settle a line-width specification?** Enter the largest acceptable bottom width and the light intensities recorded by the instrument. The program keeps every NIST profile that explains all readings within your stated error. It reports whether every survivor gives the same width decision. If the survivors disagree, it ranks the available optical settings by the number of disagreeing pairs they can separate per unit cost.

## Try the three cases

| Case | Command in this folder | What happens |
|---|---|---|
| One reading, simulated example | `python app.py example.json` | 44 profiles meet the illustrative limit and 37 exceed it. Switching polarization separates 1,560 of the 1,628 disagreeing pairs. |
| Both polarizations | `python app.py two-readings.json` | Three geometries still fit. All three bottom widths fall between 126.5 and 134 nm, below the 135 nm limit. |
| Actual measured die | `python app.py measured-die.json` | At ±0.001, all 2,566 profiles disagree with one or more readings. The best profile's largest channel discrepancy is 0.0140153. |

## Enter another job

Copy a case and provide `width_limit_nm` (inclusive), `error` (absolute normalized intensity) and `observations` (an object mapping NIST channel numbers 1–84 to normalized intensities). Alternatively provide `die` from 1–9 and `observed_channels` (a list or `"all"`) to read actual NIST measurements. `available_channels` can restrict the next available settings. Optional positive `channel_costs` use any consistent unit. An empty observations object represents the whole library before measurement. Input decimals support 13 places. The example error and width limit are illustrative choices.

Two possible profiles can share an observation on a channel precisely when their predicted intensities differ by at most twice the specified error. These closed intervals are compared using scaled integers. For each allowed channel we count all surviving opposite-grade pairs whose intervals are disjoint. Unit costs are the default. Pair counts describe a finite design library. The distribution of shapes on a manufacturing line and the calibration errors of an actual instrument require separate measurements.

## Replay the comparison

Run `python replay.py`. The exact sweep compares a fixed normal-incidence setting, adding an adjacent angle, adding the orthogonal polarization, and all 84 settings. Of 2,566 grid profiles, the respective numbers whose grade is robust to every allowed reading are 927, 951, 1,649 and 2,533. The replay also compares all nine measured dies directly to the table. None fits within the example error across all channels. Resolving at manufacturing scale requires an independently measured reference width, model discrepancy estimated on held-out wafers, and comparison at matched false-accept rates. This application makes a decision only for the stated library and error.

The three original CSV files ship gzip-compressed. Both programs check their decompressed SHA-256 digests. See the original [NIST readme](data/NIST-readme.txt) for geometry and channels, and the [NIST measurement and simulation dataset](https://doi.org/10.18434/mds2-2290) for the full resource. Python 3.11+ needs no third-party packages. Run `python -m unittest -v` in this folder.
