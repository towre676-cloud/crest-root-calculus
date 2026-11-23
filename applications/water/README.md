# Where should a leak crew listen next?

Start with pressure drops, a list of possible leak locations, and a bounded response for each available logger. The workbench returns the full set of compatible locations. It identifies the next logger that guarantees the greatest reduction in worst-case search length per unit of visit cost. It accepts both a town's own scenarios and the packaged field replays.

Run from this folder with Python 3.9 or later.

```sh
python app.py
python app.py example.json
python app.py logged.json
python app.py field-case.json
python field_replay.py
python field_replay.py --event 1 > my-event.json
python app.py my-event.json
python -m unittest -v test_app.py
```

## A first dispatch

In the two-street example the inlet records a pressure drop between 1.9 and 2.1 kPa. A 180 m reach on Orchard and a 120 m reach on Mill both fit; the 75 m Depot reach drops out. Both available loggers reduce the worst remaining search from 300 m to 180 m. North takes half the visit cost of South. Add `"North":[2.95,3.05]` to the observations and run again. Only Orchard remains. The file `logged.json` contains this illustrative follow-up. The place names, lengths and pressure responses in this teaching example are constructed.

## Enter a case

* `zones` assigns a positive search weight to each disjoint, non-overlapping unit crews can inspect. Set `extent_unit` to `m` for pipe metres or `sites` for counts.
* `sensors` defines available loggers, their maximum absolute reading error, and positive installation/visit cost. State the pressure unit for the case. An observation has the form `[measured−error, measured+error]`. Error is already included in these input bounds.
* Each candidate scenario has one `zone` and a nonnegative interval for a common leak strength q. The predicted drop at sensor s is `slope[s] × q + offset[s]` where the offset is a closed interval for local demand, model variation or baseline drift. Omitted offsets equal zero. Negative and zero slopes are supported. Supply locally derived responses from a hydraulic model or calibrated traces. Separate scenarios can describe differing demand regimes, including several scenarios per zone.
* The solver inverts each observed interval for q and intersects those intervals. The same strength must explain every reading. For an unused logger, it predicts every surviving possible measurement including sensor error. Sweeping the closed interval endpoints computes exactly the largest sum of *distinct* compatible zone weights for a possible next reading. Touching intervals remain ambiguous. Ranking divides the guaranteed improvement by visit cost.
* An empty set is a model mismatch. An absent next logger means the listed options cannot guarantee a smaller worst-case search. The conditional coverage holds for the listed scenarios and independent bounded offsets. Multiple simultaneous leaks, missed candidate locations, correlated demand and unmodeled network conditions need corresponding new scenarios or a more detailed forward model.

## Field replay on measured pressure

Mazaev and colleagues recorded deliberately released water at known points in operating networks. The included December BKTown extract uses their pressure logs, event times, released volumes and anonymized uncalibrated EPANET network. The first and second logger are selected with the location hidden. The known discharge volume and duration set the same flow hypothesis at **every** candidate. This is a known-rate test over nine published candidate sites. A real dispatch also needs uncertainty in leak rate and a broader location library.

Of nine releases, eight had a clean pre-event window. The last began as the preceding release ended and is explicitly excluded. For every usable logger we take the median pressure in minutes −8 through −3 before the event and subtract the median in minutes +10 through +29. The first complete logger supplies a common reference to reduce systemwide drift. Each differential response includes the range of same-time, same-window differences on the seven reference days December 9–15, plus an illustrative ±0.025 bar measurement error. All retained sensors have complete windows. Labels map by row order of the source logger list. Eight *events* were scored, not eight independent networks.

The hydraulic calculation predicts each site's pressure step separately by running a zero-duration snapshot without and with the known extra demand at the candidate. It uses the source demand and head patterns at the event midpoint, the original initial link states and EPANET 2.2. It does not warm-start a time-dependent state history. The provided network generated EPANET warning 6 for disconnected nodes. All reported logger reference pressures were positive and within the checked range (0,16) bar. The observed discrepancies are part of the reported result.

| Method over eight usable releases | True site retained | Total sites returned | Empty sets |
| --- | ---: | ---: | ---: |
| First logger, in addition to the reference | 5/8 | 33 | 1 |
| First + adaptive logger | 3/8 | 20 | 2 |
| First + fixed second logger | 4/8 | 24 | 1 |

The fixed baseline reads the next logger in source order. A no-improvement or singleton adaptive case uses that same fallback to keep the two-reading comparison fair. A single nearest-model answer with the selected sensors was correct **2/8**. Smaller lists hid real leaks when the pressure model was too narrow.

A second replay learns the smallest single extra model allowance that covers every available true-site logger on the four December 16 releases, then applies it to the four December 17 releases. The fitted allowance is **0.238617 bar**. Both adaptive and fixed two-logger choices retain the true site **4/4** and return **28 total sites** of the possible 36. These data show the value of allowing model discrepancy and provide no advantage for adaptive selection at the same retention on this small later-day test. The ordinary input `field-case.json` replays the first December 17 release at 07:40, using the previous day’s guard, the reference logger and two measured differential readings. Seven of nine zones survive, including Z05 where the release occurred. Before relying on a smaller search region in the field, use independent repairs in the intended network to verify site retention at a matched travel and inspection budget.

## Provenance and regeneration

The compact `field.json` contains independently simulated pressure drops, relative background ranges, measured pressure drops and event labels. The raw network map and full monitoring history are unnecessary for normal use. The original **Mazaev, Weyns, Vancoillie, Vaes, Ongenae and Van Hoecke** data are distributed under CC BY 4.0. They have been transformed as described above, including renamed site/logger labels, rounding of pressure values to six decimals and our new hydraulic simulations. [Source data and licence](https://doi.org/10.5281/zenodo.7255403) · [Authors' baseline implementation](https://github.com/predict-idlab/phys-ml-leak-localization) · [CC BY 4.0](https://creativecommons.org/licenses/by/4.0/)

To regenerate, download `data.zip` from that exact DOI and locate the native EPANET 2.2 library for your platform. The reference run used the Linux x64 `libepanet22.so` packaged with WNTR 1.5.0. The importer checks the source ZIP SHA-256 before opening it and reads CSV/text members only. Run `python build_field.py data.zip /absolute/path/to/libepanet22.so --output new-field.json`. The normal replay reads the shipped extract using Python's standard library.

The command-line reader preserves decimal input text before exact comparisons. When integrating with an app that decodes numbers as floating point, pass precision-sensitive decimal values as JSON strings.
