# Which material structures still fit the measurement?

**A powder scan can fit several crystal lattices. Even an exact list of distances between atoms can fit two different arrangements.** This application keeps every structure allowed by its stated model, and lets you check what another measurement would resolve. Two runnable modes cover a finite one-dimensional structure and the continuous cell-length parameter of a cubic powder pattern.

| Try it | Question |
| --- | --- |
| `python app.py` | Which six-site structures reproduce Bloom's fifteen pair distances? |
| `python app.py resolved.json` | Which one survives a local image showing one three-site group spanning four units or less? |
| `python app.py bounded.json` | What does the solver return when a noisy search runs out of work budget? |
| `python app.py powder.json` | Which cubic lattices can index six measured X-ray peaks from a real standard? |
| `python app.py powder_reference.json` | What changes if the specimen's independently published cell length is supplied? |
| `python -m unittest -v` | Replay the measurements, compare reconstruction with brute force and check the edge cases. |

Use Python 3.9+ and run these commands in this folder. On Windows you can use `py` in place of `python`. All examples run offline. Copy an example, edit the values and give its path to `app.py`.

## A measured powder example

The Advanced Photon Source published the GSAS-II peak-fitting tutorial's 11-BM run **3844** on lanthanum hexaboride, NIST SRM 660a. In the raw instrument file we selected the six strongest maxima between 5° and 14° in $2\theta$, requiring maxima to be more than 0.4° apart. The reported wavelength is 0.413128 Å. The selected angles are **5.6960°, 8.0595°, 9.8750°, 11.4075°, 12.7595° and 13.9830°**. `powder_source.py` downloads the pinned raw file or accepts a local copy, checks the recorded hash and regenerates the six fitted-input intervals. The original file includes 132,992 scan samples. [APS tutorial data](https://github.com/AdvancedPhotonSource/GSAS-II-tutorials/tree/494c2700cde7c92591e328fdce1515c662cf11ab/FitPeaks/data)

**The six peak positions admit three types of conventional cubic cell.** In the search range 3–9 Å, primitive indexing allows lengths near 4.157, 5.879 or 8.314 Å. Body-centered indexing allows the latter two ranges, and face-centered indexing allows the largest range. An external length measurement near 4.157 Å separates these possibilities. Applying the NIST certified interval at 22.5 °C selects the primitive lattice. The standard's identity and certificate are public; this replay demonstrates a consistency check rather than blinded identification of an unknown. A primitive-cell label alone does not determine the chemical species or atomic coordinates. [NIST SRM 660a certificate](https://tsapps.nist.gov/srmext/certificates/archives/660a.pdf)

**Scope of the measurement.** We assumed ±0.01° per peak for this exercise and held the wavelength fixed. Bragg's law gives $d=\lambda/(2\sin\theta)$; the stored $(d_{min},d_{max})$ pairs are rounded outwards at one millionth of an angstrom. This simple peak picker neither models peak shape nor supplies an instrument uncertainty budget. Overlapping peaks, multiple phases, preferred orientation, intensity and missed reflections can require a different analysis.

**For your powder data**, supply `"mode":"cubic_powder"`, `"cell_range_A":[low,high]` and a list of nonoverlapping `d_spacing_intervals_A`. Optionally restrict `centerings` to a list chosen from `P,I,F`. Use decimal strings when precise boundaries matter. In a cubic cell, $a^2=d^2(h^2+k^2+l^2)$. Body centering requires even $h+k+l$, face centering requires the three indices to have the same parity, and primitive centering permits all indices. For each line the software takes the union of the corresponding exact rational intervals in $a^2$. It intersects across all observed lines. Every possible index is covered up to the bound inferred from the maximum cell length and minimum observed $d$. The program accepts data requiring indices up to 32, and returns an input error for a broader search. Exact squared intervals are authoritative; square-root decimal displays are approximate. A permitted reflection can be physically weak or absent. Only the presence of the supplied peaks is used. At the same cell length, the primitive model permits every index admitted by body or face centering. With P allowed, positive peak positions alone therefore cannot uniquely select I or F. A larger conventional cell or a supercell can also reindex lines from the same physical phase. Distinguishing these indexings needs independently justified cell information or detectable missing-reflection and intensity constraints.

## The two six-site structures

For the fifteen distances 1 through 13, 16 and 17, the complete reconstruction gives the two lists below after accounting for shifting or reversing the entire sample. They have the same ideal diffraction intensity at every frequency for identical scatterers in this one-dimensional model, because each intensity depends on the multiset of pair displacements.

| Sites | Three-site groups spanning at most four units |
| --- | ---: |
| 0, 1, 4, 10, 12, 17 | 1 |
| 0, 1, 8, 11, 13, 17 | 0 |

The solver visits nine partial reconstructions in this example. Independent brute-force enumeration would consider 1,820 choices of four interior sites on a length-17 line. A proposed local observation must be physically resolvable for the actual material. [Bloom's original counterexample](https://doi.org/10.1016/0097-3165(77)90013-9) · [Pair distribution ambiguity in materials](https://doi.org/10.1039/D4FD00106K)

**For your distance data**, supply a complete multiset of nonnegative integer pair distances for 2–20 distinct equally weighted sites on an integer line. For exact distances the largest remaining distance is from an unplaced site to one of the two ends, giving two possible positions. Recursively subtract the newly explained distances and merge translated/reflected copies. For positive `distance_error`, every true distance must be matchable to a different measured distance within the integer error, allowing the total span to vary as well. The search enumerates sorted integer sites and prunes partial configurations that cannot be matched to the readings. These two procedures are exact for the stated model when the search completes.

Optional `probes` are objects with `sites`, `span`, and integer `cost`. Each asks for the number of distinct groups of that many sites fitting within that span. The default probes use 3–5 sites and span thresholds where a candidate prediction changes. The application chooses the smallest worst-case surviving group and breaks ties by cost, group sizes and probe dimensions. An empty list disables suggestions. Enter a known observation with `{"sites":3,"span":4,"count":1}` in `observations`.

`incompatible` means no model within the chosen bounds. `ambiguous` with no next probe means none of the offered probes separates the remaining models. `node_budget` defaults to 50,000 and can be set up to 250,000. If the search reaches that budget or the 1,024-candidate limit it returns `search_incomplete` and withholds uniqueness and probe recommendations. A conventional powder pattern has intensities, angular resolution, missing peaks and material-specific scattering factors; supplying its peak positions as a one-dimensional atom-distance multiset would misinterpret the data.
