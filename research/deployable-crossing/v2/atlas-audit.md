# Tillamook V5 source audit and crossing testbed

The recovered Atlas supplies a useful **mapped research testbed** for six downtown junctions. It supplies road centerlines, source feature identities, building and parcel outlines, and a terrain grid. It does not supply measured pavement widths, verified support locations, a vertical accuracy bound, or permission to place a crossing.

The reproducible adapter is [atlas_adapter.py](atlas_adapter.py), the extracted data is [data/atlas-subset.json](data/atlas-subset.json), and the detailed audit is [atlas-audit.json](atlas-audit.json). This audit was performed on 2026-11-24. Paths beginning `Content/` or `BuildInputs/` below refer to the original staging archive.

## Source identity and integrity

The source is `Tillamook-Atlas-V5-Staging.zip`, 35,341,779 bytes, with SHA-256:

```text
90f06d6f8bd4a67e51350f2ae849a1f16440925493c1ab2b7111241f1e1f79a8
```

Its 718 ZIP members contain 708 regular files. All 707 entries in `CONTENTS.sha256` match their file bytes; the manifest itself is the only regular file not listed in that manifest. The manifest hash is pinned in the adapter. The six OSM/parcel input hashes recorded in `Content/data/streets.json` also match the corresponding `BuildInputs/` files. The adapter rejects a different manifest, a mismatched file, or a different coordinate/grid contract before extraction.

The supplied subset is 747,966 bytes; file SHA-256:

```text
24b3a78c2e4e7cac3f55005e806d4476a39a2940cd03e27fef03ba7d6fbcc79b
```

The subset also carries a canonical content digest, source file hashes, and a separate terrain crop hash. These identify bytes and transformations; they do not certify physical conditions.

The GitHub repository [Aether-Ash-game-source](https://github.com/towre676-cloud/Aether-Ash-game-source) was checked at main commit `4c001ad283318be5b8eac84c153e1b9dce3e1836`. Complete recursive trees at all 27 branch tips (26 unique commits) contained no Tillamook-named files. Its [Site export provenance](https://github.com/towre676-cloud/Aether-Ash-game-source/blob/4c001ad283318be5b8eac84c153e1b9dce3e1836/docs/SITE-EXPORT.md) refers to Site version 45. Current Site metadata returned version 248 and source commit `527768d52282c29dd2d81fb2219e58ea8a680b39`, but source access and unauthenticated Atlas reads were unavailable. **The V5 staging package's production deployment has not been verified.** No game source or Site was changed.

## Coordinates and terrain

`Content/data/town.json` declares NAD83 / UTM zone 10N, **EPSG:26910**, with origin `[430000, 5028500]` metres. Local `x` is easting minus 430000 and local northing is northing minus 5028500. The game's horizontal `z` coordinate is therefore northing, not elevation. The context bounds are `[430000,5028500,437500,5036500]`, covering 7.5 × 8 km.

The adapter selects the local window `[3750,4550,4350,5450]`, a 600 × 900 m downtown window giving every candidate at least 250 m of map context. Vector features are selected by bounding-box intersection and retained whole. They can extend outside the window. No clipping or simplification changes their vertices, ring order, holes, or source IDs.

`compile-streets.py` projects the original OSM longitude/latitude coordinates and rounds local vertices to 0.01 m. Streaming fragments add a third decimal through interpolation. The adapter uses the source paths rather than those display fragments. **Centimetre encoding is not centimetre positional accuracy.** The build rubric says external withheld controls are still needed. No usable horizontal survey error bound is supplied.

The terrain crop comes directly from `Content/data/survey-ground.bin`, without the game's 120 m blend into the older 5 m DEM. Its source grid is 1750 × 1750 cells covering local `[2500,4500,6000,8000]`, with 2 m spacing. Signed little-endian Int16 values represent metres at scale 0.025. Rows run north to south and columns west to east; samples are at cell centres. The crop is columns 625–924 and rows 1275–1724, giving 300 × 450 samples. Its outer cell edges match the subset window.

The V2 recipe `Content/refs/repertoire/tillamook.json` asserts NAVD88 and converts source feet using 0.3048. `Content/refs/Sources.md` explicitly requires source-lineage confirmation before an independent elevation claim. Original float DEM/DSM inputs are not included in V5, although their hashes and compiler recipe are retained. The recorded DEM comparison percentiles are an internal consistency comparison, not an independent datum or error validation. The subset therefore leaves vertical datum verification and accuracy unresolved. The 0.0125 m half quantization step covers rounding alone; it is **not** a bound on physical ground elevation.

`sample_terrain` returns bilinear interpolation in metres. `sample_terrain_exact` returns the same interpolation as an exact `Fraction` relative to the encoded samples. Both reject points outside the crop; interpolation clamps to the nearest available sample within the outer half cell. Neither function turns a 2 m raster into a curb, bearing-pad, utility, or clearance survey.

## Six source-bound junctions

Each candidate is an actual shared OSM node in the compiled road paths. Latitude and longitude come from the original node record. The normalized road tangent is the secant through the immediate neighbouring road vertices, oriented toward positive northing. It is a reproducible map direction, not a surveyed road alignment.

| Candidate | Local east, north (m) | OSM node | Road spine | Cross-street ways |
|---|---|---|---|---|
| Main at 7th | 4008.95, 5161.69 | `n36984060` | `w160662295` | `w1557734242` |
| Pacific at 7th | 4072.39, 5160.73 | `n36985866` | `w160662310` | `w1557734242` |
| Main at 9th | 4007.75, 4996.25 | `n36959806` | `w160662295` | `w671168314` |
| Pacific at 9th | 4070.61, 4995.41 | `n36959805` | `w160662310` | `w671168314`, `w5246810` |
| Main at 11th | 4006.85, 4832.06 | `n36984056` | `w160662295` | `w5250328`, `w1557734251` |
| Pacific at 11th | 4068.53, 4831.62 | `n36985862` | `w160662310` | `w5250328` |

The raw OSM tags identify Main and Pacific as separate one-way US 101 ways, with two lanes, a 25 mph tagged limit, and `sidewalk:both=yes`. These are snapshot attributes, not a new traffic survey. The cross streets also carry `sidewalk:both=yes`.

All 244 selected path widths are authored. Main and Pacific display at an assumed 11 m; local street and footway widths likewise come from renderer defaults. Across the entire Atlas, 2,042 of 2,044 path widths are authored; the two OSM-attributed widths are runways. The adapter preserves `widthEvidence` and never promotes a display width to a measured span.

No separately mapped footway path lies within 25 m of any of the six junction centres. That is a mapping gap despite the sidewalk tags, not evidence that sidewalks are absent. No crossing mark is mapped within 25 m at the 7th or 9th junctions. At 11th there are three crossing nodes near Main (`n9982988085`, `n9982988086`, `n9982988087`) and two near Pacific (`n9982988083`, `n9982988084`). Their positions and original crossing tags are retained. Missing marks do not establish the absence of a legal or physical crossing.

## Retained evidence and limits

| Layer | Retained count | Evidence and limitation |
|---|---:|---|
| Streets and other linear paths | 244 | Original source paths; widths remain authored |
| Crossing, stop, signal and other marks | 58 | Original mapped points and relevant raw OSM tags |
| OSM areas | 17 | Original way outlines; not a complete land-use inventory |
| DOGAMI building footprints | 395 | Original IDs and rings; no authored façades, roof meshes or estimated building heights |
| Parcel parts | 401 | All rings retained, including 78 interior rings; outlines do not establish ownership or placement permission |
| OSM environment features | 40 | 13 points, 5 lines and 22 areas; some IDs also occur in other layers |
| Canopy proxies | 311 | Derived locations and estimated heights/radii; not surveyed individual tree envelopes |
| Terrain cells | 135,000 | Raw unblended encoded grid; unresolved vertical accuracy |

Environment records use `geometry` equal to `point`, `line`, or `area`. Their `p` field is respectively one `[east,north]` pair or a list of pairs. Buildings and parcels use `rings`; the first ring is the exterior and subsequent rings are retained holes. The selected buildings happen to have no holes; the extraction does not discard holes. OSM environment catalogue completeness, underground utilities, overhead wires, access rights, support bearing, and obstacle clearance remain separate questions.

The compact subset excludes imagery, executables, rendering libraries, authored façades, game rooms and campaign state. OSM-derived data retains **© OpenStreetMap contributors** and the [Open Database License 1.0](https://opendatacommons.org/licenses/odbl/1-0/), with [OSM's attribution guidance](https://www.openstreetmap.org/copyright). Keep that database attribution and license with redistribution. DOGAMI, ODOT and terrain derivatives retain their distinct source notices; this audit does not infer a blanket redistribution license for every third-party layer or grant new rights to them.

## Reproduce

From this directory, using the extracted V5 archive root:

```sh
python -S atlas_adapter.py /absolute/path/to/atlas-v5 --output /tmp/atlas-subset-rebuilt.json
cmp data/atlas-subset.json /tmp/atlas-subset-rebuilt.json
python -S atlas_adapter.py --check
```

The source directory is read only. The extraction is deterministic and uses only Python's standard library. Independent audit checks compared 104 selected cell centres and 22 bilinear midpoint values against the original source binary, rejected 11 damaged or out-of-bounds inputs, confirmed unique IDs within each category, and rebuilt byte-identically. These checks validate extraction and representation. Engineering site suitability is outside their scope.

### Optional independent site geometry replay

[verify_site_independent.py](verify_site_independent.py) uses installed GDAL/GEOS rather than the study's geometry helpers. Run `/usr/bin/python3 verify_site_independent.py` where the `osgeo` bindings are available. This is an optional check, not a dependency of the standard-library adapter or study.

The saved [independent receipt](verification-site-independent.json) binds the input subset and study report hashes and records GDAL 3.10.3 / GEOS 3.13.1. All 5,184 hit sets across the 648 scenarios agreed, including polygons with holes, line-distance sensitivities, points, and canopy proxies. All tested building, parcel, environment and scenario geometries were valid. Independent projection of the six original OSM node coordinates reproduced their stored local coordinates after centimetre rounding; those are internal transform checks, not withheld survey controls.

The independent check confirmed 636 layouts with a tested building or vehicle-centerline conflict and 12 without either tested conflict. Those 12 layouts still have a mapped building only 0.138805–2.138805 m away; all are within the separate 3 m building sensitivity. This finding provides no approved placement or verified clearance.
