# Official hazard query evidence for the downtown study envelope

These responses were retrieved from FEMA and DOGAMI between **2026-11-24 00:22:22 and 00:24:53 UTC**, using the expanded envelope **[433750, 5033050, 434350, 5034000] in EPSG:26910**. This 600 by 950 metre box includes the initial downtown query and the tested approach-offset sweep area. It is an investigation area, not a proposed support, bridge, property or evacuation-route footprint. The repository supplies compact clipped geometry plus the complete query provenance; original statewide-feature responses are retained externally.

## Observed results

| Layer | Envelope-query count | Attributes actually returned | Interpretation |
|---|---:|---|---|
| FEMA NFHL 28: Flood Hazard Zones | 8 | Two `AE` features; two `X` with `0.2 PCT ANNUAL CHANCE FLOOD HAZARD`; four `X` with `AREA OF MINIMAL FLOOD HAZARD` | Different mapped zones intersect the box. Intersect actual proposed footprints with each polygon. No returned feature has a floodway subtype, but this query does not determine the legal consequences of an actual proposal. |
| DOGAMI OTGC1 12: Tsunami Evacuation Zones | 2 | OBJECTID 769: `Local Tsunami Evacuation Zone`, source `AllaJC2020`; OBJECTID 1199: `Outside Hazard Area` | The expanded box includes both classifications. A feature returned by this service is not automatically a hazard-positive classification. |
| DOGAMI OTGC1 78: XXL1 local Cascadia inundation | 1 | OBJECTID 1; statewide multipart geometry, with no local depth attribute | A full source feature was returned. Compute its actual intersection with each proposed footprint. This observation does not mean the complete study box is inundated or supply an inundation depth. |
| DOGAMI SLIDO42 3: mapped deposits | 0 | ID query returns `objectIds: null`; independent count query returns `count: 0` | Successful queries with no mapped deposit feature selected. This is not a finding of absent landslide hazard or stable ground. |
| DOGAMI detailed susceptibility 0: study coverage | 1 | OBJECTID 12: `REF_ID_COD = CalhNC2020a` | A study boundary intersects the box. This is coverage evidence, not a hazard classification or proof that every point in the box is covered. |
| DOGAMI detailed susceptibility 1: deep susceptibility | 1 | OBJECTID 2: `Suscept = Low` | This layer returns a low-susceptibility feature. Preserve its actual spatial extent; low mapped susceptibility is not measured bearing capacity or a site approval. |

FEMA feature 283667 reports `STATIC_BFE = 19.5`, `LEN_UNIT = Feet`, `V_DATUM = NAVD88`. This is a mapped base-flood elevation associated with that feature, not a flood depth, ground elevation or required bridge clearance. The other seven FEMA features use `-9999` in `STATIC_BFE`; treat that sentinel as unavailable. Other attributes also contain nulls or `-9999` sentinels.

The output geometry is EPSG:26910 XY. Existing source attributes such as `Shape_Area`, `Shape_Length` and `SHAPE.STArea()` retain their source meanings and units; **do not relabel these attributes as UTM square metres or metres after requesting transformed geometry**. Compute any local projected area from the returned geometry and record that separate derivation.

## Audit and completeness

- [manifest.json](manifest.json) records the study box, clipping operations, clipped files and hashes, source feature IDs, labels, geometry validity checks and limitations. The six `*.clipped.geojson` files use EPSG:26910, explicitly declared in each file. They are not RFC 7946 longitude/latitude GeoJSON.
- [raw-manifest.json](raw-manifest.json), the six layer summaries and [query-provenance](query-provenance) preserve the original endpoint, exact parameters, request URL, HTTP status, retrieval timestamps, byte count and SHA-256. The corresponding original response bodies are external, at the directory identified by `manifest.raw_storage.external_directory`.
- Each layer has metadata, count and ID queries. Nonempty results have full-geometry requests by bounded object-ID batches. Each successful batch's returned IDs match the requested IDs; no transfer-limit flag was reported. This uses ID batches rather than assuming a single default-size query returns all features.
- FEMA's initial eight-ID geometry request and one identical retry returned HTTP 503. Both failure responses and their metadata were retained. Two smaller four-ID requests subsequently returned all eight expected features. The resolved retrieval is complete; the failed observations remain part of the audit trail.
- For the zero-deposit query only, `objectIds:null` plus a successful independent count of zero is normalized to an empty ID list. Both source responses remain unchanged. A failed or inconsistent query would retain status `unknown` instead.
- [validation.json](validation.json) records verification of **26 response hashes** and clipping status. All six successful source geometry batches passed ID controls and requested EPSG:26910. The clipping manifest records GEOS validity and projected area for each nonempty clipped feature, as well as coordinate containment within the clip window.
- [retrieve_hazards.py](retrieve_hazards.py) captures full responses externally and can invoke [clip_hazards.py](clip_hazards.py). Clipping uses GDAL 3.10.3, its ESRIJSON reader, and GEOS through `ogr2ogr -clipsrc`. It does not manually rebuild rings, skip polygon organization, simplify shapes or apply topology repair. The writer uses 15 decimal places; the process uses floating-point GIS arithmetic.

The original uncompressed capture is about **40 MB** because several selected DOGAMI records contain statewide multipart geometries. Those originals are not included in this repository. The shipped geometry is approximately **76 KB**, produced after retrieving the full selected records. Raw hashes permit auditing against the retained originals. The source query did not request geometry simplification or rounding; clipping is an explicit subsequent operation.

To reproduce both capture and clipping from this directory, use a new external raw directory and a separate output directory:

```sh
python retrieve_hazards.py \
  --output-dir /absolute/path/to/new-hazard-raw-snapshot \
  --bbox 433750 5033050 434350 5034000 \
  --clip-output-dir /absolute/path/to/new-clipped-output
```

A rerun retrieves current source data and can legitimately differ. Preserve this snapshot before overwriting it. The earlier, smaller-box capture is separately retained at `/workspace/scratch/9f2debab16eb/tillamook-source/hazard-raw/2026-11-24-initial-envelope`; it is not the coverage basis for this version.

Count and ID queries are sequential observations, not an atomic source snapshot. Mapping can change. These layers do not collectively constitute a complete geotechnical, seismic, liquefaction, utility, flood-design or evacuation investigation. The [Oregon design basis](../../oregon-design-basis.md) describes their proper use and the remaining evidence gates.
