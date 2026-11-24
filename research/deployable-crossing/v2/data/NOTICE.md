# Data attribution and evidence scope

## Saved Tillamook Atlas V5

The source is the user's existing `Tillamook-Atlas-V5-Staging.zip`, pinned by SHA-256 `90f06d6f8bd4a67e51350f2ae849a1f16440925493c1ab2b7111241f1e1f79a8`. `atlas-subset.json` records the consumed file hashes and original source IDs. `atlas_adapter.py` verifies the source manifest before extraction.

- **OpenStreetMap:** © OpenStreetMap contributors. OSM extracts and derived OSM database content are redistributed under the [Open Database License 1.0](https://opendatacommons.org/licenses/odbl/1-0/). [Copyright and attribution](https://www.openstreetmap.org/copyright). This includes OSM streets, marks, areas and mapped environmental objects; the adapter retains relevant original tags. Road widths labelled `authored` are display defaults, not measured roads.
- **Oregon DOGAMI:** Building footprint source O-23-01, with original inventory IDs retained. [Source release](https://pubs.oregon.gov/dogami/ofr/O-23-01/p-O-23-01.htm). The selected building geometry excludes authored façade, room and exterior-appearance refinements. Terrain/canopy source evidence comes through the saved V2/V5 Atlas pipeline; source lineage and independent vertical accuracy remain unresolved, as recorded per attribute.
- **ODOT Tillamook County Taxlots:** Original parcel-part IDs and every ring are retained. [Source service](https://gis.odot.state.or.us/arcgis1006/rest/services/ames/ames/MapServer/39). Parcel boundaries are a captured map layer, not a legal boundary survey, ownership determination or access permission. No ownership/address attributes are added.
- **Other Atlas sources:** The original source recipe identifies USGS 3DEP, Oregon OSIP2024 imagery, OWRD/ODOT road centerlines and Census TIGER water. This compact subset does not redistribute the original aerial textures or full-town game. Its terrain crop uses the unblended DOGAMI-derived survey grid, not the game's mixed-resolution blending.

The Atlas data and external GIS layers do not acquire a new blanket software license through their inclusion here. Retain their attribution, source conditions and OSM's database license when redistributing. The source audit does not independently establish rights beyond the recorded source terms.

## Added official hazard layers

FEMA NFHL and Oregon DOGAMI layers retain original attributes and object IDs. Exact services, query parameters, UTC timestamps, source response hashes and GDAL/GEOS clipping steps are in [hazards/README.md](hazards/README.md), `hazards/raw-manifest.json`, `hazards/manifest.json` and `hazards/query-provenance/`.

The GeoJSON explicitly uses projected EPSG:26910 coordinates; it is not longitude/latitude RFC7946 output. Zero returned features, no geometric intersection, a mapped “Outside Hazard Area,” a susceptibility category and unavailable evidence have different meanings. None is a ground-capacity or hazard-free certificate.
