#!/usr/bin/python3
"""Overlay recorded site footprints on attributed official GIS snapshots.

Generate with /usr/bin/python3 (optional osgeo dependency). Validate the saved
offline report with any Python 3 using --check; that path is standard library
only. Spatial intersections preserve feature attributes and never yield a
hazard-free, structurally adequate, or approved-site classification.
"""
from __future__ import annotations

import argparse
from collections import Counter, defaultdict
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).parent
FORMAT = "crest-root-tillamook-hazard-screen/2"
LAYER_ROLES = {
    "fema_nfhl_28": "FEMA flood-zone classification; retain zone, subtype, units and sentinel attributes",
    "dogami_tsunami_12": "Evacuation-zone category; Outside Hazard Area is a returned category, not a safety finding",
    "dogami_tsunami_78": "XXL1 modeled local-source tsunami inundation extent; distinct from the evacuation-zone layer",
    "dogami_landslide_3": "Mapped landslide-deposit inventory; no returned deposit does not establish absence",
    "dogami_ls_coverage_0": "Reference-map coverage only; this is not a landslide-susceptibility class",
    "dogami_ls_deep_1": "Deep-landslide susceptibility category; Low is not no hazard or a ground-capacity assessment",
}


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _canonical(data):
    return json.dumps(data, sort_keys=True, ensure_ascii=False, separators=(",", ":"), allow_nan=False).encode("utf-8")


def _digest(report):
    return _sha(_canonical({k: v for k, v in report.items() if k != "content_sha256"}))


def resolve_feature(report, intersection):
    """Resolve one normalized intersection to its complete source record."""
    try:
        return report["feature_catalog"][intersection["feature_id"]]
    except (KeyError, TypeError) as exc:
        raise ValueError("intersection refers to an unknown feature") from exc


def _covered(ring, bounds):
    return bool(ring) and all(len(p) == 2 and all(math.isfinite(x) for x in p)
                              and bounds[0] <= p[0] <= bounds[2] and bounds[1] <= p[1] <= bounds[3]
                              for p in ring)


def _absolute_footprints(scenario, origin):
    return [{"id": f["id"], "kind": f["kind"],
             "ring": [[p[0] + origin[0], p[1] + origin[1]] for p in f["ring"]]}
            for f in scenario["footprints"]]


def _label(key, properties):
    if key == "fema_nfhl_28":
        return " | ".join(str(properties[k]) for k in ("FLD_ZONE", "ZONE_SUBTY") if properties.get(k) is not None) or "unlabelled flood feature"
    if key == "dogami_tsunami_12":
        return str(properties.get("Tsunami_Region", "unlabelled evacuation feature"))
    if key == "dogami_tsunami_78":
        return "XXL1 modeled inundation extent"
    if key == "dogami_ls_coverage_0":
        return str(properties.get("REF_ID_COD", "unlabelled reference coverage"))
    if key == "dogami_ls_deep_1":
        return str(properties.get("Suscept", "unlabelled susceptibility feature"))
    return "mapped landslide deposit"


def _inputs(atlas_path, site_path, manifest_path):
    paths = {"atlas_subset": Path(atlas_path), "site_report": Path(site_path), "hazard_manifest": Path(manifest_path)}
    blobs = {key: path.read_bytes() for key, path in paths.items()}
    data = {key: json.loads(value) for key, value in blobs.items()}
    atlas, site, manifest = (data[key] for key in ("atlas_subset", "site_report", "hazard_manifest"))
    if atlas["crs"] != "EPSG:26910" or atlas["origin"] != [430000, 5028500]:
        raise ValueError("expected the pinned Atlas local EPSG:26910 frame")
    if site["atlas_content_sha256"] != atlas["content_sha256"]:
        raise ValueError("site report and Atlas identities differ")
    if manifest["schema_version"] != 2 or manifest["study_envelope"]["srid"] != 26910:
        raise ValueError("expected hazard manifest schema 2 in EPSG:26910")
    if {layer["key"] for layer in manifest["layers"]} != set(LAYER_ROLES):
        raise ValueError("expected all six declared hazard/reference layers")
    folder = paths["hazard_manifest"].parent
    layers = []
    hashes = {key: _sha(value) for key, value in blobs.items()}
    for entry in manifest["layers"]:
        relative = entry["clipped_file"]
        path = (folder / relative).resolve()
        if not path.is_relative_to(folder.resolve()):
            raise ValueError("hazard path lies outside its data directory")
        raw = path.read_bytes()
        if _sha(raw) != entry["clipped_sha256"]:
            raise ValueError(f"hazard clip hash mismatch: {relative}")
        collection = json.loads(raw)
        crs = collection.get("crs", {}).get("properties", {}).get("name")
        if crs not in ("urn:ogc:def:crs:EPSG::26910", "EPSG:26910"):
            raise ValueError(f"wrong hazard clip CRS: {relative}")
        if len(collection["features"]) != entry["clipped_feature_count"]:
            raise ValueError(f"hazard clip count mismatch: {relative}")
        summary_path = folder / entry["raw_summary"]
        summary_bytes = summary_path.read_bytes()
        summary = json.loads(summary_bytes)
        if summary["query_srid"] != 26910:
            raise ValueError("source query uses an unsupported CRS")
        bounds = summary["query_envelope"]
        clip = manifest["processing"]["clip_bounds"]
        available_bounds = [max(bounds[0], clip[0]), max(bounds[1], clip[1]),
                            min(bounds[2], clip[2]), min(bounds[3], clip[3])]
        if available_bounds[0] >= available_bounds[2] or available_bounds[1] >= available_bounds[3]:
            raise ValueError("source query and clip bounds do not overlap")
        known = entry["clip_status"] == "complete" and entry["query_status"] in {"features_returned", "no_features_returned"}
        hashes[relative] = _sha(raw)
        hashes[entry["raw_summary"]] = _sha(summary_bytes)
        paths[relative], blobs[relative] = path, raw
        paths[entry["raw_summary"]], blobs[entry["raw_summary"]] = summary_path, summary_bytes
        layers.append({"entry": entry, "features": collection["features"], "available_bounds": available_bounds,
                       "source_complete": known, "summary": summary})
    return atlas, site, manifest, layers, hashes, paths, blobs


def run_screen(atlas_path=ROOT / "data/atlas-subset.json", site_path=ROOT / "site-report.json",
               manifest_path=ROOT / "data/hazards/manifest.json"):
    try:
        from osgeo import gdal, ogr
    except ImportError as exc:
        raise SystemExit("Generation requires optional GDAL/GEOS: use /usr/bin/python3 hazard_screen.py; --check uses only stdlib") from exc
    gdal.UseExceptions()
    ogr.UseExceptions()
    atlas, site, manifest, layers, hashes, paths, blobs = _inputs(atlas_path, site_path, manifest_path)
    origin = atlas["origin"]
    code_sha = _sha(Path(__file__).read_bytes())

    def polygon(ring):
        boundary = ogr.Geometry(ogr.wkbLinearRing)
        for x, y in ring:
            boundary.AddPoint_2D(x, y)
        boundary.CloseRings()
        geometry = ogr.Geometry(ogr.wkbPolygon)
        geometry.AddGeometry(boundary)
        if not geometry.IsValid():
            raise ValueError("invalid scenario footprint")
        return geometry

    catalog, layer_metadata = {}, []
    for layer in layers:
        entry, key = layer["entry"], layer["entry"]["key"]
        features, seen = [], set()
        for feature in layer["features"]:
            if feature["geometry"]["type"] not in {"Polygon", "MultiPolygon"}:
                raise ValueError("hazard source must contain polygons or multipolygons")
            geometry = ogr.CreateGeometryFromJson(json.dumps(feature["geometry"]))
            if geometry is None or not geometry.IsValid():
                raise ValueError(f"invalid hazard polygon in {key}")
            properties = feature["properties"]
            object_id = properties[entry["object_id_field"]]
            identity = f"{key}:{object_id}"
            if identity in seen:
                raise ValueError(f"duplicate clipped hazard feature: {identity}")
            seen.add(identity)
            record = {"feature_id": identity, "source_object_id": object_id,
                      "label": _label(key, properties), "properties": properties,
                      "geometry_type": feature["geometry"]["type"]}
            catalog[identity] = record
            features.append((record, geometry))
        layer["ogr_features"] = features
        layer_metadata.append({"key": key, "label": entry["label"], "role": LAYER_ROLES[key],
                               "endpoint": entry["endpoint"], "query_status": entry["query_status"],
                               "clip_status": entry["clip_status"], "source_complete": layer["source_complete"],
                               "query_envelope": layer["summary"]["query_envelope"],
                               "available_clip_bounds": layer["available_bounds"], "srid": 26910,
                               "clipped_file": entry["clipped_file"], "clipped_sha256": entry["clipped_sha256"],
                               "feature_ids": [record["feature_id"] for record, _ in features]})
    rows, statuses, labels = [], defaultdict(Counter), defaultdict(Counter)
    comparisons = 0
    for scenario in site["scenarios"]:
        footprints = _absolute_footprints(scenario, origin)
        geometries = [(f, polygon(f["ring"])) for f in footprints]
        results = []
        for layer in layers:
            key = layer["entry"]["key"]
            coverage = {f["id"]: _covered(f["ring"], layer["available_bounds"]) for f in footprints}
            intersections = []
            for record, hazard in layer["ogr_features"]:
                hits = []
                for footprint, geometry in geometries:
                    comparisons += 1
                    if geometry.Intersects(hazard):
                        area = geometry.Intersection(hazard).Area()
                        hits.append({"footprint_id": footprint["id"], "intersection_area_m2": area,
                                     "relationship": "area_overlap" if area > 0 else "boundary_contact"})
                if hits:
                    intersections.append({"feature_id": record["feature_id"],
                                          "footprint_ids": [h["footprint_id"] for h in hits],
                                          "footprint_intersections": hits})
            if not all(coverage.values()):
                status = "unknown_footprint_outside_query_bounds"
            elif not layer["source_complete"]:
                status = "unknown_incomplete_source"
            elif intersections:
                status = "intersects_returned_map_features"
            else:
                status = "no_intersection_with_returned_features"
            statuses[key][status] += 1
            for label in {catalog[r["feature_id"]]["label"] for r in intersections}:
                labels[key][label] += 1
            results.append({"key": key, "status": status,
                            "footprint_query_coverage": coverage, "all_footprints_within_query_bounds": all(coverage.values()),
                            "intersections": intersections})
        rows.append({"scenario_id": scenario["id"], "site_id": scenario["site_id"],
                     "architecture": scenario["architecture"], "footprint_count": len(footprints),
                     "layers": results, "site_approval": "unresolved"})
    site_summaries = []
    for candidate in atlas["candidates"]:
        selected = [row for row in rows if row["site_id"] == candidate["id"]]
        entries = []
        for metadata in layer_metadata:
            results = [layer for row in selected for layer in row["layers"] if layer["key"] == metadata["key"]]
            count = Counter(label for layer in results for label in {catalog[hit["feature_id"]]["label"] for hit in layer["intersections"]})
            entries.append({"key": metadata["key"], "scenario_status_counts": dict(Counter(layer["status"] for layer in results)),
                            "scenario_counts_by_intersecting_label": dict(sorted(count.items()))})
        site_summaries.append({"site_id": candidate["id"], "label": candidate["label"],
                               "scenario_count": len(selected), "layers": entries, "site_approval": "unresolved"})
    if any(path.read_bytes() != blobs[key] for key, path in paths.items()) or _sha(Path(__file__).read_bytes()) != code_sha:
        raise ValueError("hazard inputs or generator changed during screening")
    report = {
        "format": FORMAT, "crs": "EPSG:26910", "origin_added_to_local_footprints": origin,
        "engines": {"gdal": gdal.VersionInfo("RELEASE_NAME"),
                    "geos": f"{ogr.GetGEOSVersionMajor()}.{ogr.GetGEOSVersionMinor()}.{ogr.GetGEOSVersionMicro()}"},
        "input_hashes": hashes, "generator_sha256": code_sha,
        "retrieval_interval_utc": manifest["retrieval_interval_utc"],
        "layers": layer_metadata, "feature_catalog": catalog, "sites": site_summaries, "scenarios": rows,
        "intersection_reference_contract": "Resolve intersections[].feature_id through feature_catalog; complete source properties, label and identity are stored once in that catalogue.",
        "result_interpretation": "Feature labels and source scope govern every result; no hazard absence or site approval inferred.",
        "summary": {"scenario_count": len(rows), "footprints_per_scenario": 4,
                    "layer_count": len(layers), "scenario_layer_results": len(rows) * len(layers),
                    "geometry_pair_tests": comparisons,
                    "all_footprints_within_every_query": all(layer["all_footprints_within_query_bounds"] for row in rows for layer in row["layers"]),
                    "scenario_status_counts_by_layer": {key: dict(value) for key, value in statuses.items()},
                    "scenario_counts_by_intersecting_label": {key: dict(sorted(value.items())) for key, value in labels.items()},
                    "approved_sites": 0},
        "scope": "Intersects each recorded base/access footprint separately, including boundary contact, with Polygon/MultiPolygon geometry and holes preserved. No new placement is proposed.",
        "limitations": ["No returned feature or no intersection is not evidence that a site is hazard-free.",
                        "Outside Hazard Area, Low and reference-map coverage retain their source-specific meanings.",
                        "An incomplete query or any footprint outside available query/clip bounds remains unknown.",
                        "Source FEMA -9999 values remain raw sentinels; no flood-elevation or ground-elevation comparison is performed.",
                        "Source area/length properties retain source units and describe the original features, not these intersections.",
                        "Mapped polygons, dates, model scenarios and applicability need agency review; no foundation, evacuation or construction approval is inferred.",
                        "This is floating-point GIS overlay, not an exact physical hazard boundary certificate."],
    }
    report["content_sha256"] = _digest(report)
    return report


def check_report(report_path=ROOT / "hazard-report.json", atlas_path=ROOT / "data/atlas-subset.json",
                 site_path=ROOT / "site-report.json", manifest_path=ROOT / "data/hazards/manifest.json"):
    """Offline stdlib integrity/coverage replay; does not rerun GEOS predicates."""
    report = json.loads(Path(report_path).read_bytes())
    atlas, site, _, layers, hashes, _, _ = _inputs(atlas_path, site_path, manifest_path)
    if report.get("format") != FORMAT or report.get("content_sha256") != _digest(report):
        raise ValueError("hazard report format or content digest mismatch")
    if report["input_hashes"] != hashes or report["generator_sha256"] != _sha(Path(__file__).read_bytes()):
        raise ValueError("hazard report does not bind current input/generator bytes")
    by_id = {scenario["id"]: scenario for scenario in site["scenarios"]}
    by_layer = {layer["entry"]["key"]: layer for layer in layers}
    expected_catalog = {}
    for layer in layers:
        entry = layer["entry"]
        for feature in layer["features"]:
            properties = feature["properties"]
            object_id = properties[entry["object_id_field"]]
            identity = f'{entry["key"]}:{object_id}'
            expected_catalog[identity] = {"feature_id": identity, "source_object_id": object_id,
                                          "label": _label(entry["key"], properties),
                                          "properties": properties, "geometry_type": feature["geometry"]["type"]}
    if report["feature_catalog"] != expected_catalog:
        raise ValueError("feature catalogue does not preserve the complete source properties and identities")
    if {row["scenario_id"] for row in report["scenarios"]} != set(by_id) or len(report["scenarios"]) != len(by_id):
        raise ValueError("hazard report scenario set differs from study")
    for row in report["scenarios"]:
        footprints = _absolute_footprints(by_id[row["scenario_id"]], atlas["origin"])
        if len(footprints) != 4 or row["footprint_count"] != 4 or row["site_approval"] != "unresolved":
            raise ValueError("unexpected footprint or approval contract")
        if {layer["key"] for layer in row["layers"]} != set(by_layer) or len(row["layers"]) != len(by_layer):
            raise ValueError("incomplete scenario layer set")
        for result in row["layers"]:
            layer = by_layer[result["key"]]
            coverage = {f["id"]: _covered(f["ring"], layer["available_bounds"]) for f in footprints}
            if result["footprint_query_coverage"] != coverage or result["all_footprints_within_query_bounds"] != all(coverage.values()):
                raise ValueError("incorrect footprint query coverage")
            if not all(coverage.values()) and result["status"] != "unknown_footprint_outside_query_bounds":
                raise ValueError("out-of-query footprint must remain unknown")
            if all(coverage.values()) and not layer["source_complete"] and result["status"] != "unknown_incomplete_source":
                raise ValueError("incomplete layer must remain unknown")
            seen = set()
            for hit in result["intersections"]:
                if set(hit) != {"feature_id", "footprint_ids", "footprint_intersections"}:
                    raise ValueError("intersection must use the normalized catalogue-reference schema")
                record = resolve_feature(report, hit)
                identity = record["feature_id"]
                if not identity.startswith(result["key"] + ":") or identity in seen:
                    raise ValueError("wrong-layer or repeated feature reference")
                seen.add(identity)
                footprint_ids = [item["footprint_id"] for item in hit["footprint_intersections"]]
                if (not footprint_ids or hit["footprint_ids"] != footprint_ids
                        or len(set(footprint_ids)) != len(footprint_ids) or not set(footprint_ids) <= set(coverage)):
                    raise ValueError("intersection footprint references are inconsistent")
                for item in hit["footprint_intersections"]:
                    area = item["intersection_area_m2"]
                    if type(area) not in (int, float) or not math.isfinite(area) or area < 0:
                        raise ValueError("intersection area must be finite and nonnegative")
                    if item["relationship"] != ("area_overlap" if area > 0 else "boundary_contact"):
                        raise ValueError("intersection relationship disagrees with recorded area")
            if all(coverage.values()) and layer["source_complete"]:
                expected_status = "intersects_returned_map_features" if result["intersections"] else "no_intersection_with_returned_features"
                if result["status"] != expected_status:
                    raise ValueError("result status disagrees with the recorded intersections")
    return report["summary"]


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=ROOT / "data/atlas-subset.json")
    parser.add_argument("--site-report", type=Path, default=ROOT / "site-report.json")
    parser.add_argument("--manifest", type=Path, default=ROOT / "data/hazards/manifest.json")
    parser.add_argument("--output", type=Path, default=ROOT / "hazard-report.json")
    parser.add_argument("--check", action="store_true")
    args = parser.parse_args()
    if args.check:
        summary = check_report(args.output, args.atlas, args.site_report, args.manifest)
    else:
        report = run_screen(args.atlas, args.site_report, args.manifest)
        args.output.write_bytes(_canonical(report) + b"\n")
        summary = report["summary"]
    print(json.dumps(summary, indent=2, ensure_ascii=False))


if __name__ == "__main__":
    main()
