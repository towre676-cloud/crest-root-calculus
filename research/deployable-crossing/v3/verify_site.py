#!/usr/bin/python3
"""Independent GDAL/GEOS replay of v3 site and hazard intersection records.

This optional verifier reads the stored polygons and source snapshots. It does
not import site.py, site_geometry.py, the Atlas interpolator or motion solvers.
It checks their recorded GIS conclusions with GEOS. It does not establish
source accuracy, constructability, site permission or physical clearance.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
from pathlib import Path


ROOT = Path(__file__).resolve().parent
V2 = ROOT.parent / "v2"
VEHICLES = {"trunk", "primary", "secondary", "tertiary", "residential", "unclassified", "service", "track"}
FIELDS = ("building_ids", "building_ids_at_3m_sensitivity", "vehicle_centerline_ids", "other_vehicle_centerline_ids", "parcel_ids")


def digest(data):
    return hashlib.sha256(data).hexdigest()


def canonical(value):
    return digest(json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False).encode("utf-8"))


def verify():
    try:
        from osgeo import ogr, gdal
    except ImportError as exc:
        raise SystemExit("Independent GIS verification needs optional GDAL Python bindings; use /usr/bin/python3 verify_site.py") from exc
    ogr.UseExceptions()
    gdal.UseExceptions()
    paths, blobs = {}, {}

    def load(path):
        path = Path(path)
        raw = path.read_bytes()
        paths[str(path)] = path
        blobs[str(path)] = raw
        return json.loads(raw)

    report = load(ROOT / "site-report.json")
    atlas = load(V2 / "data/atlas-subset.json")
    design = load(ROOT / "design.json")
    manifest_path = V2 / "data/hazards/manifest.json"
    manifest = load(manifest_path)
    if report["design_sha256"] != canonical(design) or report["atlas_content_sha256"] != atlas["content_sha256"]:
        raise ValueError("The site report does not identify these design/Atlas inputs")
    if atlas["crs"] != "EPSG:26910" or atlas["origin"] != [430000, 5028500]:
        raise ValueError("Unexpected Atlas reference frame")
    if manifest["schema_version"] != 2 or manifest["study_envelope"]["srid"] != 26910:
        raise ValueError("Unexpected hazard manifest or CRS")

    def polygon(rings):
        geometry = ogr.Geometry(ogr.wkbPolygon)
        for points in rings:
            boundary = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in points:
                boundary.AddPoint_2D(x, y)
            boundary.CloseRings()
            geometry.AddGeometry(boundary)
        return geometry

    def line(points):
        geometry = ogr.Geometry(ogr.wkbLineString)
        for x, y in points:
            geometry.AddPoint_2D(x, y)
        return geometry

    def near(a, b, pad=0):
        x, y = a.GetEnvelope(), b.GetEnvelope()
        return not (x[1] + pad < y[0] or y[1] + pad < x[0] or x[3] + pad < y[2] or y[3] + pad < x[2])

    def intersects(footprints, geometry, distance=0):
        return any(near(foot, geometry, distance) and
                   (foot.Intersects(geometry) if distance == 0 else foot.Distance(geometry) <= distance + 1e-8)
                   for foot in footprints)

    invalid = []

    def checked(identity, geometry):
        if geometry is None or not geometry.IsValid():
            invalid.append(identity)
        return geometry

    buildings = [(row["id"], checked("building:" + row["id"], polygon(row["rings"]))) for row in atlas["buildings"]]
    parcels = [(row["id"], checked("parcel:" + row["id"], polygon(row["rings"]))) for row in atlas["parcels"]]
    roads = [(row, checked("road:" + row["id"], line(row["p"]))) for row in atlas["streets"] if row["kind"] in VEHICLES]
    sites = {row["id"]: row for row in atlas["candidates"]}
    layers = []
    hazard_hashes = {"manifest.json": digest(blobs[str(manifest_path)])}
    for entry in manifest["layers"]:
        source_path = manifest_path.parent / entry["clipped_file"]
        data = load(source_path)
        actual_hash = digest(blobs[str(source_path)])
        if actual_hash != entry["clipped_sha256"] or len(data["features"]) != entry["clipped_feature_count"]:
            raise ValueError("The clipped hazard source differs: " + entry["key"])
        crs = data.get("crs", {}).get("properties", {}).get("name")
        if crs not in ("urn:ogc:def:crs:EPSG::26910", "EPSG:26910"):
            raise ValueError("Unexpected clipped hazard CRS")
        summary_path = manifest_path.parent / entry["raw_summary"]
        summary = load(summary_path)
        if summary["query_srid"] != 26910:
            raise ValueError("Unexpected source query CRS")
        a, b = manifest["processing"]["clip_bounds"], summary["query_envelope"]
        bounds = [max(a[0], b[0]), max(a[1], b[1]), min(a[2], b[2]), min(a[3], b[3])]
        if bounds[0] >= bounds[2] or bounds[1] >= bounds[3]:
            raise ValueError("Query and clip bounds do not overlap")
        coverage = polygon([[[bounds[0], bounds[1]], [bounds[2], bounds[1]], [bounds[2], bounds[3]], [bounds[0], bounds[3]]]])
        features = [(str(row["properties"][entry["object_id_field"]]),
                     checked(entry["key"] + ":" + str(row["properties"][entry["object_id_field"]]),
                             ogr.CreateGeometryFromJson(json.dumps(row["geometry"])))) for row in data["features"]]
        layers.append({"key": entry["key"], "features": features, "coverage": coverage,
                       "complete": entry["clip_status"] == "complete" and entry["query_status"] in
                           ("features_returned", "no_features_returned")})
        hazard_hashes[entry["clipped_file"]] = actual_hash
        hazard_hashes[entry["raw_summary"]] = digest(blobs[str(summary_path)])
    if hazard_hashes != report["hazard_source_sha256"]:
        raise ValueError("The report's hazard identity differs from the independently loaded snapshot")

    mismatches, static_results, sweep_results = [], {}, {}
    totals = Counter()
    records = [(row, "footprints", "direct_hits", "hazards", "static") for row in report["scenarios"]]
    records += [(row, "sweep_footprints", "sweep_possible_hits", "sweep_hazards", "sweep") for row in report["placements"]]
    for row, footprints_key, hit_key, hazard_key, kind in records:
        footprints = [checked(row["id"] + "/" + foot["id"], polygon([foot["ring"]])) for foot in row[footprints_key]]
        expected = {key: [] for key in FIELDS}
        for key, geometry in buildings:
            if intersects(footprints, geometry):
                expected["building_ids"].append(key)
            if intersects(footprints, geometry, 3):
                expected["building_ids_at_3m_sensitivity"].append(key)
        for road, geometry in roads:
            if intersects(footprints, geometry):
                expected["vehicle_centerline_ids"].append(road["id"])
                if road["id"] != sites[row["site_id"]]["road_way_id"]:
                    expected["other_vehicle_centerline_ids"].append(road["id"])
        for key, geometry in parcels:
            if intersects(footprints, geometry):
                expected["parcel_ids"].append(key)
        for key, values in expected.items():
            values = sorted(set(values))
            expected[key] = values
            totals[kind + "_hit_sets"] += 1
            if values != row[hit_key][key]:
                mismatches.append({"record": row["id"], "kind": kind, "field": key,
                                   "reported": row[hit_key][key], "geos": values})
        (static_results if kind == "static" else sweep_results)[row["id"]] = expected
        absolute = [polygon([[[p[0] + atlas["origin"][0], p[1] + atlas["origin"][1]]
                              for p in foot["ring"]]]) for foot in row[footprints_key]]
        reported_layers = {value["layer"]: value for value in row[hazard_key]}
        if set(reported_layers) != {layer["key"] for layer in layers}:
            raise ValueError("Missing or unexpected hazard layer in " + row["id"])
        for layer in layers:
            # GEOS set difference checks coverage, independently of the producer's
            # coordinate-by-coordinate rectangle membership implementation.
            covered = all(geometry.Difference(layer["coverage"]).IsEmpty() for geometry in absolute)
            status = "intersections_evaluated" if covered and layer["complete"] else "unknown_source_coverage"
            ids = sorted(key for key, geometry in layer["features"] if status == "intersections_evaluated" and intersects(absolute, geometry))
            recorded = reported_layers[layer["key"]]
            totals[kind + "_hazard_sets"] += 1
            totals[kind + "_hazard_coverage"] += 1
            totals[kind + "_hazard_status"] += 1
            if covered != recorded["all_footprints_within_query"] or status != recorded["status"] or ids != sorted(value["id"] for value in recorded["features"]):
                mismatches.append({"record": row["id"], "kind": kind, "layer": layer["key"],
                    "geos": {"covered": covered, "status": status, "feature_ids": ids}, "reported": recorded})

    clear = [row for row in report["scenarios"] if not static_results[row["id"]]["building_ids"]
             and not static_results[row["id"]]["vehicle_centerline_ids"]]
    both = [row for row in clear if not sweep_results[row["placement_id"]]["building_ids"]]
    counts = {"placement_count": len(report["placements"]), "scenario_count": len(report["scenarios"]),
              "without_direct_building_or_road_conflict": len(clear),
              "also_without_sweep_bound_building_intersection": len(both), "approved_sites": 0,
              "possible_sweep_building_conflicts": sum(bool(values["building_ids"]) for values in sweep_results.values())}
    if counts != report["summary"]:
        mismatches.append({"field": "summary", "geos": counts, "reported": report["summary"]})
    if len(report["scenarios"]) != 216 or len(report["placements"]) != 18:
        mismatches.append({"field": "expected_scope", "expected": {"scenarios": 216, "placements": 18}})
    if report.get("physical_approval") is not False or any(row["site_approval"] != "unresolved" for row in report["scenarios"]):
        mismatches.append({"field": "physical_approval_semantics"})
    for key, path in paths.items():
        if path.read_bytes() != blobs[key]:
            raise ValueError("Input changed during independent GIS verification: " + key)
    return {
        "format": "crest-root-integrated-site-independent/3", "design_sha256": report["design_sha256"],
        "passed": not mismatches and not invalid,
        "engines": {"gdal": gdal.VersionInfo("RELEASE_NAME"),
                    "geos": f"{ogr.GetGEOSVersionMajor()}.{ogr.GetGEOSVersionMinor()}.{ogr.GetGEOSVersionMicro()}"},
        "comparison": {"static_scenarios": len(report["scenarios"]), "sweep_placements": len(report["placements"]),
                       **dict(totals), "mismatch_count": len(mismatches), "mismatches": mismatches,
                       "invalid_geometries": invalid, "independent_summary": counts},
        "input_sha256": {str(path.relative_to(ROOT.parent)): digest(blobs[key]) for key, path in paths.items()},
        "source_code_sha256": {name: digest((ROOT / name).read_bytes()) for name in ("site.py", "verify_site.py")},
        "method": "GDAL/GEOS intersections with holes and multipart geometry preserved; exact GEOS distance for the 3 m sensitivity; independent set-difference source-coverage checks. Source and GIS calculations use binary64.",
        "scope": "Recorded 2D predicates and coverage only. Does not replay layout construction, terrain extrema, continuous motion or physical clearance. Conservative sweep hits are possible obstructions; no hazard-free or approved classification follows.",
        "physical_approval": False,
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=ROOT / "verification-site.json")
    args = parser.parse_args()
    result = verify()
    args.output.write_text(json.dumps(result, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"passed": result["passed"], "comparison": result["comparison"]}, indent=2))
    if not result["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
