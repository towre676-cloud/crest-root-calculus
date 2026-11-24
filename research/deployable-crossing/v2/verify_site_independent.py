#!/usr/bin/python3
"""Optional independent GDAL/GEOS replay of the recorded site overlap results.

Run with /usr/bin/python3 on a machine with osgeo installed. This dependency
is deliberately outside the standard-library verification path. The verifier
does not import site_geometry, site_study, or the Atlas interpolation code.
It reads the existing report, independently reconstructs its polygon/line/
point predicates, and binds the result to the exact input file bytes.
"""
from __future__ import annotations

import argparse
from collections import Counter
import hashlib
import json
import math
from pathlib import Path


ROOT = Path(__file__).parent
VEHICLES = {"trunk", "primary", "secondary", "tertiary", "residential", "unclassified", "service", "track"}
FOOTWAYS = {"footway", "path", "steps", "pedestrian", "cycleway"}
FIELDS = ("building_ids", "building_ids_at_3m_sensitivity", "vehicle_centerline_ids",
          "vehicle_width_proxy_ids", "mapped_footway_ids", "parcel_ids", "environment_ids", "canopy_proxy_ids")


def sha(data):
    return hashlib.sha256(data).hexdigest()


def verify(atlas_path, report_path):
    try:
        from osgeo import gdal, ogr, osr
    except ImportError as exc:
        raise SystemExit("Optional verifier requires GDAL Python bindings; run /usr/bin/python3 verify_site_independent.py") from exc
    gdal.UseExceptions()
    ogr.UseExceptions()
    osr.UseExceptions()
    atlas_path, report_path = Path(atlas_path), Path(report_path)
    atlas_bytes, report_bytes = atlas_path.read_bytes(), report_path.read_bytes()
    data, report = json.loads(atlas_bytes), json.loads(report_bytes)
    source_paths = [ROOT / name for name in ("site_geometry.py", "site_study.py", "atlas_adapter.py", "verify_site_independent.py")]
    code_hashes = {p.name: sha(p.read_bytes()) for p in source_paths}
    if data["crs"] != "EPSG:26910" or data["origin"] != [430000, 5028500]:
        raise ValueError("independent check is scoped to the declared Tillamook frame")
    if report["atlas_content_sha256"] != data["content_sha256"]:
        raise ValueError("recorded study identifies a different Atlas payload")

    def polygon(rings):
        out = ogr.Geometry(ogr.wkbPolygon)
        for points in rings:
            ring = ogr.Geometry(ogr.wkbLinearRing)
            for x, y in points:
                ring.AddPoint_2D(x, y)
            ring.CloseRings()
            out.AddGeometry(ring)
        return out

    def line(points):
        out = ogr.Geometry(ogr.wkbLineString if len(points) > 1 else ogr.wkbPoint)
        for x, y in points:
            out.AddPoint_2D(x, y)
        return out

    def nearby(g, h, pad=0):
        a, b = g.GetEnvelope(), h.GetEnvelope()  # min-x, max-x, min-y, max-y
        return not (a[1] + pad < b[0] or b[1] + pad < a[0] or a[3] + pad < b[2] or b[3] + pad < a[2])

    def touches(polygons, geometry, pad=0):
        # GEOS Distance implements the circular distance sensitivity directly;
        # no polygon approximation to a round buffer is introduced.
        return any(nearby(p, geometry, pad) and
                   (p.Intersects(geometry) if pad == 0 else p.Distance(geometry) <= pad + 1e-8)
                   for p in polygons)

    buildings = [(r["id"], polygon(r["rings"])) for r in data["buildings"]]
    roads = [(r, line(r["p"])) for r in data["streets"]]
    parcels = [(r["id"], polygon(r["rings"])) for r in data["parcels"]]
    environment = [(r["id"], polygon([r["p"]]) if r["geometry"] == "area" else
                    line([r["p"]] if r["geometry"] == "point" else r["p"])) for r in data["environment"]]
    canopy = [(r["id"], line([r["p"]]), r["radius_estimate_m"]) for r in data["canopy"]]
    sites = {r["id"]: r for r in data["candidates"]}
    invalid = {
        "buildings": [key for key, geometry in buildings if not geometry.IsValid()],
        "parcels": [key for key, geometry in parcels if not geometry.IsValid()],
        "environment": [key for key, geometry in environment if not geometry.IsValid()],
        "scenario_footprints": [],
    }
    mismatches, unconflicted = [], []
    independent_screens = Counter()
    for row in report["scenarios"]:
        polygons = [polygon([f["ring"]]) for f in row["footprints"]]
        invalid["scenario_footprints"].extend(f'{row["id"]}/{f["id"]}' for f, p in zip(row["footprints"], polygons) if not p.IsValid())
        expected = {key: [] for key in FIELDS}
        for key, geometry in buildings:
            if touches(polygons, geometry):
                expected["building_ids"].append(key)
            if touches(polygons, geometry, 3):
                expected["building_ids_at_3m_sensitivity"].append(key)
        for road, geometry in roads:
            if road["kind"] in VEHICLES:
                if touches(polygons, geometry):
                    expected["vehicle_centerline_ids"].append(road["id"])
                width = row["width_hypothesis_m"] if road["id"] == sites[row["site_id"]]["road_way_id"] else road["width"]
                if touches(polygons, geometry, width / 2):
                    expected["vehicle_width_proxy_ids"].append(road["id"])
            elif road["kind"] in FOOTWAYS and touches(polygons, geometry):
                expected["mapped_footway_ids"].append(road["id"])
        for key, geometry in parcels:
            if touches(polygons, geometry):
                expected["parcel_ids"].append(key)
        for key, geometry in environment:
            if touches(polygons, geometry):
                expected["environment_ids"].append(key)
        for key, geometry, radius in canopy:
            if touches(polygons, geometry, radius):
                expected["canopy_proxy_ids"].append(key)
        for key, values in expected.items():
            expected[key] = sorted(set(values))
            if expected[key] != row["hits"][key]:
                mismatches.append({"scenario": row["id"], "field": key,
                                   "reported": row["hits"][key], "geos": expected[key]})
        clear = not expected["building_ids"] and not expected["vehicle_centerline_ids"]
        screen = "no_tested_building_or_road_centerline_conflict" if clear else "conflict_in_declared_layout"
        independent_screens[screen] += 1
        if screen != row["mapped_obstacle_screen"] or row.get("site_approval") != "unresolved":
            mismatches.append({"scenario": row["id"], "field": "screen_or_approval_semantics"})
        if clear:
            distance, building, footprint_index = min(
                (p.Distance(geometry), key, index)
                for index, p in enumerate(polygons) for key, geometry in buildings)
            unconflicted.append({"scenario": row["id"], "architecture": row["architecture"],
                                 "nearest_building_id": building,
                                 "nearest_footprint_id": row["footprints"][footprint_index]["id"],
                                 "nearest_building_distance_m": distance,
                                 "building_within_3m_sensitivity": bool(expected["building_ids_at_3m_sensitivity"]),
                                 "site_approval": "unresolved"})

    geographic, projected = osr.SpatialReference(), osr.SpatialReference()
    geographic.ImportFromEPSG(4326)
    projected.ImportFromEPSG(26910)
    geographic.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    projected.SetAxisMappingStrategy(osr.OAMS_TRADITIONAL_GIS_ORDER)
    transform = osr.CoordinateTransformation(geographic, projected)
    controls = []
    for candidate in data["candidates"]:
        lat, lon = candidate["latlon"]
        east, north, _ = transform.TransformPoint(lon, lat)
        local = [east - 430000, north - 5028500]
        rounded = [round(v, 2) for v in local]
        tangent = candidate["road_tangent"]
        control = {"candidate": candidate["id"], "node_id": candidate["node_id"],
                   "source_latitude_longitude": [lat, lon], "gdal_local_metres": local,
                   "rounded_local_metres": rounded, "stored_local_metres": candidate["point"],
                   "rounding_residual_m": math.dist(local, candidate["point"]),
                   "rounded_coordinates_match": rounded == candidate["point"],
                   "tangent_norm": math.hypot(*tangent), "tangent_points_north": tangent[1] > 0,
                   "tangent_azimuth_clockwise_from_north_degrees": math.degrees(math.atan2(*tangent))}
        controls.append(control)
        if rounded != candidate["point"] or abs(math.hypot(*tangent) - 1) > 1e-12 or tangent[1] <= 0:
            mismatches.append({"candidate": candidate["id"], "field": "coordinate_transform_or_tangent"})
    if report["summary"]["scenario_count"] != len(report["scenarios"]) or report["summary"]["counts_by_screen"] != dict(independent_screens):
        mismatches.append({"field": "summary_counts"})
    if atlas_path.read_bytes() != atlas_bytes or report_path.read_bytes() != report_bytes:
        raise ValueError("checked input changed during verification")
    if any(sha(p.read_bytes()) != code_hashes[p.name] for p in source_paths):
        raise ValueError("associated source code changed during verification")
    distances = [r["nearest_building_distance_m"] for r in unconflicted]
    return {
        "format": "crest-root-site-independent-verification/1",
        "passed": not mismatches and not any(invalid.values()),
        "engines": {"gdal": gdal.VersionInfo("RELEASE_NAME"),
                    "geos": f"{ogr.GetGEOSVersionMajor()}.{ogr.GetGEOSVersionMinor()}.{ogr.GetGEOSVersionMicro()}"},
        "input_hashes": {"atlas_subset_sha256": sha(atlas_bytes), "site_report_sha256": sha(report_bytes),
                         "atlas_content_sha256": data["content_sha256"]},
        "associated_source_code_sha256": code_hashes,
        "comparison": {"scenarios": len(report["scenarios"]), "hit_fields_per_scenario": len(FIELDS),
                       "hit_sets_compared": len(report["scenarios"]) * len(FIELDS),
                       "mismatch_count": len(mismatches), "mismatches": mismatches,
                       "invalid_geometries": invalid, "independent_screen_counts": dict(independent_screens),
                       "unconflicted_count": len(unconflicted),
                       "unconflicted_nearest_building_distance_range_m": [min(distances), max(distances)] if distances else None,
                       "all_unconflicted_cases_have_building_within_3m": bool(unconflicted) and all(r["building_within_3m_sensitivity"] for r in unconflicted)},
        "cases_without_tested_building_or_road_centerline_conflict": unconflicted,
        "candidate_crs_controls": controls,
        "coordinate_control_scope": "Internal OSM-to-Atlas transform checks, not independent withheld survey controls or an accuracy estimate.",
        "method": "GDAL/GEOS polygon intersection and distance for each separately filled footprint; all rings retained; separate line/point tests.",
        "scope": "Independently checks stored overlap results and coordinate transforms. Does not approve a site, certify source accuracy, rerun deployment mechanics, or prove terrain extrema.",
        "replay": "/usr/bin/python3 verify_site_independent.py",
    }


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--atlas", type=Path, default=ROOT / "data/atlas-subset.json")
    parser.add_argument("--report", type=Path, default=ROOT / "site-report.json")
    parser.add_argument("--output", type=Path, default=ROOT / "verification-site-independent.json")
    args = parser.parse_args()
    receipt = verify(args.atlas, args.report)
    args.output.write_text(json.dumps(receipt, indent=2, ensure_ascii=False, allow_nan=False) + "\n", encoding="utf-8")
    print(json.dumps({"passed": receipt["passed"], "engines": receipt["engines"],
                      "comparison": receipt["comparison"], "input_hashes": receipt["input_hashes"]}, indent=2))
    if not receipt["passed"]:
        raise SystemExit(1)


if __name__ == "__main__":
    main()
