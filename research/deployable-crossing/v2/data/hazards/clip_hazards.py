#!/usr/bin/env python3
"""Clip audited ArcGIS ESRIJSON responses with GDAL/GEOS; keep raw data external.

No manual ring interpretation, topology repair, generalization or source-attribute
reclassification is performed. All output coordinates use the source query's CRS.
"""

import argparse
from datetime import datetime, timezone
import hashlib
import json
import math
from pathlib import Path
import shutil
import subprocess
import tempfile


def write_json(path, obj):
    path.write_text(json.dumps(obj, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def sha(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def main():
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument("--manifest", type=Path, required=True)
    ap.add_argument("--output-dir", type=Path, required=True)
    args = ap.parse_args()
    source = args.manifest.resolve().parent
    dest = args.output_dir.resolve()
    if source == dest:
        ap.error("Raw source and clipped output directories must differ")
    dest.mkdir(parents=True, exist_ok=True)
    provenance = dest / "query-provenance"
    provenance.mkdir(exist_ok=True)
    raw = json.loads(args.manifest.read_text(encoding="utf-8"))
    bounds = raw["study_envelope"]["bbox"]
    srid = raw["study_envelope"]["srid"]
    gdal_version = subprocess.check_output(["ogr2ogr", "--version"], text=True).strip()
    shutil.copyfile(args.manifest, dest / "raw-manifest.json")
    hash_checks = []
    for request in sorted(source.glob("*.request.json")):
        record = json.loads(request.read_text(encoding="utf-8"))
        shutil.copyfile(request, provenance / request.name)
        if record.get("response_file"):
            body = source / record["response_file"]
            valid = sha(body) == record["sha256"]
            hash_checks.append({"request": "query-provenance/" + request.name, "response_file": record["response_file"], "sha256": record["sha256"], "verified": valid})
            if not valid:
                raise ValueError("Raw response hash mismatch: " + str(body))

    layers = []
    for layer in raw["layers"]:
        key = layer["key"]
        output = {"key": key, "label": layer.get("label"), "endpoint": layer.get("endpoint"),
                  "query_status": layer["status"], "query_feature_count": layer.get("reported_count"),
                  "raw_object_ids": layer.get("object_ids"), "object_id_field": layer.get("object_id_field"),
                  "raw_summary": key + ".summary.json", "clip_status": "unknown", "operations": []}
        if (source / (key + ".summary.json")).exists():
            shutil.copyfile(source / (key + ".summary.json"), dest / (key + ".summary.json"))
        clipped_path = dest / (key + ".clipped.geojson")
        if layer["status"] == "unknown":
            output["reason"] = layer.get("reason", "Source query incomplete or failed")
            if clipped_path.exists():
                clipped_path.unlink()
            layers.append(output)
            continue
        all_features = []
        collection_template = None
        failed = False
        with tempfile.TemporaryDirectory(prefix="hazard-clip-") as temp:
            for i, batch in enumerate(layer.get("features", []), 1):
                input_path = source / batch["file"]
                if sha(input_path) != batch["sha256"]:
                    raise ValueError("Feature response hash mismatch: " + str(input_path))
                interim = Path(temp) / ("batch-%03d.geojson" % i)
                command = ["ogr2ogr", "-if", "ESRIJSON", "-f", "GeoJSON", str(interim), str(input_path),
                           "-clipsrc", *map(str, bounds), "-s_srs", "EPSG:" + str(srid), "-t_srs", "EPSG:" + str(srid),
                           "-nln", key, "-lco", "RFC7946=NO", "-lco", "COORDINATE_PRECISION=15"]
                process = subprocess.run(command, text=True, capture_output=True)
                operation = {"source_file": batch["file"], "source_sha256": batch["sha256"],
                             "command": command, "returncode": process.returncode,
                             "stdout": process.stdout, "stderr": process.stderr,
                             "requested_source_ids": batch["requested_object_ids"]}
                output["operations"].append(operation)
                if process.returncode or not interim.exists():
                    failed = True
                    output["reason"] = "GDAL conversion/clipping failed; retained diagnostic output."
                    break
                dataset = json.loads(interim.read_text(encoding="utf-8"))
                collection_template = {k: v for k, v in dataset.items() if k != "features"}
                features = dataset.get("features", [])
                all_features.extend(features)
                operation["clipped_feature_count"] = len(features)
                operation["clipped_feature_ids"] = [f.get("properties", {}).get(layer["object_id_field"]) for f in features]
        if failed:
            if clipped_path.exists():
                clipped_path.unlink()
            layers.append(output)
            continue
        if collection_template is None:
            collection_template = {"type": "FeatureCollection", "name": key,
                                   "crs": {"type": "name", "properties": {"name": "urn:ogc:def:crs:EPSG::" + str(srid)}}}
        dataset = {**collection_template, "features": all_features}
        write_json(clipped_path, dataset)
        clipped_ids = [f.get("properties", {}).get(layer.get("object_id_field")) for f in all_features]
        if not set(clipped_ids).issubset(set(layer.get("object_ids") or [])):
            raise ValueError("Clipped IDs not in original ID query: " + key)
        output.update({"clip_status": "complete", "clipped_file": clipped_path.name,
                       "clipped_sha256": sha(clipped_path), "clipped_bytes": clipped_path.stat().st_size,
                       "clipped_feature_count": len(all_features), "clipped_object_ids": clipped_ids,
                       "properties": [f.get("properties", {}) for f in all_features]})
        coordinates = []
        def collect_coordinates(value):
            if isinstance(value, list) and len(value) >= 2 and all(isinstance(v, (int, float)) for v in value[:2]):
                coordinates.append(value[:2])
            elif isinstance(value, list):
                for item in value:
                    collect_coordinates(item)
        for feature in all_features:
            collect_coordinates((feature.get("geometry") or {}).get("coordinates", []))
        output["coordinates_within_clip_bounds"] = all(math.isfinite(x) and math.isfinite(y) and bounds[0] - 1e-8 <= x <= bounds[2] + 1e-8 and bounds[1] - 1e-8 <= y <= bounds[3] + 1e-8 for x, y in coordinates)
        if not output["coordinates_within_clip_bounds"]:
            output["clip_status"] = "unknown"
            output["reason"] = "Coordinate-bound validation failed."
        if all_features:
            sql = 'SELECT "' + layer["object_id_field"] + '" AS source_object_id, ST_IsValid(geometry) AS is_valid, ST_Area(geometry) AS area_square_metres FROM "' + key + '"'
            check_command = ["ogr2ogr", "-f", "GeoJSON", "/vsistdout/", str(clipped_path), "-dialect", "SQLite", "-sql", sql]
            check = subprocess.run(check_command, text=True, capture_output=True)
            output["geometry_validation"] = {"command": check_command, "returncode": check.returncode, "stderr": check.stderr}
            try:
                checks = [f["properties"] for f in json.loads(check.stdout)["features"]]
                valid = check.returncode == 0 and len(checks) == len(all_features) and all(v.get("is_valid") == 1 for v in checks)
                output["geometry_validation"].update({"features": checks, "all_geometries_valid": valid})
                if not valid:
                    output["clip_status"] = "unknown"
                    output["reason"] = "GEOS geometry-validity check failed or incomplete."
            except Exception as exc:
                output["clip_status"] = "unknown"
                output["reason"] = "Geometry-validation output could not be interpreted: " + str(exc)
        layers.append(output)
        print(json.dumps({"key": key, "query_count": layer.get("reported_count"), "clipped_count": len(all_features), "bytes": clipped_path.stat().st_size}), flush=True)

    manifest = {"schema_version": 2, "created_at_utc": datetime.now(timezone.utc).isoformat(timespec="seconds"),
                "study_envelope": raw["study_envelope"],
                "retrieval_interval_utc": [raw["started_at_utc"], raw["completed_at_utc"]],
                "raw_storage": {"external_directory": str(source), "included_in_repository": False,
                                "manifest_file": "raw-manifest.json", "request_provenance_directory": "query-provenance"},
                "processing": {"tool": gdal_version, "input_driver": "ESRIJSON", "output_driver": "GeoJSON",
                               "clipping": "GDAL/GEOS rectangle intersection via ogr2ogr -clipsrc; preserves holes and multipart topology through the geometry engine.",
                               "clip_bounds": bounds, "input_srid": srid, "output_srid": srid,
                               "RFC7946": False, "coordinate_precision": 15, "simplification": None,
                               "topology_repair": None, "raw_hash_checks": hash_checks},
                "limitations": raw["limitations"] + ["Clipped data are valid only inside the declared study envelope.",
                              "Returned feature labels must be interpreted; Outside Hazard Area and Low are not hazard-positive classes.",
                              "An empty clipped result is not evidence of complete mapping or a hazard-free site.",
                              "Clipping and CRS conversion use floating-point GIS, not exact geodetic arithmetic.",
                              "Source area/length attributes keep original units and are not recomputed by clipping."],
                "layers": layers}
    write_json(dest / "manifest.json", manifest)
    validation = {"raw_response_hash_count": len(hash_checks), "all_raw_hashes_verified": all(x["verified"] for x in hash_checks),
                  "all_layer_clips_complete": all(x["clip_status"] == "complete" for x in layers),
                  "unknown_layers": [x["key"] for x in layers if x["clip_status"] != "complete"]}
    write_json(dest / "validation.json", validation)


if __name__ == "__main__":
    main()
