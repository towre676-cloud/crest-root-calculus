#!/usr/bin/env python3
"""Retrieve official hazard subsets with byte-level provenance; no site certification.

Run from any working directory with Python 3. The original HTTP response bytes,
request parameters, timestamps and SHA-256 hashes are retained beside this file.
The envelope is a study-area filter, not a surveyed bridge/support footprint.
"""

import argparse
import concurrent.futures
import hashlib
import json
from pathlib import Path
import sys
import subprocess
import urllib.error
import urllib.parse
import urllib.request
from datetime import datetime, timezone

ROOT = Path(__file__).resolve().parent
BBOX = [433750, 5033050, 434350, 5034000]
WKID = 26910
LAYERS = [
    ("fema_nfhl_28", "FEMA flood hazard zones", "https://hazards.fema.gov/arcgis/rest/services/public/NFHL/MapServer/28"),
    ("dogami_tsunami_12", "DOGAMI tsunami evacuation zones", "https://gis.dogami.oregon.gov/arcgis/rest/services/Public/OTGC1/MapServer/12"),
    ("dogami_tsunami_78", "DOGAMI XXL1 Cascadia local-source inundation", "https://gis.dogami.oregon.gov/arcgis/rest/services/Public/OTGC1/MapServer/78"),
    ("dogami_landslide_3", "DOGAMI mapped landslide deposits", "https://gis.dogami.oregon.gov/arcgis/rest/services/Public/SLIDO42/MapServer/3"),
    ("dogami_ls_coverage_0", "DOGAMI detailed susceptibility reference-map coverage", "https://gis.dogami.oregon.gov/arcgis/rest/services/Public/SLIDO42_DetailedLSSus/MapServer/0"),
    ("dogami_ls_deep_1", "DOGAMI detailed deep-landslide susceptibility", "https://gis.dogami.oregon.gov/arcgis/rest/services/Public/SLIDO42_DetailedLSSus/MapServer/1"),
]


def utc_now():
    return datetime.now(timezone.utc).isoformat(timespec="seconds")


def save_json(path, value):
    path.write_text(json.dumps(value, indent=2, ensure_ascii=False) + "\n", encoding="utf-8")


def request_json(prefix, endpoint, params):
    query = urllib.parse.urlencode(params)
    url = endpoint + "?" + query
    record = {"endpoint": endpoint, "method": "GET", "parameters": params,
              "request_url": url, "started_at_utc": utc_now()}
    body = None
    payload = None
    try:
        req = urllib.request.Request(url, headers={"User-Agent": "Crest-Root-research/2.0 (public-data-audit)", "Accept": "application/json"})
        with urllib.request.urlopen(req, timeout=45) as response:
            body = response.read()
            record.update({"http_status": response.status, "response_url": response.url,
                           "response_headers": dict(response.headers.items())})
    except urllib.error.HTTPError as exc:
        body = exc.read()
        record.update({"http_status": exc.code, "error": str(exc),
                       "response_headers": dict(exc.headers.items()) if exc.headers else {}})
    except Exception as exc:
        record.update({"http_status": None, "error": type(exc).__name__ + ": " + str(exc)})
    record["retrieved_at_utc"] = utc_now()
    if body is not None:
        try:
            payload = json.loads(body)
            suffix = ".raw.json"
        except Exception as exc:
            suffix = ".raw.txt"
            record["parse_error"] = type(exc).__name__ + ": " + str(exc)
        path = ROOT / (prefix + suffix)
        path.write_bytes(body)
        record.update({"response_file": path.name, "response_bytes": len(body),
                       "sha256": hashlib.sha256(body).hexdigest()})
        if isinstance(payload, dict) and "error" in payload:
            record["service_error"] = payload["error"]
    record["success"] = record.get("http_status") == 200 and isinstance(payload, dict) and "error" not in payload
    save_json(ROOT / (prefix + ".request.json"), record)
    return payload if record["success"] else None, record


def retrieve_layer(spec):
    key, label, endpoint = spec
    result = {"key": key, "label": label, "endpoint": endpoint,
              "query_envelope": BBOX, "query_srid": WKID,
              "status": "unknown", "requests": [], "features": [],
              "interpretation": "A spatial query is screening evidence, not a hazard-free or buildable-site finding."}

    def get(tag, suffix, params):
        payload, record = request_json(key + "." + tag, endpoint + suffix, params)
        result["requests"].append(key + "." + tag + ".request.json")
        if payload is None and record.get("http_status") in (None, 502, 503, 504):
            retry_prefix = key + "." + tag + ".retry1"
            payload, record = request_json(retry_prefix, endpoint + suffix, params)
            result["requests"].append(retry_prefix + ".request.json")
        return payload, record

    metadata, _ = get("layer", "", {"f": "json"})
    if metadata is None:
        result["reason"] = "Layer metadata could not be retrieved. Exposure and coverage unknown."
        return result
    result["layer_metadata"] = {k: metadata.get(k) for k in ["name", "type", "geometryType", "objectIdField", "maxRecordCount", "extent", "editingInfo"]}
    base = {"f": "json", "where": "1=1", "geometryType": "esriGeometryEnvelope",
            "geometry": json.dumps({"xmin": BBOX[0], "ymin": BBOX[1], "xmax": BBOX[2], "ymax": BBOX[3], "spatialReference": {"wkid": WKID}}, separators=(",", ":")),
            "inSR": str(WKID), "spatialRel": "esriSpatialRelIntersects"}
    ids_response, _ = get("ids", "/query", {**base, "returnIdsOnly": "true"})
    count_response, _ = get("count", "/query", {**base, "returnCountOnly": "true"})
    if ids_response is None or count_response is None:
        result["reason"] = "Count or ID query failed. No empty-result inference is allowed."
        return result
    ids = ids_response.get("objectIds")
    count = count_response.get("count")
    result["reported_count"] = count
    result["object_id_field"] = ids_response.get("objectIdFieldName") or metadata.get("objectIdField")
    if ids is None and "objectIds" in ids_response and count == 0:
        ids = []
        result["empty_id_normalization"] = "The successful ID query returned objectIds:null and the independent successful count query returned count:0. Both original responses are retained."
    result["object_ids"] = ids
    if not isinstance(ids, list) or not isinstance(count, int) or len(ids) != count or len(set(ids)) != count or ids_response.get("exceededTransferLimit") or count_response.get("exceededTransferLimit"):
        result["reason"] = "Count/ID control is incomplete or inconsistent; possible transfer limit or changing source."
        return result
    result["count_id_control_passed"] = True
    if not ids:
        result["status"] = "no_features_returned"
        result["reason"] = "Successful envelope query returned zero features in this layer; this does not establish hazard absence, mapping completeness or suitability."
        return result
    # Small batches also limit the cost of retrying a transient feature-service failure.
    batch_size = min(int(metadata.get("maxRecordCount") or 4), 4)
    fetched = []
    oid_field = result["object_id_field"]
    for batch_number, start in enumerate(range(0, len(ids), batch_size), 1):
        selected = sorted(ids)[start:start + batch_size]
        params = {"f": "json", "objectIds": ",".join(map(str, selected)), "outFields": "*",
                  "returnGeometry": "true", "outSR": str(WKID), "returnZ": "false", "returnM": "false"}
        payload, record = get("features_%03d" % batch_number, "/query", params)
        if payload is None:
            result["reason"] = "Full-geometry batch failed. Any retained partial features do not complete the query."
            return result
        features = payload.get("features", [])
        batch_ids = [feature.get("attributes", {}).get(oid_field) for feature in features]
        batch_valid = not payload.get("exceededTransferLimit") and set(batch_ids) == set(selected) and len(batch_ids) == len(selected)
        result["features"].append({"file": record["response_file"], "sha256": record["sha256"],
                                   "count": len(features), "spatialReference": payload.get("spatialReference"),
                                   "requested_object_ids": selected, "returned_object_ids": batch_ids,
                                   "exceededTransferLimit": payload.get("exceededTransferLimit", False),
                                   "batch_id_control_passed": batch_valid})
        if not batch_valid:
            result["reason"] = "Feature-batch IDs or transfer-limit check failed; retained response is incomplete."
            return result
        fetched.extend(batch_ids)
    result["returned_feature_count"] = len(fetched)
    if set(fetched) != set(ids) or len(fetched) != len(ids):
        result["reason"] = "Returned feature IDs do not match the envelope ID query."
        return result
    result["status"] = "features_returned"
    result["reason"] = "All features selected by this envelope query were retrieved with full XY geometry. Assess their attributes and actual proposed footprints; the study envelope is not a bridge site."
    return result


def main():
    global ROOT, BBOX
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output-dir", type=Path, required=True, help="External directory for original HTTP response bodies and request audit metadata.")
    parser.add_argument("--bbox", type=float, nargs=4, default=BBOX, metavar=("XMIN", "YMIN", "XMAX", "YMAX"))
    parser.add_argument("--clip-output-dir", type=Path, help="Optional repo directory for clipped GeoJSON and copied raw/query provenance, without raw statewide geometries.")
    args = parser.parse_args()
    ROOT = args.output_dir.resolve()
    BBOX = args.bbox
    if BBOX[0] >= BBOX[2] or BBOX[1] >= BBOX[3]:
        parser.error("bbox must have positive width and height")
    ROOT.mkdir(parents=True, exist_ok=True)
    started = utc_now()
    results = []
    with concurrent.futures.ThreadPoolExecutor(max_workers=3) as executor:
        future_map = {executor.submit(retrieve_layer, spec): spec[0] for spec in LAYERS}
        for future in concurrent.futures.as_completed(future_map):
            key = future_map[future]
            try:
                result = future.result()
            except Exception as exc:
                result = {"key": key, "status": "unknown", "reason": type(exc).__name__ + ": " + str(exc)}
            results.append(result)
            save_json(ROOT / (key + ".summary.json"), result)
            print(json.dumps({"key": key, "status": result["status"], "count": result.get("reported_count"), "reason": result.get("reason")}), flush=True)
    ordered = {r["key"]: r for r in results}
    manifest = {"schema_version": 1, "started_at_utc": started, "completed_at_utc": utc_now(),
                "study_envelope": {"bbox": BBOX, "srid": WKID, "units": "metres", "source": "Parent task's downtown candidate corridor, not a surveyed site boundary"},
                "query_method": "Envelope-intersection count and complete-ID controls; full XY geometry by bounded object-ID batches; no generalization or coordinate rounding requested.",
                "limitations": ["Successful zero-feature queries are not findings of hazard absence.",
                                "Service failures and incomplete batches remain unknown.",
                                "Full feature geometries may extend beyond the study envelope.",
                                "Count/ID queries are sequential observations, not an atomic source snapshot.",
                                "Coverage, source dates, vertical datum and site applicability require review.",
                                "No ground capacity, seismic class, road clearance or construction approval is inferred."],
                "layers": [ordered[s[0]] for s in LAYERS]}
    save_json(ROOT / "manifest.json", manifest)
    if args.clip_output_dir is not None:
        clip_script = Path(__file__).resolve().with_name("clip_hazards.py")
        subprocess.run([sys.executable, str(clip_script), "--manifest", str(ROOT / "manifest.json"), "--output-dir", str(args.clip_output_dir)], check=True)
    return 0


if __name__ == "__main__":
    sys.exit(main())
