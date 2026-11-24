#!/usr/bin/env python3
"""Extract a pinned, attributed downtown testbed from Tillamook Atlas V5.

Coordinates are local east/north metres in EPSG:26910. Encoding precision
is not a positional accuracy claim. This module does not infer road width,
public access, survey error bounds, ground bearing capacity, or clearance.
All operations use the Python standard library and leave the source intact.
"""
from __future__ import annotations

import argparse
import base64
from collections import Counter
from fractions import Fraction
from functools import lru_cache
import hashlib
import json
import math
from pathlib import Path
import struct
import xml.etree.ElementTree as ET


FORMAT = "crest-root-tillamook-subset/1"
ARCHIVE_SHA256 = "90f06d6f8bd4a67e51350f2ae849a1f16440925493c1ab2b7111241f1e1f79a8"
LEDGER_SHA256 = "91b1a7064e7922f033a2591b1b0c0a591d86b7336a7c7e3244bfbdcb23b0adbc"
DEFAULT_BOUNDS = [3750, 4550, 4350, 5450]
DEFAULT_SUBSET = Path(__file__).parent / "data" / "atlas-subset.json"
OSM_LICENSE = "https://opendatacommons.org/licenses/odbl/1-0/"
OSM_ATTRIBUTION = "© OpenStreetMap contributors"
_TAGS = {"name", "ref", "highway", "lanes", "maxspeed", "oneway", "access",
         "motor_vehicle", "vehicle", "sidewalk", "sidewalk:both", "footway",
         "crossing", "crossing:markings", "width", "surface", "bridge", "layer",
         "kerb", "lit", "wheelchair", "tactile_paving"}
_SITES = (
    ("main-7th", "Main Avenue at 7th Street", "w160662295", ["w1557734242"], "36984060"),
    ("pacific-7th", "Pacific Avenue at 7th Street", "w160662310", ["w1557734242"], "36985866"),
    ("main-9th", "Main Avenue at 9th Street", "w160662295", ["w671168314"], "36959806"),
    ("pacific-9th", "Pacific Avenue at 9th Street", "w160662310", ["w671168314", "w5246810"], "36959805"),
    ("main-11th", "Main Avenue at 11th Street", "w160662295", ["w5250328", "w1557734251"], "36984056"),
    ("pacific-11th", "Pacific Avenue at 11th Street", "w160662310", ["w5250328"], "36985862"),
)


def _sha(data):
    return hashlib.sha256(data).hexdigest()


def _canonical(value):
    return json.dumps(value, sort_keys=True, ensure_ascii=False,
                      separators=(",", ":"), allow_nan=False).encode("utf-8")


def _digest(subset):
    return _sha(_canonical({k: v for k, v in subset.items() if k != "content_sha256"}))


def _json(path):
    return json.loads(path.read_text(encoding="utf-8"))


def verify_source(source_root):
    """Check the pinned V5 manifest and all 707 listed regular files."""
    root = Path(source_root).resolve()
    ledger = (root / "CONTENTS.sha256").read_bytes()
    if _sha(ledger) != LEDGER_SHA256:
        raise ValueError("source is not the pinned Tillamook V5 manifest")
    hashes = {}
    for line in ledger.decode("utf-8").splitlines():
        expected, relative = line.split("  ", 1)
        path = root / relative
        if not path.resolve().is_relative_to(root) or not path.is_file():
            raise ValueError(f"source file missing or outside root: {relative}")
        if _sha(path.read_bytes()) != expected:
            raise ValueError(f"source hash mismatch: {relative}")
        hashes[Path(relative).as_posix()] = expected
    town = _json(root / "Content/data/town.json")
    if town["crs"] != "EPSG:26910" or town["origin"] != [430000, 5028500]:
        raise ValueError("unsupported Atlas CRS or origin")
    survey = town["survey"]
    if (survey["extent"], survey["shape"], survey["scale"]) != ([2500, 4500, 6000, 8000], [1750, 1750], .025):
        raise ValueError("unexpected survey grid or quantization")
    return hashes


def _bbox(points):
    if not points:
        raise ValueError("empty geometry")
    return [min(p[0] for p in points), min(p[1] for p in points),
            max(p[0] for p in points), max(p[1] for p in points)]


def _overlaps(a, b):
    return a[0] <= b[2] and b[0] <= a[2] and a[1] <= b[3] and b[1] <= a[3]


def _inside(p, bounds):
    return bounds[0] <= p[0] <= bounds[2] and bounds[1] <= p[1] <= bounds[3]


def _tags(element):
    return {t.get("k"): t.get("v") for t in element.findall("tag") if t.get("k") in _TAGS}


def _osm(root):
    # This duplicates the original compiler's sorted-file, last-record-wins rule.
    nodes, ways = {}, {}
    for file in sorted((root / "BuildInputs").glob("*.osm")):
        xml = ET.parse(file).getroot()
        nodes.update({n.get("id"): n for n in xml.findall("node")})
        ways.update({w.get("id"): w for w in xml.findall("way")})
    return nodes, ways


def _crop_terrain(root, town, bounds):
    survey = town["survey"]
    west, south, east, north = survey["extent"]
    nx, nz = survey["shape"]
    dx, dz = (east - west) / nx, (north - south) / nz
    if not (west <= bounds[0] < bounds[2] <= east and south <= bounds[1] < bounds[3] <= north):
        raise ValueError("subset must lie inside the unblended survey grid")
    c0, c1 = math.floor((bounds[0] - west) / dx), math.ceil((bounds[2] - west) / dx)
    r0, r1 = math.floor((north - bounds[3]) / dz), math.ceil((north - bounds[1]) / dz)
    data = (root / "Content/data/survey-ground.bin").read_bytes()
    if len(data) != 2 * nx * nz:
        raise ValueError("incorrect source terrain byte count")
    crop = b"".join(data[2 * (r * nx + c0):2 * (r * nx + c1)] for r in range(r0, r1))
    return {"source_file": "Content/data/survey-ground.bin",
            "extent": [west + c0 * dx, north - r1 * dz, west + c1 * dx, north - r0 * dz],
            "shape": [c1 - c0, r1 - r0], "source_window": [c0, r0, c1, r1],
            "cell_size_m": [dx, dz], "scale": survey["scale"],
            "encoding": "base64 of signed little-endian Int16", "row_order": "north-to-south",
            "sample_location": "cell centres", "blended": False,
            "vertical_datum": "unresolved; V2 recipe asserts NAVD88, lineage not independently verified",
            "vertical_accuracy_m": None,
            "quantization_only_half_step_m": .0125,
            "data_sha256": _sha(crop), "data_base64": base64.b64encode(crop).decode("ascii")}


def extract_atlas(source_root, bounds=None):
    """Return a deterministic subset; bounds select complete vector features.

    Vector selection uses bounding-box intersection, so retained geometries
    can extend outside the window. Ring order, holes and all vertices survive.
    No vector clipping, simplification or new accuracy estimate is performed.
    """
    root = Path(source_root).resolve()
    hashes = verify_source(root)
    bounds = list(DEFAULT_BOUNDS if bounds is None else bounds)
    if (len(bounds) != 4 or any(type(x) not in (int, float) or not math.isfinite(x) for x in bounds)
            or not bounds[0] < bounds[2] or not bounds[1] < bounds[3]):
        raise ValueError("bounds must be four finite ordered east/north values")
    consumed = {"Content/refs/Sources.md", "Content/refs/repertoire/compile-streets.py",
                "Content/refs/repertoire/compile-town.py", "Content/refs/repertoire/tillamook.json"}

    def read(relative):
        consumed.add(relative)
        return _json(root / relative)

    town = read("Content/data/town.json")
    features = read("Content/data/streets/source-features.json")
    metadata = read("Content/data/streets.json")
    for name, expected in metadata["inputSha256"].items():
        relative = f"BuildInputs/{name}"
        consumed.add(relative)
        if hashes.get(relative) != expected:
            raise ValueError(f"street input provenance disagrees: {name}")
    nodes, ways = _osm(root)
    paths_by_id = {f["id"]: f for f in features["paths"]}
    streets = [dict(f, osm_tags=_tags(ways[f["id"][1:]])) for f in features["paths"]
               if _overlaps(_bbox(f["p"]), bounds)]
    marks = [dict(f, osm_tags=_tags(nodes[f["id"][1:]])) for f in features["marks"] if _inside(f["p"], bounds)]
    areas = [f for f in features["areas"] if _overlaps(_bbox(f["p"]), bounds)]
    candidates = []
    for key, label, spine_id, cross_ids, node_id in _SITES:
        spine = paths_by_id[spine_id]
        refs = [nd.get("ref") for nd in ways[spine_id[1:]].findall("nd")]
        if len(refs) != len(spine["p"]):
            raise ValueError("source road nodes and compiled vertices differ")
        index = refs.index(node_id)
        point = spine["p"][index]
        if not _inside(point, bounds):
            raise ValueError(f"bounds exclude required candidate {key}")
        for cross_id in cross_ids:
            cross_refs = [nd.get("ref") for nd in ways[cross_id[1:]].findall("nd")]
            if node_id not in cross_refs or point not in paths_by_id[cross_id]["p"]:
                raise ValueError(f"candidate is not a shared source node: {key}")
        a, b = spine["p"][max(0, index - 1)], spine["p"][min(len(refs) - 1, index + 1)]
        delta = [b[i] - a[i] for i in (0, 1)]
        length = math.hypot(*delta)
        if not length or delta[1] == 0:
            raise ValueError("candidate tangent cannot be oriented north")
        sign = 1 if delta[1] > 0 else -1
        node = nodes[node_id]
        candidates.append({"id": key, "label": label, "point": point,
                           "node_id": f"n{node_id}", "road_way_id": spine_id,
                           "cross_street_way_ids": cross_ids,
                           "latlon": [float(node.get("lat")), float(node.get("lon"))],
                           "osm_node_version": int(node.get("version")),
                           "osm_node_timestamp": node.get("timestamp"),
                           "road_tangent": [round(sign * x / length, 15) for x in delta],
                           "road_tangent_points": [a, b],
                           "road_tangent_method": "unit secant through immediate road neighbours, oriented north",
                           "status": "mapped junction for research screening; no crossing design selected"})
    buildings, canopy = [], []
    for file in sorted((root / "Content/data").glob("buildings-*.json")):
        relative = file.relative_to(root).as_posix()
        for feature in _json(file):
            if _overlaps(_bbox([p for ring in feature["rings"] for p in ring]), bounds):
                consumed.add(relative)
                buildings.append({"id": feature["id"], "rings": feature["rings"], "source_file": relative})
    for file in sorted((root / "Content/data").glob("canopy-*.json")):
        relative = file.relative_to(root).as_posix()
        for index, (x, n, h, radius) in enumerate(_json(file)):
            if _overlaps([x - radius, n - radius, x + radius, n + radius], bounds):
                consumed.add(relative)
                canopy.append({"id": f"{file.stem}:{index}", "p": [x, n],
                               "height_estimate_m": h, "radius_estimate_m": radius,
                               "source_file": relative, "source_index": index,
                               "evidence": "lidar/imagery-derived display proxy; species and precise envelope unverified"})
    parcels = [f for f in read("Content/data/streets/parcels.json")
               if _overlaps(_bbox([p for ring in f["rings"] for p in ring]), bounds)]
    environment = []
    for f in read("Content/data/environment/catalogue.json"):
        points = [f["p"]] if f["geometry"] == "point" else f["p"]
        if _overlaps(_bbox(points), bounds):
            environment.append(f)
    consumed.add("Content/data/survey-ground.bin")
    terrain = _crop_terrain(root, town, bounds)
    payload = {
        "format": FORMAT, "crs": "EPSG:26910", "origin": [430000, 5028500],
        "axis_order": ["east", "north"], "units": "metres", "bounds": bounds,
        "vector_selection": "complete feature whose bounding box intersects bounds; no clipping",
        "provenance": {
            "source_archive": "Tillamook-Atlas-V5-Staging.zip", "archive_sha256": ARCHIVE_SHA256,
            "manifest_sha256": LEDGER_SHA256, "manifest_files_verified": len(hashes),
            "source_file_sha256": {name: hashes[name] for name in sorted(consumed)},
            "osm_snapshot_date_as_documented": "2026-11-20",
            "geometry_encoding_step_m": .01, "osm_latlon_encoding_degrees": .0000001,
            "horizontal_accuracy_m": None,
            "osm_attribution": OSM_ATTRIBUTION, "osm_database_license": OSM_LICENSE,
            "osm_license_scope": "streets, marks, areas, environment and candidate OSM source records",
            "other_sources": {
                "buildings": "Oregon DOGAMI O-23-01 footprint IDs and rings; no authored facades or roof meshes",
                "parcels": "ODOT Tillamook County Taxlots snapshot; contours are not ownership/access authorization",
                "terrain_canopy": "paired DOGAMI DEM/DSM derivatives; heights are not independent survey controls"},
            "other_source_license_status": "retain source attribution; this subset grants no new rights to third-party inputs",
            "restrictions": ["No encoding step is a survey accuracy bound.",
                             "Authored widths remain authored and must not become engineering span evidence.",
                             "Missing mapped features do not prove clear ground, public access, or absent utilities.",
                             "Vertical datum and source epochs need independent engineering confirmation."]},
        "candidates": candidates, "streets": streets, "marks": marks, "areas": areas,
        "buildings": buildings, "parcels": parcels, "environment": environment,
        "canopy": canopy, "terrain": terrain,
    }
    payload["counts"] = {name: len(payload[name]) for name in
                         ("candidates", "streets", "marks", "areas", "buildings", "parcels", "environment", "canopy")}
    payload["width_evidence_counts"] = dict(sorted(Counter(x["widthEvidence"] for x in streets).items()))
    payload["content_sha256"] = _digest(payload)
    validate_subset(payload)
    return payload


@lru_cache(maxsize=4)
def _terrain_bytes(encoded):
    return base64.b64decode(encoded, validate=True)


def validate_subset(subset):
    """Reject a wrong frame or damaged payload; a digest is not a survey certificate."""
    if subset.get("format") != FORMAT or subset.get("crs") != "EPSG:26910" or subset.get("origin") != [430000, 5028500]:
        raise ValueError("wrong subset format, CRS or origin")
    provenance = subset.get("provenance", {})
    if provenance.get("archive_sha256") != ARCHIVE_SHA256 or provenance.get("manifest_sha256") != LEDGER_SHA256:
        raise ValueError("wrong source archive or manifest identity")
    if subset.get("content_sha256") != _digest(subset):
        raise ValueError("subset content hash mismatch")
    terrain = subset["terrain"]
    raw = _terrain_bytes(terrain["data_base64"])
    nx, nz = terrain["shape"]
    if (type(nx) is not int or type(nz) is not int or min(nx, nz) <= 0 or len(raw) != 2 * nx * nz
            or _sha(raw) != terrain["data_sha256"]):
        raise ValueError("invalid terrain shape or byte hash")
    if (terrain["scale"] != .025 or terrain["row_order"] != "north-to-south"
            or terrain["sample_location"] != "cell centres" or terrain["blended"] is not False
            or terrain["cell_size_m"] != [2.0, 2.0]):
        raise ValueError("unexpected terrain contract")
    for name, count in subset["counts"].items():
        if len(subset[name]) != count:
            raise ValueError(f"subset count mismatch: {name}")
    if len(subset["candidates"]) != 6:
        raise ValueError("six mapped junctions are required")
    return True


def load_subset(path=DEFAULT_SUBSET):
    subset = _json(Path(path))
    validate_subset(subset)
    return subset


def sample_terrain_exact(subset, x, northing):
    """Bilinear interpolation of encoded heights, returned as an exact Fraction.

    This is exact relative to stored samples, not to physical ground. Values
    within a half cell of the crop edge use its nearest sample; outside the
    crop is rejected. The function never samples the game's blended terrain.
    """
    t = subset["terrain"]
    west, south, east, north = map(lambda z: Fraction(str(z)), t["extent"])
    if isinstance(x, bool) or isinstance(northing, bool):
        raise ValueError("terrain coordinates must be finite numbers")
    x, n = Fraction(str(x)), Fraction(str(northing))
    if not west <= x <= east or not south <= n <= north:
        raise ValueError("terrain coordinate lies outside the crop")
    nx, nz = t["shape"]
    u = min(Fraction(nx - 1), max(Fraction(0), (x - west) * nx / (east - west) - Fraction(1, 2)))
    v = min(Fraction(nz - 1), max(Fraction(0), (north - n) * nz / (north - south) - Fraction(1, 2)))
    i, j = math.floor(u), math.floor(v)
    k, l = min(i + 1, nx - 1), min(j + 1, nz - 1)
    a, b = u - i, v - j
    raw = _terrain_bytes(t["data_base64"])
    value = lambda col, row: struct.unpack_from("<h", raw, 2 * (row * nx + col))[0]
    result = (1 - b) * ((1 - a) * value(i, j) + a * value(k, j)) + b * ((1 - a) * value(i, l) + a * value(k, l))
    return result * Fraction(str(t["scale"]))


def sample_terrain(subset, x, northing):
    """Return the unblended encoded-grid interpolation in metres as a float."""
    return float(sample_terrain_exact(subset, x, northing))


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("source_root", nargs="?", type=Path)
    parser.add_argument("--output", type=Path, default=DEFAULT_SUBSET)
    parser.add_argument("--check", action="store_true", help="validate the existing subset without the source archive")
    args = parser.parse_args()
    if args.check:
        subset = load_subset(args.output)
    else:
        if args.source_root is None:
            parser.error("source_root is required unless --check is used")
        subset = extract_atlas(args.source_root)
        args.output.parent.mkdir(parents=True, exist_ok=True)
        args.output.write_bytes(_canonical(subset) + b"\n")
    print(json.dumps({"counts": subset["counts"], "content_sha256": subset["content_sha256"],
                      "terrain_shape": subset["terrain"]["shape"]}, sort_keys=True))


if __name__ == "__main__":
    main()
