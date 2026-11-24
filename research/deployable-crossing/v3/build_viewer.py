#!/usr/bin/env python3
"""Build one offline review surface from the integrated, hash-bound v3 reports.

Only Python's standard library is required. No page is produced before all five
reports and the reference geometry are available. The viewer does not solve or
invent a second design: it draws the samples exported by the motion report.
"""
from __future__ import annotations

import argparse
import hashlib
import json
import math
from fractions import Fraction
from pathlib import Path


ROOT = Path(__file__).resolve().parent
COMPONENTS = ("assembly", "motion", "site", "evidence", "operation")
SOURCE_REPORTS = COMPONENTS + ("outcome",)


def read_json(path: Path):
    return json.loads(path.read_text(encoding="utf-8"))


def sha256(path: Path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def content_sha256(value):
    encoded = json.dumps(value, sort_keys=True, separators=(",", ":"), ensure_ascii=False, allow_nan=False)
    return hashlib.sha256(encoded.encode("utf-8")).hexdigest()


def number(value):
    if isinstance(value, str):
        return float(Fraction(value))
    return float(value)


def integrated_components(report):
    """Accept the orchestration wrapper without silently accepting absent reports."""
    container = report.get("reports", report)
    result = {}
    for name in SOURCE_REPORTS:
        value = container.get(name, container.get(name + "_report"))
        if not isinstance(value, dict):
            raise ValueError("Integrated report lacks the " + name + " report")
        result[name] = value
    expected = report.get("design_sha256")
    if not isinstance(expected, str) or len(expected) != 64:
        raise ValueError("A complete canonical design SHA-256 is required")
    for name, value in result.items():
        if value.get("design_sha256") != expected:
            raise ValueError(name + " belongs to a different or unbound design")
        if report["report_content_sha256"].get(name) != content_sha256(value):
            raise ValueError(name + " content differs from the integrated receipt")
    assembly, motion = result["assembly"], result["motion"]
    if motion.get("geometry_and_mass_source_sha256") != content_sha256(assembly["assembly"]):
        raise ValueError("Motion geometry/mass do not bind to the supplied assembly content")
    if motion.get("loaded_capture_source_sha256") != content_sha256(assembly["capture_reference"]):
        raise ValueError("Loaded capture does not bind to the supplied structural case")
    if result["site"].get("geometry_and_mass_source_sha256") != content_sha256(assembly["assembly"]):
        raise ValueError("Site placement does not bind to the supplied assembly content")
    if result["site"].get("motion_report_sha256") != content_sha256(motion):
        raise ValueError("Site sweep does not bind to the supplied motion report")
    if result["outcome"].get("site_report_sha256") != content_sha256(result["site"]):
        raise ValueError("Outcome record does not bind to this site report")
    nodes = set(assembly["assembly"]["physical_nodes"])
    for frame in motion["frames"]:
        if set(frame["world_nodes_m"]) != nodes:
            raise ValueError("Motion frame changes the physical node identities: " + frame["id"])
        total_mass = sum(Fraction(half["mass_kg"]) for half in frame["halves"])
        if total_mass != Fraction(assembly["assembly"]["mass_ledger"]["bridge_mass_kg"]):
            raise ValueError("Motion frame changes the declared moving mass: " + frame["id"])
    if result["operation"]["public_opening"]["allowed"] is not False:
        raise ValueError("This reference workbench cannot imply physical public opening")
    return result


def scene_payload(report, parts):
    """Normalize only explicit canonical geometry, without a second motion solver."""
    structure, motion, site = parts["assembly"], parts["motion"], parts["site"]
    assembly = structure["assembly"]
    dimensions, design = assembly["dimensions"], assembly["design"]
    mapping = site_mapping(site)
    frames = []
    half_length = Fraction(dimensions["half_length_m"])
    floor_half_width = Fraction(dimensions["floor_width_m"]) / 2
    walking_z = Fraction(dimensions["walking_surface_offset_m"])
    member_radius = Fraction(design["geometry"]["member_envelope_radius_m"])
    local_floor = [[x, y, walking_z] for x, y in (
        (0, -floor_half_width), (half_length, -floor_half_width),
        (half_length, floor_half_width), (0, floor_half_width))]
    support_inside = motion["checks"]["continuous_support_inside"]
    for source in motion["frames"]:
        frame = {key: source[key] for key in ("id", "stage_id", "label", "progress", "world_nodes_m")}
        frame["stage_label"] = source["label"]
        frame["scope"] = "exported pose; continuous conclusions come from receipts"
        frame["centres_of_mass_m"] = {h["id"]: h["center_of_mass_m"] for h in source["halves"]}
        frame["ideal_fixture_axes_m"] = [
            {"half": half["id"], "start": ["0" if half["id"] == "L" else dimensions["span_m"], "0", "0"],
             "end": half["transform"]["translation"]} for half in source["halves"]]
        decks = []
        for half in source["halves"]:
            rotation = [[Fraction(v) for v in row] for row in half["transform"]["rotation"]]
            translation = [Fraction(v) for v in half["transform"]["translation"]]
            corners = [[str(translation[i] + sum(rotation[i][j] * point[j] for j in range(3)))
                        for i in range(3)] for point in local_floor]
            decks.append({"half": half["id"], "corners_m": corners})
        left_transform = source["halves"][0]["transform"]
        angle = abs(math.degrees(math.atan2(number(left_transform["rotation"][1][0]),
                                           number(left_transform["rotation"][0][0]))))
        frame["deck_surfaces"] = decks
        centre_heights = [number(half["center_of_mass_m"][2]) for half in source["halves"]]
        frame["display"] = {
            "clearance_m": min(number(p[2]) for p in source["world_nodes_m"].values()) - float(member_radius),
            "angle_degrees": angle,
            "cg_label": f"{centre_heights[0]:.3f} m (both)" if abs(centre_heights[0] - centre_heights[1]) < 1e-12
                else f"{centre_heights[0]:.2f} / {centre_heights[1]:.2f} m",
            "support_label": "inside declared support bounds" if support_inside else "support condition not established",
            "clearance_status": "model plane only; site unresolved",
            "coupling_label": "nominal positions coincide; no lock evidence"
                if source["stage_id"] == "approach" and Fraction(source["progress"]) == 1
                else "independently held halves",
        }
        frames.append(frame)
    bases = design["base_interfaces"]
    pad = bases["pad_rectangle_m"]
    pivot_y = Fraction(bases["pivot_position_in_base_m"][1])
    supports = [[float(root + Fraction(x)), float(Fraction(y) - pivot_y), 0]
                for root in (Fraction(0), Fraction(dimensions["span_m"]))
                for x in (pad["x_min"], pad["x_max"])
                for y in (pad["y_min"], pad["y_max"])]
    return {
        "frames": frames,
        "members": [{key: m[key] for key in ("id", "half", "start", "end", "role")}
                    for m in assembly["physical_members"]],
        "span_m": dimensions["span_m"],
        "assembly_mass_kg": assembly["mass_ledger"]["bridge_mass_kg"],
        "walking_height_m": dimensions["installed_walking_surface_elevation_m"],
        "truss_depth_m": dimensions["truss_depth_m"],
        "member_envelope_radius_m": str(member_radius),
        "joint_contact_zone_length_per_half_m": design["geometry"]["joint_contact_zone_length_m"],
        "support_points_m": supports,
        "clearance_reference_m": motion["reference_clearance"]["source_reference_m"],
        "road_width_m": mapping["road_width_m"],
        "map_transform": mapping["map_transform"],
        "footprints": mapping["footprints"],
        "building_hit_ids": mapping.get("building_hit_ids", []),
        "direct_building_hit_ids": mapping["direct_building_hit_ids"],
        "sweep_building_hit_ids": mapping["sweep_building_hit_ids"],
        "building_sensitivity_ids": mapping["building_sensitivity_ids"],
        "site_label": mapping["site_label"],
        "site_scenario_id": mapping["scenario_id"],
        "site_architecture": mapping["architecture"],
        "site_architecture_label": mapping["architecture_label"],
        "station_offset_m": mapping["station_offset_m"],
        "site_hazards": mapping["hazards"],
        "site_summary": site["summary"],
        "swept_boxes_model": motion["swept_boxes_model"],
        "scope": motion["scope"],
        "coordinates": motion["coordinate_convention"],
        "source_report": "motion-report.json frames; assembly-report.json physical members and floor dimensions",
        "interface_cards": interface_cards(parts, mapping),
        "integration_findings": [
            ["Walking level, shared depth", f"{number(dimensions['installed_walking_surface_elevation_m']):.5f} m"],
            ["Structure below walking plane", f"{number(dimensions['truss_depth_m']):.4f} m truss + 0.1 m floor"],
            ["Nearby buildings at 3 m sensitivity", str(len(mapping["building_sensitivity_ids"]))],
            ["Real independent field measurements", "0 in the synthetic evidence fixtures"],
            ["Physical opening", "unresolved"],
        ],
    }


def site_mapping(site):
    selected = next(row for row in site["scenarios"] if row["id"] == site["selected_scenario_id"])
    placement = next(row for row in site["placements"] if row["id"] == selected["placement_id"])
    frame = placement["scene_mapping"]
    midpoint, normal, tangent = frame["origin_at_span_midpoint"], frame["cross_road_unit"], frame["along_road_unit"]
    origin = [midpoint[i] - frame["model_midpoint_x_m"] * normal[i] for i in (0, 1)]
    direct = selected["direct_hits"]["building_ids"]
    swept = placement["sweep_possible_hits"]["building_ids"]
    sensitivity = sorted(set(selected["direct_hits"]["building_ids_at_3m_sensitivity"]) |
                         set(placement["sweep_possible_hits"]["building_ids_at_3m_sensitivity"]))
    architecture = {
        "straight-ramp": "straight ramp", "switchback-ramp": "unstacked switchback",
        "elevator-and-stair-allowance": "elevator + stair allowance",
    }[selected["architecture"]]
    return {
        "scenario_id": selected["id"], "architecture": selected["architecture"],
        "architecture_label": architecture, "station_offset_m": placement["station_offset_m"],
        "site_label": placement["label"], "road_width_m": selected["road_width_hypothesis_m"],
        "clearance_reference_m": "5.2832",
        "map_transform": {"origin": origin, "x_axis": normal, "y_axis": tangent},
        "footprints": selected["footprints"],
        "direct_building_hit_ids": direct, "sweep_building_hit_ids": swept,
        "building_sensitivity_ids": sensitivity,
        "building_hit_ids": sorted(set(direct) | set(swept)),
        "summary_rows": [
            ["Selected research access", architecture],
            ["Station offset", f"{placement['station_offset_m']:+g} m"],
            ["Base/access building hits", str(len(direct))],
            ["Base/access road-centreline hits", str(len(selected["direct_hits"]["vehicle_centerline_ids"]))],
            ["Possible swept-bound building hits", str(len(swept))],
            ["Buildings at 3 m sensitivity", str(len(sensitivity))],
        ],
        "hazards": selected["hazards"],
    }


def interface_cards(parts, mapping):
    assembly, motion = parts["assembly"], parts["motion"]
    evidence, operation = parts["evidence"], parts["operation"]
    cases = assembly["cases"]
    stable = sum(row["status"] == "unique_linear_solution" for row in cases)
    mechanisms = [row for row in cases if row["status"] == "mechanism"]
    crowd = next(row for row in cases if row["id"] == "joined_service/dead_plus_symmetric_crowd")
    displacement = abs(number(crowd["scalar_extrema"]["largest_absolute_vertical_displacement"]["signed_value_m"])) * 1000
    capture_mismatch = max(abs(number(value)) for row in motion["loaded_capture"]["connectors"]
                           for value in row["nominal_difference_m"]) * 1000
    lowest_loaded_envelope = min(number(row["lowest_declared_member_surface_m"])
                                 for row in motion["reference_clearance"]["cases"])
    return [
        {"module": "01 / assembly", "title": "One mass ledger, several support states",
         "metric": f"{assembly['assembly']['physical_member_count']} physical bars",
         "description": "The same three-dimensional assembly is used before joining, during controlled restraint transfer and in service load cases.",
         "rows": [["Stable linear cases", f"{stable} / {len(cases)}"],
                  ["Deliberate release controls", f"{len(mechanisms)} mechanisms detected"],
                  ["Symmetric crowd response", f"{displacement:.3f} mm maximum vertical displacement"]],
         "limit": "Response is not a capacity check. Member buckling, connections, floor/rail capacity and dynamic response remain outside this model."},
        {"module": "02 / motion", "title": "Lift, slew, then short slide",
         "metric": f"{len(motion['frames'])} exported poses",
         "description": "A prescribed fixture carries unchanged rigid halves. Reported continuous bounds, not the animation samples, support the scoped conclusions.",
         "rows": [["Support certificates", str(len(motion["support_certificates"]))],
                  ["Full FE base-wrench checks", str(len(motion["assembled_support"]["case_support"]))],
                  ["Continuous release-support paths", str(len(motion["assembled_support"]["release_support"]))],
                  ["Swept interval boxes", str(len(motion["swept_boxes_model"]))],
                  ["Loaded seam mismatch, maximum axis", f"{capture_mismatch:.3f} mm"],
                  ["Lowest loaded envelope, flat model", f"{lowest_loaded_envelope:.4f} m"]],
         "limit": "Half separation excludes the 0.15 m reserved joint-contact zone per half. Actual engagement, omitted hardware, actuator capacities and deployment time remain unresolved."},
        {"module": "03 / site", "title": mapping["site_label"],
         "metric": f"{len(mapping['building_sensitivity_ids'])} buildings at 3 m sensitivity",
         "description": "Actual mapped outlines remain visible at the chosen reference placement. A best-ranked screen is not an approved site.",
         "rows": mapping.get("summary_rows", []),
         "limit": "Source widths are authored. Building heights, survey control, utilities, permissions and a complete accessible route remain unresolved."},
        {"module": "04 / evidence", "title": "A command is not a measurement",
         "metric": f"{len(evidence['cases'])} synthetic cases",
         "description": "One compatible candidate must satisfy all observed channels. Missing, inconsistent and ambiguous evidence retain their distinct outcomes.",
         "rows": [["Physical measurement count", str(evidence.get("physical_measurement_count", 0))],
                  ["Outcome vocabulary", "established / fails / ambiguous / mismatch"],
                  ["Evidence binding", "design + model + source"]],
         "limit": "Synthetic fixtures exercise inference. Physical candidate completeness, calibration and independent sensing are unverified."},
        {"module": "05 / operation", "title": "Model transitions keep their guards",
         "metric": f"{len(operation['reachable_modes'])} reachable model modes",
         "description": "The finite graph consumes bound assembly and motion receipts and replays the selected synthetic evidence case.",
         "rows": [["Declared graph invariant", "holds" if operation["safety"]["declared_graph_invariants_hold"] else "fails"],
                  ["Model service reachable", str(operation["liveness"]["model_service_state_reached"]).lower()],
                  ["Public-use authorization", "not established"]],
         "limit": "Existence of a model path is not guaranteed completion. Physical interlocks, loss-of-power behavior and recovery remain separate obligations."},
    ]


def build_payload():
    report_path = ROOT / "integrated-report.json"
    report = read_json(report_path)
    parts = integrated_components(report)
    atlas_path = ROOT.parent / "v2" / "data" / "atlas-subset.json"
    atlas = read_json(atlas_path)
    scene = scene_payload(report, parts)
    return {
        "designSha256": report["design_sha256"],
        "reports": {
            "evidence": {"design_sha256": report["design_sha256"],
                         "cases": [{"case_id": row["case_id"], "description": row.get("description"),
                                    "result": row["result"]} for row in parts["evidence"]["cases"]]},
            "operation": {key: parts["operation"][key] for key in (
                "design_sha256", "guards", "reachable_modes", "public_opening", "safety", "liveness")},
            "outcome": {key: parts["outcome"][key] for key in (
                "design_sha256", "source_model", "site_results", "summary")},
        },
        "scene": scene,
        "atlas": {key: atlas[key] for key in (
            "origin", "bounds", "counts", "buildings", "parcels", "streets",
            "canopy", "areas", "provenance")},
        "sourceHashes": {
            "integrated-report.json": sha256(report_path),
            "../v2/data/atlas-subset.json": sha256(atlas_path),
        },
    }


HTML = r'''<!doctype html>
<html lang="en"><head><meta charset="utf-8"><meta name="viewport" content="width=device-width,initial-scale=1">
<meta http-equiv="Content-Security-Policy" content="default-src 'none'; script-src 'unsafe-inline'; style-src 'unsafe-inline'; img-src data:; connect-src 'none'; base-uri 'none'; form-action 'none'">
<meta name="description" content="One shared 23 metre reference crossing: actual Tillamook mapping, prescribed deployment, exact assembly checks and explicit evidence limits. Offline engineering research workbench.">
<title>Tillamook · One design, connected evidence</title>
<style>
:root{--ink:#183331;--muted:#61716d;--paper:#f2f2ec;--card:#fffef9;--line:#d3ddd5;--teal:#006d77;--aqua:#abd9d4;--amber:#97551c;--red:#b64c39;--navy:#28557a;--soft:#e9f1ec}*{box-sizing:border-box}body{margin:0;background:var(--paper);color:var(--ink);font:14px/1.5 system-ui,-apple-system,BlinkMacSystemFont,"Segoe UI",sans-serif}button,input,select{font:inherit}button{cursor:pointer;border:1px solid var(--line);background:var(--card);color:var(--ink);border-radius:7px;padding:8px 11px;min-height:38px}button:hover{border-color:var(--teal);background:#eaf4f0}button:focus-visible,a:focus-visible,input:focus-visible,select:focus-visible,canvas:focus-visible{outline:3px solid #578eb5;outline-offset:3px}button[aria-pressed=true]{color:white;background:var(--teal);border-color:var(--teal)}a{color:var(--teal);text-underline-offset:3px}h1,h2,h3,p{margin:0}h1{font-size:clamp(29px,4vw,49px);font-weight:650;letter-spacing:-.045em;line-height:1.1}h2{font-size:24px;line-height:1.2;font-weight:650;letter-spacing:-.025em}h3{font-size:15px;font-weight:700}small{font-size:11px}.wrap{max-width:1560px;padding:0 28px;margin:auto}.topbar{min-height:66px;display:flex;justify-content:space-between;align-items:center;gap:18px;border-bottom:1px solid var(--line)}.brand{font-size:13px;letter-spacing:.12em;font-weight:750}.brand span{letter-spacing:0;color:var(--muted);font-weight:400;margin-left:10px}.nav{display:flex;gap:22px;font-size:12px;flex-wrap:wrap}.nav a{color:var(--ink);text-decoration:none}.hero{padding:28px 0 22px;display:grid;grid-template-columns:1.35fr 1fr;gap:32px;align-items:end}.eyebrow{font-size:10px;letter-spacing:.16em;font-weight:750;color:var(--teal);text-transform:uppercase}.hero h1{margin:8px 0 12px}.hero p{color:var(--muted);font-size:14px;max-width:760px}.scope-note{border-left:3px solid #ae762c;padding:8px 0 8px 18px;color:#67542f;font-size:13px}.scope-note strong{display:block;color:var(--ink);font-size:15px;margin-bottom:5px}.metrics{display:grid;grid-template-columns:repeat(4,1fr);gap:1px;background:var(--line);border:1px solid var(--line);border-radius:10px;overflow:hidden;margin-bottom:20px}.metric{padding:15px 19px;background:var(--card)}.metric b{display:block;font-size:26px;font-weight:650;letter-spacing:-.035em;line-height:1.2}.metric span{font-size:11px;color:var(--muted)}.metric.warn b{color:var(--amber)}.panel{border:1px solid var(--line);border-radius:12px;background:var(--card);overflow:hidden}.panel-head{padding:14px 18px;border-bottom:1px solid var(--line);display:flex;justify-content:space-between;align-items:center;gap:12px}.panel-head p{color:var(--muted);font-size:11px;margin-top:3px}.tag{padding:4px 8px;border:1px solid #d8c49d;color:#8a5922;background:#fff4dc;border-radius:30px;font-size:10px;font-weight:700;white-space:nowrap}.workspace{display:grid;grid-template-columns:minmax(0,1.65fr) minmax(320px,1fr);gap:18px;align-items:start}.scene-shell{position:relative;background:#eaf0eb}.scene-shell canvas{display:block;width:100%;height:470px;touch-action:none;cursor:grab}.canvas-note{position:absolute;top:12px;left:14px;right:14px;pointer-events:none;color:var(--muted);font-size:10px;display:flex;justify-content:space-between;gap:15px}.canvas-note span{padding:5px 8px;background:#fffef8e8;border-radius:5px;border:1px solid #dde4dc}.tools{position:absolute;bottom:13px;right:13px;display:flex;gap:5px}.tools button{padding:6px 9px;font-size:11px;background:#fffef9ed}.legend{display:flex;flex-wrap:wrap;gap:7px 16px;padding:10px 17px;font-size:10px;color:var(--muted);border-top:1px solid var(--line)}.key{display:inline-flex;align-items:center;gap:5px}.swatch{width:12px;height:8px;display:inline-block;border:1px solid #7d9288;background:#cbd2c8;border-radius:2px}.swatch.west{background:var(--teal);border-color:var(--teal)}.swatch.east{background:var(--navy);border-color:var(--navy)}.swatch.obstacle{background:#dca08d;border-color:var(--red)}.swatch.access{background:#adcfc5;border-style:dashed}.swatch.pad{background:#b88b4b;border-color:#9b7038}.controls{padding:15px 18px;border-top:1px solid var(--line)}.phase-top{display:flex;justify-content:space-between;gap:12px;align-items:center;margin-bottom:8px}.phase-title b{font-size:15px;display:block}.phase-title span{color:var(--muted);font-size:11px}input[type=range]{width:100%;accent-color:var(--teal);height:30px;cursor:pointer}.phase-buttons{display:flex;gap:6px;flex-wrap:wrap;margin-top:7px}.phase-buttons button{font-size:10px;padding:5px 8px;min-height:32px;flex:1 1 auto}.helper{font-size:10px;color:var(--muted);margin-top:7px}.right-stack{display:grid;gap:18px}.map-shell{position:relative;background:#eaf0e6}.map-shell canvas{width:100%;height:343px;display:block;touch-action:none;cursor:grab}.toggle-row{display:flex;gap:13px;flex-wrap:wrap;padding:10px 17px;font-size:10px;border-top:1px solid var(--line);color:var(--muted)}.toggle-row label{display:flex;align-items:center;gap:5px;cursor:pointer}.toggle-row input{accent-color:var(--teal)}.phase-status{padding:15px 17px}.status-line{display:flex;justify-content:space-between;gap:15px;padding:7px 0;border-bottom:1px solid #e8ece5;font-size:12px}.status-line:last-child{border-bottom:0}.status-line span:first-child{color:var(--muted)}.status-line b{font-size:11px;text-align:right}.bad{color:var(--red)}.uncertain{color:var(--amber)}.computed{color:var(--navy)}.side-panel{margin-top:18px}.side-panel canvas{width:100%;height:220px;display:block;background:#eef2ec}.strip{display:grid;grid-template-columns:repeat(4,1fr);padding:12px 18px;border-top:1px solid var(--line);gap:18px}.strip b{font-size:18px;font-weight:650;display:block}.strip span{font-size:10px;color:var(--muted)}section.content{padding-top:30px}.section-head{display:flex;justify-content:space-between;gap:24px;align-items:end;margin-bottom:15px}.section-head p{color:var(--muted);font-size:12px;max-width:610px}.interface-grid{display:grid;grid-template-columns:repeat(3,minmax(0,1fr));gap:15px}.interface{padding:18px}.interface .eyebrow{font-size:9px;margin-bottom:6px}.interface p{font-size:12px;color:var(--muted);margin-top:8px}.interface .figure{font-size:28px;letter-spacing:-.035em;font-weight:650;margin:10px 0 0}.mini{display:grid;gap:3px;margin-top:12px}.mini div{display:flex;gap:14px;justify-content:space-between;font-size:11px;padding:4px 0;border-bottom:1px solid #e9ece7}.mini span{color:var(--muted)}.mini b{font-weight:600;text-align:right}.callout{margin-top:12px;padding:10px 12px;background:#f6eedf;border-left:3px solid #b27a35;font-size:11px;color:#6f5c3a}.evidence-grid{display:grid;grid-template-columns:1fr 1.1fr;gap:18px}.table-wrap{overflow:auto}.table{border-collapse:collapse;width:100%;font-size:11px}.table th{text-align:left;color:var(--muted);font-weight:650;border-bottom:1px solid var(--line);padding:10px 12px;white-space:nowrap}.table td{padding:10px 12px;border-bottom:1px solid #e9ede6;vertical-align:top}.table tr:last-child td{border-bottom:none}.badge{display:inline-block;font-size:10px;padding:3px 6px;border-radius:4px;line-height:1.25;background:#e7eef5;color:#2f5877}.badge.ambiguous,.badge.unresolved,.badge.unknown{background:#fcf0d6;color:#8a5b21}.badge.condition_fails,.badge.failed,.badge.model_mismatch{background:#f8e3dc;color:#a44231}.badge.condition_established,.badge.satisfied{background:#e3efed;color:#276762}.journey{padding:18px}.journey p{font-size:12px;color:var(--muted);margin:8px 0}.formula{padding:13px;background:#edf2ee;border:1px solid #dce4dc;border-radius:7px;font:13px/1.7 ui-monospace,SFMono-Regular,Consolas,monospace;overflow:auto;margin:13px 0}.mode-list{display:flex;gap:7px;flex-wrap:wrap;padding:16px}.mode-chip{border:1px solid #c8d6cf;padding:6px 9px;border-radius:6px;font-size:11px}.mode-chip.blocked{color:#955320;border-style:dashed;background:#fcf2df}.open-boundary{margin:0 16px 16px;padding:11px 12px;background:#fcf0dd;border:1px solid #e2caa2;border-radius:7px;font-size:12px;color:#79531f}.reading{display:grid;grid-template-columns:repeat(4,1fr);gap:12px}.reading a{padding:16px;border:1px solid var(--line);border-radius:9px;text-decoration:none;background:var(--card);font-size:12px}.reading b{display:block;margin-bottom:5px}.reading span{color:var(--muted);font-size:10px}.provenance{margin-top:16px;padding:17px;font-size:11px;color:var(--muted)}.provenance p{margin:6px 0}.hash{font:10px/1.6 ui-monospace,Consolas,monospace;word-break:break-all}.provenance details{margin-top:10px}.provenance summary{cursor:pointer}.provenance pre{white-space:pre-wrap;word-break:break-word;font-size:10px;max-height:220px;overflow:auto;margin-top:9px}.footer{padding:24px 0 30px;font-size:10px;color:var(--muted)}.fail{padding:18px;background:#f7ddd5;color:#8b3529;margin:20px 0;border:1px solid #c76550;border-radius:8px}canvas:focus{outline-offset:-4px}.muted{color:var(--muted)}.sr-only{position:absolute;width:1px;height:1px;padding:0;margin:-1px;overflow:hidden;clip:rect(0,0,0,0);white-space:nowrap;border:0}@media(max-width:1100px){.workspace{grid-template-columns:minmax(0,1.4fr) minmax(300px,1fr)}.scene-shell canvas{height:430px}.interface-grid{grid-template-columns:repeat(2,minmax(0,1fr))}.hero{gap:24px}.reading{grid-template-columns:repeat(2,1fr)}}@media(max-width:800px){.wrap{padding:0 16px}.hero{grid-template-columns:1fr;gap:18px;padding-top:23px}.metrics{grid-template-columns:repeat(2,1fr)}.metric{padding:13px 15px}.workspace{grid-template-columns:1fr}.scene-shell canvas{height:380px}.right-stack{grid-template-columns:1fr}.map-shell canvas{height:330px}.interface-grid,.evidence-grid{grid-template-columns:1fr}.strip{grid-template-columns:repeat(2,1fr);gap:10px}.section-head{display:block}.section-head p{margin-top:8px}.topbar{min-height:70px}.nav{gap:12px}.brand span{display:none}.canvas-note{font-size:9px}.side-panel canvas{height:210px}.reading{grid-template-columns:repeat(2,1fr)}.table td,.table th{padding:8px}.phase-buttons button{font-size:9px}.scope-note{font-size:12px}.panel-head{padding:13px 15px}.panel-head h3{font-size:14px}.tag{font-size:9px}.phase-top{gap:6px}}@media(prefers-reduced-motion:reduce){*{scroll-behavior:auto!important}}
</style></head>
<body><div class="wrap"><header class="topbar"><div class="brand">CREST–ROOT <span>Research workbench · v3</span></div><nav class="nav" aria-label="Page sections"><a href="#model">Model</a><a href="#interfaces">Interfaces</a><a href="#evidence">Evidence</a><a href="#sources">Sources</a></nav></header>
<main><section class="hero"><div><div class="eyebrow">Tillamook / one shared reference design</div><h1>Connect the design.<br>Keep the evidence visible.</h1><p>A 23 m deck-truss reference links the same assembly, prescribed lift/slew/short-slide movement, mapped site, evidence tests and operating states. Move through its phases to see what has been computed and what still needs physical evidence.</p></div><div class="scope-note"><strong>Integrated software demonstration.</strong>Geometry, masses and loads are declared reference inputs. Actual mapped obstacles stay visible. A computed model result does not authorize deployment or public use.</div></section>
<div class="metrics" aria-label="Shared reference summary"><div class="metric"><b id="spanMetric">23 m</b><span>shared reference base-centre span</span></div><div class="metric"><b id="massMetric">—</b><span>bridge assembly mass · bases excluded</span></div><div class="metric"><b id="phaseMetric">—</b><span>reference walking surface height</span></div><div class="metric warn"><b id="openingMetric">Unresolved</b><span>physical public-use requirements</span></div></div>
<section id="model" class="workspace" aria-label="Connected reference model"><div><article class="panel"><header class="panel-head"><div><h3>Shared assembly · prescribed deployment</h3><p id="sceneDescription">Rigid motion from the canonical report; member geometry stays bound to the design.</p></div><span class="tag">Reference mechanism</span></header><div class="scene-shell"><canvas id="scene" tabindex="0" aria-label="Three dimensional reference bridge. Drag to orbit; arrow keys rotate; plus and minus zoom."></canvas><div class="canvas-note"><span>Member centrelines + floor reference · schematic ground</span><span id="scenePhase">—</span></div><div class="tools"><button id="sceneReset" aria-label="Reset three dimensional view">Reset view</button><button id="sceneZoomOut" aria-label="Zoom out three dimensional view">−</button><button id="sceneZoomIn" aria-label="Zoom in three dimensional view">+</button></div></div><div class="legend"><span class="key"><i class="swatch west"></i>Half A</span><span class="key"><i class="swatch east"></i>Half B</span><span class="key"><i class="swatch pad"></i>Support locations</span><span class="key"><i class="swatch obstacle"></i>Mapped obstacle</span></div><div class="controls"><div class="phase-top"><div class="phase-title"><b id="phaseLabel">—</b><span id="phaseDetail">—</span></div><button id="play" aria-label="Play reference phases" aria-pressed="false">Play phases</button></div><label class="sr-only" for="phase">Reference phase sample</label><input id="phase" type="range" min="0" max="0" value="0" step="1"><div id="phaseButtons" class="phase-buttons"></div><p class="helper">Playback steps through exported reference poses. Playback speed is not a deployment-time prediction. The slider is keyboard accessible.</p></div></article>
<article class="panel side-panel"><header class="panel-head"><div><h3>Side elevation · same phase</h3><p>Same coordinates and levels. This projection collapses the along-road axis; overlap here is not a collision result.</p></div><span class="tag">Clearance is scoped</span></header><canvas id="side" aria-label="Side elevation of the current canonical reference pose"></canvas><div class="strip"><div><b id="clearanceMetric">—</b><span>lowest member envelope above model plane</span></div><div><b id="angleMetric">—</b><span>rotation remaining to final alignment</span></div><div><b id="cgMetric">—</b><span>moving mass-centre elevation</span></div><div><b id="supportMetric">—</b><span>reported support condition</span></div></div></article></div>
<div class="right-stack"><article class="panel"><header class="panel-head"><div><h3 id="siteTitle">Actual Tillamook map</h3><p>Fixed reference placement · no hidden relocation</p></div><span class="tag">Mapped source geometry</span></header><div class="map-shell"><canvas id="map" tabindex="0" aria-label="Tillamook reference site map. Drag or use arrow keys to pan; plus and minus zoom."></canvas><div class="canvas-note"><span id="mapNote">Source widths are authored, not measured</span></div><div class="tools"><button id="mapFit">Fit site</button><button id="mapDistrict">District</button><button id="mapZoomOut" aria-label="Zoom out map">−</button><button id="mapZoomIn" aria-label="Zoom in map">+</button></div></div><div class="toggle-row"><label><input id="parcels" type="checkbox" checked>Parcels</label><label><input id="canopy" type="checkbox">Canopy proxies</label><label><input id="widths" type="checkbox">Authored road widths</label><label><input id="access" type="checkbox" checked>Access allowances</label><label><input id="swept" type="checkbox">Continuous swept boxes</label></div><div class="legend"><span class="key"><i class="swatch"></i>395 source buildings</span><span class="key"><i class="swatch obstacle"></i>Red: footprint · ochre: sweep / 3 m flag</span><span class="key"><i class="swatch access"></i>Access allowance</span></div></article><article class="panel"><header class="panel-head"><div><h3>At this phase</h3><p>Report values and operation guards remain distinct.</p></div></header><div id="phaseStatus" class="phase-status"></div></article></div></section>
<section id="interfaces" class="content"><div class="section-head"><h2>Five modules. One design identity.</h2><p>These results share a design hash. That checks interface consistency; it does not establish that a parameter matches the physical world.</p></div><div id="interfaceCards" class="interface-grid"></div></section>
<section id="evidence" class="content"><div class="section-head"><h2>Evidence controls what can be concluded.</h2><p>Missing measurements, inconsistent models and ambiguous states are explicit outcomes. No unknown physical gate becomes a pass.</p></div><div class="evidence-grid"><article class="panel"><header class="panel-head"><div><h3>Evidence cases</h3><p>Synthetic fixtures exercise the inference contract.</p></div><span class="tag">Not sensor trials</span></header><div class="table-wrap"><table class="table"><thead><tr><th>Case</th><th>Outcome</th><th>Meaning</th></tr></thead><tbody id="evidenceRows"></tbody></table></div></article><article class="panel"><header class="panel-head"><div><h3>Operating-state contract</h3><p>Reachable model modes and blocked transitions.</p></div></header><div id="operationModes" class="mode-list"></div><div id="openingBoundary" class="open-boundary"></div><div class="table-wrap"><table class="table"><thead><tr><th>Guard</th><th>State</th><th>Evidence boundary</th></tr></thead><tbody id="guardRows"></tbody></table></div></article></div></section>
<section class="content"><div class="section-head"><h2>The complete journey is the outcome.</h2><p>The recovered Safe Roads objective connects the mechanism to usable access, uptake and local exposure. Geometry alone supplies none of those effects.</p></div><div class="evidence-grid"><article class="panel journey"><h3>Keep the counterfactual explicit</h3><div class="formula">H₀ − H₁ = N₀r₀ − N₁[au·rᵦ + (1 − au)r𝒻] − D − S</div><p>N₀ and N₁ are person-crossings. Availability a is weighted by crossing demand; uptake u is conditional on accessible availability. Baseline, facility and fallback risks remain separate. D and S cover nonoverlapping deployment and displacement harm.</p><p>The reference bridge does not supply measured local risk, demand or benefit. Input ranges describe scenarios; they are not confidence intervals or a lives-saved forecast.</p></article><article class="panel journey"><h3>Field outcome inputs remain unknown</h3><p id="outcomeSummary"></p><div class="table-wrap"><table class="table"><thead><tr><th>Mapped site</th><th>Deaths</th><th>Serious injuries</th></tr></thead><tbody id="outcomeRows"></tbody></table></div><p id="outcomeMissing"></p><p><a href="outcome-report.json">Read the bound outcome report</a></p><div id="integrationFindings" class="mini"></div><div class="callout">No site benefit is estimated. The next physical test must supply evidence for the ideal actuator, lock, support and access interfaces.</div></article></div></section>
<section id="sources" class="content"><div class="section-head"><h2>Inspect the assumptions and receipts.</h2><p>The page embeds its report and map inputs and makes no network requests. Supporting local files are available in the research packet.</p></div><div class="reading"><a href="integration-review.md"><b>Integration review</b><span>What now connects; physical work still required.</span></a><a href="integrated-report.json"><b>Integrated report</b><span>All five component results and their shared identity.</span></a><a href="../v2/oregon-design-basis.md"><b>Oregon design basis</b><span>Source applicability, clearance and access scope.</span></a><a href="../v2/design-review.md"><b>Earlier site study</b><span>648 hypotheses and the original evidence limits.</span></a></div><article class="panel provenance"><p><b>Design SHA-256</b></p><p id="designHash" class="hash"></p><p id="sourceSummary"></p><p>Mapped outlines retain their source geometry. Parcel boundaries do not grant access. Authored road widths are display proxies. Source coordinate encoding is not survey accuracy; vertical datum and source epochs require engineering confirmation.</p><p>Reference nodes, rational results and motion parameters are converted to browser floating point for drawing. The underlying reports retain their numerical meaning. No approved site, load rating, deployment speed or public opening is established.</p><details><summary>Embedded input hashes</summary><pre id="sourceHashes"></pre></details><details><summary>Canonical scene export</summary><pre id="sceneReceipt"></pre></details></article></section></main>
<footer class="footer">Crest–Root / Tillamook integration study · © OpenStreetMap contributors · DOGAMI · ODOT · All physical readiness statements remain scoped to the evidence supplied.</footer></div>
<script id="payload" type="application/json">__DATA__</script>
<script>
'use strict';
const D=JSON.parse(document.getElementById('payload').textContent),A=D.atlas,R=D.reports,S=D.scene;
const $=id=>document.getElementById(id),esc=x=>String(x??'').replace(/[&<>"']/g,c=>({'&':'&amp;','<':'&lt;','>':'&gt;','"':'&quot;',"'":'&#39;'}[c]));
const q=v=>typeof v==='string'&&v.includes('/')?v.split('/').map(Number).reduce((a,b)=>a/b):Number(v);
const nf=(v,d=2)=>Number.isFinite(q(v))?q(v).toLocaleString('en-US',{maximumFractionDigits:d,minimumFractionDigits:0}):'unreported';
const human=x=>String(x??'unreported').replace(/_/g,' '),arr=x=>Array.isArray(x)?x:[],fmt=(v,u='',d=2)=>v==null?'unreported':nf(v,d)+(u?' '+u:'');
let frameIndex=0,playing=false,timer=null,orbit={yaw:-.70,pitch:.52,zoom:1},sceneDrag=null,mapDrag=null;
const mapView={cx:0,cy:0,scale:3},frames=S.frames,scene=$('scene'),map=$('map'),side=$('side');
const activeFrame=()=>frames[frameIndex];
function canvasSize(c){const rect=c.getBoundingClientRect(),w=rect.width||c.clientWidth||600,h=rect.height||c.clientHeight||300,d=Math.min(window.devicePixelRatio||1,2);if(c.width!==Math.round(w*d)||c.height!==Math.round(h*d)){c.width=Math.round(w*d);c.height=Math.round(h*d);}const ctx=c.getContext('2d');ctx.setTransform(d,0,0,d,0,0);return {ctx,w,h};}
function path(ctx,points,project,close=true){if(!points.length)return;ctx.beginPath();points.forEach((p,i)=>{let z=project(p);i?ctx.lineTo(z[0],z[1]):ctx.moveTo(z[0],z[1]);});if(close)ctx.closePath();}
function strokeLine(ctx,a,b,color,width=1,dash=[]){ctx.beginPath();ctx.strokeStyle=color;ctx.lineWidth=width;ctx.setLineDash(dash);ctx.moveTo(...a);ctx.lineTo(...b);ctx.stroke();ctx.setLineDash([]);}
function nodeMap(){return activeFrame().world_nodes_m||activeFrame().nodes||{};}
function xyz(value){return arr(value).map(q);}
function frameNodes(){const nodes=nodeMap();return Object.fromEntries(Object.entries(nodes).map(([k,v])=>[k,xyz(v)]));}
function members(){return S.members||[];}
function memberEnds(m){return [m.start??m.a??m.nodes?.[0],m.end??m.b??m.nodes?.[1]];}
function halfColor(m){return /^(r|b|right|east|2|half_b)$/i.test(String(m.half??m.half_id??m.side??''))?'#28557a':'#006d77';}
function mapPoint(p){const t=S.map_transform||{},o=t.origin||[0,0],x=t.x_axis||[1,0],y=t.y_axis||[0,1];return [o[0]+p[0]*x[0]+p[1]*y[0],o[1]+p[0]*x[1]+p[1]*y[1]];}
function localPoint(p){const t=S.map_transform||{},o=t.origin||[0,0],x=t.x_axis||[1,0],y=t.y_axis||[0,1],dx=p[0]-o[0],dy=p[1]-o[1];return [dx*x[0]+dy*x[1],dx*y[0]+dy*y[1],0];}
function mapXY(p,w=map.clientWidth,h=map.clientHeight){return [w/2+(p[0]-mapView.cx)*mapView.scale,h/2-(p[1]-mapView.cy)*mapView.scale];}
function mapFit(){const pts=(S.footprints||[]).flatMap(f=>f.ring||[]);if(!pts.length)pts.push(mapPoint([-25,-25,0]),mapPoint([25,25,0]));const min=[Math.min(...pts.map(p=>p[0]))-15,Math.min(...pts.map(p=>p[1]))-15],max=[Math.max(...pts.map(p=>p[0]))+15,Math.max(...pts.map(p=>p[1]))+15];mapView.cx=(min[0]+max[0])/2;mapView.cy=(min[1]+max[1])/2;mapView.scale=Math.min((map.clientWidth||400)/(max[0]-min[0]),(map.clientHeight||343)/(max[1]-min[1]));drawMap();}
function mapZoom(f){mapView.scale=Math.min(25,Math.max(.15,mapView.scale*f));drawMap();}
function mapDistrict(){const b=A.bounds;mapView.cx=(b[0]+b[2])/2;mapView.cy=(b[1]+b[3])/2;mapView.scale=Math.min((map.clientWidth||400)/(b[2]-b[0]),(map.clientHeight||343)/(b[3]-b[1]));drawMap();}
function sourceBuildingHits(){return new Set(S.building_hit_ids||[]);}
function drawMap(){const {ctx,w,h}=canvasSize(map),p=v=>mapXY(v,w,h);ctx.clearRect(0,0,w,h);ctx.fillStyle='#edf1e8';ctx.fillRect(0,0,w,h);for(const area of A.areas){path(ctx,area.p,p);ctx.fillStyle=area.kind==='parking'?'#e2e6dc':'#e5eade';ctx.fill();}if($('parcels').checked){ctx.lineWidth=.7;ctx.strokeStyle='#c1b895';for(const parcel of A.parcels){for(const ring of parcel.rings){path(ctx,ring,p);ctx.stroke();}}}for(const road of A.streets){path(ctx,road.p,p,false);ctx.lineWidth=$('widths').checked?Math.max(1,q(road.width)*mapView.scale):1.1;ctx.strokeStyle=$('widths').checked?'#d2d8ce':'#87968b';ctx.stroke();if($('widths').checked){path(ctx,road.p,p,false);ctx.strokeStyle='#96a197';ctx.lineWidth=.6;ctx.stroke();}}
const hits=sourceBuildingHits(),direct=new Set(S.direct_building_hit_ids||[]),sensitive=new Set(S.building_sensitivity_ids||[]);for(const b of A.buildings){ctx.beginPath();for(const ring of b.rings){ring.forEach((v,i)=>{const xy=p(v);i?ctx.lineTo(...xy):ctx.moveTo(...xy);});ctx.closePath();}ctx.fillStyle=direct.has(b.id)?'#dfa18e':hits.has(b.id)?'#ddbf86':sensitive.has(b.id)?'#daceaf':'#c8cec3';ctx.strokeStyle=direct.has(b.id)?'#b64c39':hits.has(b.id)||sensitive.has(b.id)?'#a27632':'#98a38f';ctx.lineWidth=hits.has(b.id)||sensitive.has(b.id)?1.3:.7;ctx.setLineDash(sensitive.has(b.id)&&!hits.has(b.id)?[3,2]:[]);ctx.fill('evenodd');ctx.stroke();ctx.setLineDash([]);}if($('canopy').checked){ctx.fillStyle='#7a9d6b33';ctx.strokeStyle='#76936d88';for(const t of A.canopy){let xy=p(t.p);ctx.beginPath();ctx.arc(xy[0],xy[1],Math.max(1,t.radius_estimate_m*mapView.scale),0,2*Math.PI);ctx.fill();ctx.stroke();}}
if($('swept').checked){for(const box of S.swept_boxes_model||[]){if(box.stage_id!==activeFrame().stage_id)continue;const a=box.minimum_m.map(q),b=box.maximum_m.map(q);path(ctx,[[a[0],a[1],0],[b[0],a[1],0],[b[0],b[1],0],[a[0],b[1],0]].map(mapPoint),p);ctx.fillStyle='#70a3c010';ctx.strokeStyle='#5c91a447';ctx.lineWidth=.7;ctx.fill();ctx.stroke();}}for(const f of S.footprints||[]){if(f.kind==='access'&&!$('access').checked)continue;path(ctx,f.ring,p);ctx.fillStyle=f.kind==='base'?'#28557ab0':'#abd9d481';ctx.strokeStyle=f.kind==='base'?'#224c6b':'#3f8d81';ctx.lineWidth=1.4;ctx.setLineDash(f.kind==='access'?[4,3]:[]);ctx.fill();ctx.stroke();ctx.setLineDash([]);}const nodes=frameNodes();for(const m of members()){const [a,b]=memberEnds(m);if(!nodes[a]||!nodes[b])continue;strokeLine(ctx,p(mapPoint(nodes[a])),p(mapPoint(nodes[b])),halfColor(m),Math.max(1.3,Math.min(2.5,mapView.scale*.35)));}const o=p((S.map_transform||{}).origin||[0,0]);ctx.fillStyle='#2b5147';ctx.beginPath();ctx.arc(...o,2.6,0,Math.PI*2);ctx.fill();drawScale(ctx,w,h,mapView.scale);ctx.save();ctx.translate(w-25,42);ctx.fillStyle='#375d4f';ctx.font='bold 11px system-ui';ctx.textAlign='center';ctx.fillText('N',0,-9);ctx.beginPath();ctx.moveTo(0,-3);ctx.lineTo(-4,10);ctx.lineTo(0,7);ctx.lineTo(4,10);ctx.closePath();ctx.fill();ctx.restore();}
function drawScale(ctx,w,h,scale){let length=scale>6?5:scale>2?10:scale>.8?25:100;ctx.fillStyle='#fffef8e8';ctx.fillRect(12,h-40,length*scale+19,28);strokeLine(ctx,[21,h-22],[21+length*scale,h-22],'#375649',2);ctx.fillStyle='#375649';ctx.font='10px system-ui';ctx.textAlign='left';ctx.fillText(length+' m',21,h-28);}
function project3D(p,w,h){const y=orbit.yaw,tilt=orbit.pitch,c=Math.cos(y),s=Math.sin(y),dx=p[0]-q(S.span_m??23)/2,x=dx*c-p[1]*s,depth=dx*s+p[1]*c,scale=Math.min(w/42,h/31)*orbit.zoom;return [w/2+x*scale,h*.64+(depth*Math.sin(tilt)-p[2]*Math.cos(tilt))*scale,depth*Math.cos(tilt)+p[2]*Math.sin(tilt)];}
function drawScene(){const {ctx,w,h}=canvasSize(scene),p=v=>project3D(v,w,h);ctx.clearRect(0,0,w,h);ctx.fillStyle='#eaf0e9';ctx.fillRect(0,0,w,h);const r=q(S.road_width_m??12)/2,mid=q(S.span_m??23)/2;path(ctx,[[mid-r,-22,0],[mid+r,-22,0],[mid+r,22,0],[mid-r,22,0]],p);ctx.fillStyle='#d5ddd3';ctx.fill();strokeLine(ctx,p([mid,-22,0]),p([mid,22,0]),'#a9b5a6',1,[5,5]);for(let g=-20;g<=20;g+=5){strokeLine(ctx,p([g+mid,-20,0]),p([g+mid,20,0]),'#d4dfd388',.6);strokeLine(ctx,p([mid-20,g,0]),p([mid+20,g,0]),'#d4dfd388',.6);}const hits=sourceBuildingHits(),direct=new Set(S.direct_building_hit_ids||[]),sensitive=new Set(S.building_sensitivity_ids||[]);for(const b of A.buildings){const rings=b.rings.map(r=>r.map(localPoint));const flat=rings.flat();if(!flat.length||Math.max(...flat.map(v=>v[0]))<mid-32||Math.min(...flat.map(v=>v[0]))>mid+32||Math.max(...flat.map(v=>v[1]))<-25||Math.min(...flat.map(v=>v[1]))>25)continue;ctx.beginPath();for(const ring of rings){ring.forEach((v,i)=>{const xy=p(v);i?ctx.lineTo(xy[0],xy[1]):ctx.moveTo(xy[0],xy[1]);});ctx.closePath();}ctx.fillStyle=hits.has(b.id)?'#d99b8766':'#c0cbbb88';ctx.strokeStyle=hits.has(b.id)?'#b87662':'#a3b39b';ctx.lineWidth=.7;ctx.fill('evenodd');ctx.stroke();}for(const f of S.footprints||[]){if(f.kind==='access'&&!$('access').checked)continue;path(ctx,f.ring.map(localPoint),p);ctx.fillStyle=f.kind==='base'?'#375b7777':'#8abcb033';ctx.strokeStyle=f.kind==='base'?'#45697c':'#74a096';ctx.lineWidth=1;ctx.setLineDash(f.kind==='access'?[4,3]:[]);ctx.fill();ctx.stroke();ctx.setLineDash([]);}for(const v of S.support_points_m||[]){const xy=p(xyz(v));ctx.fillStyle='#a47a42';ctx.fillRect(xy[0]-3,xy[1]-3,6,6);}for(const axis of activeFrame().ideal_fixture_axes_m||[]){strokeLine(ctx,p(xyz(axis.start)),p(xyz(axis.end)),'#647f79',1.3,[4,4]);}const nodes=frameNodes();const edges=members().map(m=>{const [a,b]=memberEnds(m);return nodes[a]&&nodes[b]?{m,a:nodes[a],b:nodes[b]}:null;}).filter(Boolean).sort((x,y)=>(p(x.a)[2]+p(x.b)[2])-(p(y.a)[2]+p(y.b)[2]));for(const m of edges){let a=p(m.a),b=p(m.b);strokeLine(ctx,[a[0],a[1]],[b[0],b[1]],halfColor(m.m),String(m.m.role||'').includes('diagonal')?1.2:2.2);}for(const deck of activeFrame().deck_surfaces||[]){path(ctx,deck.corners_m.map(xyz),p);ctx.fillStyle=halfColor(deck)+'55';ctx.strokeStyle=halfColor(deck);ctx.lineWidth=1.1;ctx.fill();ctx.stroke();}for(const v of Object.values(nodes)){const xy=p(v);ctx.fillStyle='#164944';ctx.beginPath();ctx.arc(xy[0],xy[1],1.7,0,2*Math.PI);ctx.fill();}for(const [id,v] of Object.entries(activeFrame().centres_of_mass_m||{})){const xy=p(xyz(v));ctx.fillStyle='#d49636';ctx.beginPath();ctx.arc(xy[0],xy[1],4,0,2*Math.PI);ctx.fill();ctx.strokeStyle='#fffaf0';ctx.lineWidth=1.3;ctx.stroke();}ctx.fillStyle='#648072';ctx.font='10px system-ui';ctx.fillText('Dashed axes: ideal fixture constraints, not lifting hardware. Ground plane is schematic.',16,h-16);}
function drawSide(){const {ctx,w,h}=canvasSize(side),nodes=frameNodes(),pts=Object.values(nodes),baseExtent=q(S.span_m??23)/2,top=Math.max(10,...pts.map(p=>p[2]))+1,scale=Math.min((w-100)/(2*baseExtent+12),(h-42)/top),bottom=h-25,p=v=>[w/2+(v[0]-baseExtent)*scale,bottom-v[2]*scale];ctx.clearRect(0,0,w,h);ctx.fillStyle='#eef2ec';ctx.fillRect(0,0,w,h);const r=q(S.road_width_m??12)/2;path(ctx,[[baseExtent-r,0,0],[baseExtent+r,0,0],[baseExtent+r,0,-.3],[baseExtent-r,0,-.3]],p);ctx.fillStyle='#cbd5c8';ctx.fill();strokeLine(ctx,[15,bottom],[w-15,bottom],'#9db19c',1.2);const clear=S.clearance_reference_m;if(clear!=null){strokeLine(ctx,p([baseExtent-r,0,q(clear)]),p([baseExtent+r,0,q(clear)]),'#ae7a38',1,[4,3]);ctx.fillStyle='#8e672e';ctx.font='10px system-ui';ctx.textAlign='center';ctx.fillText(nf(clear,4)+' m clearance reference',w/2,p([0,0,q(clear)])[1]-6);}for(const m of members()){const [a,b]=memberEnds(m);if(nodes[a]&&nodes[b])strokeLine(ctx,p(nodes[a]),p(nodes[b]),halfColor(m),1.5);}for(const deck of activeFrame().deck_surfaces||[]){const corners=deck.corners_m.map(xyz);strokeLine(ctx,p(corners[0]),p(corners[1]),halfColor(deck),3);}for(const v of S.support_points_m||[]){const xy=p(xyz(v));ctx.fillStyle='#967044';ctx.beginPath();ctx.moveTo(xy[0],xy[1]);ctx.lineTo(xy[0]-4,xy[1]+8);ctx.lineTo(xy[0]+4,xy[1]+8);ctx.closePath();ctx.fill();}strokeLine(ctx,p([0,0,-.8]),p([2*baseExtent,0,-.8]),'#78917d',1);ctx.font='10px system-ui';ctx.fillStyle='#567260';ctx.textAlign='center';ctx.fillText(nf(S.span_m??23)+' m base-centre reference',w/2,h-4);}
function badge(status){return `<span class="badge ${esc(status)}">${esc(human(status))}</span>`;}
function statusRow(label,value,cls=''){return `<div class="status-line"><span>${esc(label)}</span><b class="${esc(cls)}">${esc(value)}</b></div>`;}
function selectedPhaseInfo(){return activeFrame().display||{};}
function updateFrame(){const f=activeFrame(),d=selectedPhaseInfo();$('phase').value=String(frameIndex);$('phaseLabel').textContent=f.label||f.stage_label||human(f.stage_id||f.id);$('scenePhase').textContent=$('phaseLabel').textContent;$('phaseDetail').textContent=`Sample ${frameIndex+1} of ${frames.length} · ${f.scope||'prescribed rigid pose'}`;$('clearanceMetric').textContent=fmt(d.clearance_m??f.clearance_m,'m',3);$('angleMetric').textContent=fmt(d.angle_degrees??f.angle_degrees,'°',1);$('cgMetric').textContent=d.cg_label||'shared mass model';$('supportMetric').textContent=d.support_label||human(f.support_state||'unreported');$('phaseStatus').innerHTML=statusRow('Phase',f.label||f.stage_label||human(f.stage_id||f.id))+statusRow('Support model',d.support_label||human(f.support_state||'unreported'))+statusRow('Clearance screen',d.clearance_status||human(f.clearance_status||'unreported'),d.clearance_status==='failed'?'bad':'uncertain')+statusRow('Coupling model',d.coupling_label||human(f.coupling_state||'unreported'))+statusRow('Physical use', 'Not authorized by this model','uncertain');document.querySelectorAll('[data-phase-index]').forEach(b=>b.setAttribute('aria-pressed',String(Number(b.dataset.phaseIndex)===frameIndex)));drawScene();drawSide();drawMap();window.__workbenchState={designSha256:D.designSha256,frameIndex,frameCount:frames.length,stageId:f.stage_id||f.id,memberCount:members().length,nodeCount:Object.keys(nodeMap()).length};}
function stopPlay(){playing=false;if(timer)clearInterval(timer);timer=null;$('play').textContent='Play phases';$('play').setAttribute('aria-pressed','false');}
function setFrame(i){frameIndex=Math.max(0,Math.min(frames.length-1,Math.round(i)));updateFrame();}
function play(){if(playing){stopPlay();return;}playing=true;$('play').textContent='Pause';$('play').setAttribute('aria-pressed','true');timer=setInterval(()=>{if(frameIndex>=frames.length-1){stopPlay();return;}setFrame(frameIndex+1);},Math.max(80,8000/frames.length));}
function reportText(value){if(value==null)return'unreported';if(typeof value==='string'||typeof value==='number'||typeof value==='boolean')return human(value);return JSON.stringify(value);}
function miniRows(rows){return rows.map(([a,b])=>`<div><span>${esc(a)}</span><b>${esc(b)}</b></div>`).join('');}
function renderInterfaces(){const cards=S.interface_cards||[];$('interfaceCards').innerHTML=cards.map(c=>`<article class="panel interface"><div class="eyebrow">${esc(c.module)}</div><h3>${esc(c.title)}</h3><div class="figure">${esc(c.metric||'Shared input')}</div><p>${esc(c.description)}</p><div class="mini">${miniRows(c.rows||[])}</div><div class="callout">${esc(c.limit)}</div></article>`).join('');}
function renderEvidence(){const cases=arr(R.evidence.cases);$('evidenceRows').innerHTML=cases.map(c=>{const r=c.result||{},text=arr(r.reasons).join('; ')||c.description||r.scope||'No explanation supplied';return `<tr><td>${esc(c.case_id)}</td><td>${badge(r.outcome||'unknown')}</td><td>${esc(text)}${arr(r.missing_measurement_ids).length?'<br><span class="muted">Missing: '+esc(r.missing_measurement_ids.join(', '))+'</span>':''}</td></tr>`;}).join('');const guards=R.operation.guards||{};$('guardRows').innerHTML=Object.entries(guards).map(([name,g])=>`<tr><td>${esc(human(name))}</td><td>${badge(g.status||'unresolved')}</td><td>${esc(reportText(g.evidence??g.reason??g.scope))}</td></tr>`).join('');const modes=arr(R.operation.reachable_modes);$('operationModes').innerHTML=modes.map(m=>`<span class="mode-chip">${esc(human(m))}</span>`).join('')+'<span class="mode-chip blocked">Physical public opening unresolved</span>';const opening=R.operation.public_opening||{};$('openingBoundary').textContent='Public opening is not authorized. '+(arr(opening.unresolved_gates).length?'Unresolved gates: '+opening.unresolved_gates.map(human).join('; ')+'.':'Model guards do not establish physical approval.');}
function renderOutcomes(){const sites=R.outcome.site_results;let total=0,unknown=0;const missing=new Set();$('outcomeRows').innerHTML=sites.map(site=>{let values={};for(const row of site.outcomes){total++;if(row.result.status==='unknown')unknown++;arr(row.result.missing_inputs).forEach(x=>missing.add(typeof x==='string'?x:JSON.stringify(x)));values[row.outcome]=row.result.status;}return `<tr><td>${esc(site.label||site.site_id)}</td><td>${badge(values.deaths||'unknown')}</td><td>${badge(values.serious_nonfatal_injuries||'unknown')}</td></tr>`;}).join('');$('outcomeSummary').textContent=`${unknown} of ${total} site/outcome results are unknown. No missing exposure or risk value is replaced with zero.`;$('outcomeMissing').textContent='Required field inputs: '+[...missing].map(human).join('; ')+'.';}
function bindInteractions(){$('phase').max=String(frames.length-1);$('phase').addEventListener('input',()=>{stopPlay();setFrame(Number($('phase').value));});$('play').addEventListener('click',play);const groups=[];frames.forEach((f,i)=>{const key=f.stage_id||f.id||String(i);if(!groups.some(g=>g.key===key))groups.push({key,i,label:f.stage_label||f.label||human(key)});});$('phaseButtons').innerHTML=groups.map(g=>`<button data-phase-index="${g.i}" aria-pressed="false">${esc(g.label)}</button>`).join('');document.querySelectorAll('[data-phase-index]').forEach(b=>b.addEventListener('click',()=>{stopPlay();setFrame(Number(b.dataset.phaseIndex));}));for(const id of ['parcels','canopy','widths','access','swept'])$(id).addEventListener('change',()=>{drawMap();drawScene();});$('mapFit').addEventListener('click',mapFit);$('mapDistrict').addEventListener('click',mapDistrict);$('mapZoomOut').addEventListener('click',()=>mapZoom(.8));$('mapZoomIn').addEventListener('click',()=>mapZoom(1.25));$('sceneReset').addEventListener('click',()=>{orbit={yaw:-.70,pitch:.52,zoom:1};drawScene();});$('sceneZoomOut').addEventListener('click',()=>{orbit.zoom=Math.max(.55,orbit.zoom*.85);drawScene();});$('sceneZoomIn').addEventListener('click',()=>{orbit.zoom=Math.min(2,orbit.zoom*1.15);drawScene();});scene.addEventListener('pointerdown',e=>{sceneDrag={x:e.clientX,y:e.clientY,yaw:orbit.yaw,pitch:orbit.pitch};scene.setPointerCapture(e.pointerId);});scene.addEventListener('pointermove',e=>{if(sceneDrag){orbit.yaw=sceneDrag.yaw+(e.clientX-sceneDrag.x)*.009;orbit.pitch=Math.max(.18,Math.min(1.2,sceneDrag.pitch+(e.clientY-sceneDrag.y)*.006));drawScene();}});for(const name of ['pointerup','pointercancel'])scene.addEventListener(name,()=>sceneDrag=null);scene.addEventListener('keydown',e=>{let changed=true;if(e.key==='ArrowLeft')orbit.yaw-=.1;else if(e.key==='ArrowRight')orbit.yaw+=.1;else if(e.key==='ArrowUp')orbit.pitch=Math.min(1.2,orbit.pitch+.08);else if(e.key==='ArrowDown')orbit.pitch=Math.max(.18,orbit.pitch-.08);else if(e.key==='+'||e.key==='=')orbit.zoom=Math.min(2,orbit.zoom*1.1);else if(e.key==='-')orbit.zoom=Math.max(.55,orbit.zoom/1.1);else changed=false;if(changed){e.preventDefault();drawScene();}});map.addEventListener('pointerdown',e=>{mapDrag={x:e.clientX,y:e.clientY,cx:mapView.cx,cy:mapView.cy};map.setPointerCapture(e.pointerId);});map.addEventListener('pointermove',e=>{if(mapDrag){mapView.cx=mapDrag.cx-(e.clientX-mapDrag.x)/mapView.scale;mapView.cy=mapDrag.cy+(e.clientY-mapDrag.y)/mapView.scale;drawMap();}});for(const name of ['pointerup','pointercancel'])map.addEventListener(name,()=>mapDrag=null);map.addEventListener('wheel',e=>{e.preventDefault();mapZoom(e.deltaY<0?1.12:.89);},{passive:false});map.addEventListener('keydown',e=>{let changed=true;if(e.key==='ArrowLeft')mapView.cx-=25/mapView.scale;else if(e.key==='ArrowRight')mapView.cx+=25/mapView.scale;else if(e.key==='ArrowUp')mapView.cy+=25/mapView.scale;else if(e.key==='ArrowDown')mapView.cy-=25/mapView.scale;else if(e.key==='+'||e.key==='=')mapZoom(1.2);else if(e.key==='-')mapZoom(1/1.2);else changed=false;if(changed){e.preventDefault();drawMap();}});window.addEventListener('resize',()=>{drawScene();drawMap();drawSide();});}
function initialize(){if(!frames.length||!members().length)throw Error('Canonical scene export is incomplete');$('spanMetric').textContent=fmt(S.span_m,'m');$('massMetric').textContent=fmt(S.assembly_mass_kg,'kg',0);$('phaseMetric').textContent=fmt(S.walking_height_m,'m',3);$('openingMetric').textContent='Unresolved';$('siteTitle').textContent=S.site_label||'Tillamook reference placement';$('mapNote').textContent=`${nf(S.station_offset_m,0)} m station offset · ${S.site_architecture_label}. 3 m is sensitivity, not survey error.`;$('designHash').textContent=D.designSha256;$('sourceHashes').textContent=JSON.stringify(D.sourceHashes,null,2);$('sceneReceipt').textContent=JSON.stringify({scope:S.scope,coordinates:S.coordinates,frameCount:frames.length,memberCount:members().length,mapTransform:S.map_transform,sourceReport:S.source_report},null,2);$('sourceSummary').textContent=`Embedded source geometry: ${A.counts.buildings} buildings, ${A.counts.parcels} parcel parts and ${A.counts.streets} street/path records. Atlas source ${A.provenance.source_archive}; ${A.provenance.manifest_files_verified} source files were verified in v2.`;$('integrationFindings').innerHTML=miniRows(S.integration_findings||[]);renderInterfaces();renderEvidence();renderOutcomes();bindInteractions();mapFit();updateFrame();window.__workbenchReady=true;}
try{initialize();}catch(e){const box=document.createElement('div');box.className='fail';box.textContent='The workbench cannot render this report: '+e.message;document.querySelector('main').prepend(box);window.__workbenchError=e.message;throw e;}
</script></body></html>'''


def main():
    parser = argparse.ArgumentParser(description=__doc__)
    parser.parse_args()
    data = build_payload()
    encoded = json.dumps(data, ensure_ascii=False, separators=(",", ":"))
    encoded = encoded.replace("<", "\\u003c").replace(">", "\\u003e").replace("&", "\\u0026")
    output = ROOT / "workbench.html"
    output.write_text(HTML.replace("__DATA__", encoded), encoding="utf-8")
    print(json.dumps({"path": str(output), "bytes": output.stat().st_size,
                      "design_sha256": data["designSha256"],
                      "frames": len(data["scene"]["frames"])}))


if __name__ == "__main__":
    main()
