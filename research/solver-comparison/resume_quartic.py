"""Resume the published interrupted quartic study without repeating measured rows.

The original protocol and adapter remain immutable. A separate resume record
stores the exact run boundary and hashes of the imported research code.
"""
from datetime import datetime, timezone
import json
from pathlib import Path
import random
import argparse

from benchmark import METHODS, run_one, sha256, summarize, platform_record

HERE = Path(__file__).resolve().parent


def dependencies():
    names = ["chart.py", "proof.py", "solve.py"]
    return {"chart-one-cofactors/"+name: sha256(HERE.parent/"chart-one-cofactors"/name)
            for name in names}


def main(output, cpu):
    panel = HERE/"quartic-panel.json"
    protocol_path = output.with_name(output.stem+"-protocol.json")
    protocol = json.loads(protocol_path.read_text())
    data = json.loads(panel.read_text())
    payload = json.loads(output.read_text())
    if payload["completed"]:
        raise ValueError("This study is already complete.")
    if payload["protocol_sha256"] != sha256(protocol_path):
        raise ArithmeticError("The original protocol changed.")
    if protocol["worker_sha256"] != sha256(HERE/"worker.py") or protocol["panel_sha256"] != sha256(panel):
        raise ArithmeticError("The original input or adapter changed.")
    if platform_record() != protocol["platform"]:
        raise ArithmeticError("The study must resume on the same platform and software.")
    rng = random.Random(76123)
    schedule = []
    for rep in (1, 2):
        order = [(case["id"], method) for case in data["cases"] for method in METHODS]
        rng.shuffle(order)
        schedule += [(case, method, rep) for case, method in order]
    rows, supports = payload["runs"], payload["supports"]
    if [(r["case"],r["method"],r["replicate"]) for r in rows] != schedule[:len(rows)]:
        raise ArithmeticError("Recorded rows are not the exact frozen schedule prefix.")
    checkpoint = {"resumed_utc": datetime.now(timezone.utc).isoformat(),
                  "completed_rows_preserved": len(rows), "worker_cpu": cpu,
                  "research_source_sha256": dependencies(),
                  "reason": "Execution session interrupted after the last durable row; resume continues the exact seeded schedule and preserves every measured row.",
                  "hash_scope": "Research dependency hashes recorded on resumption and checked at completion. Original protocol already froze the input and adapter."}
    payload.setdefault("resumptions", []).append(checkpoint)
    output.write_text(json.dumps(payload, indent=2)+"\n")
    for case, method, rep in schedule[len(rows):]:
        row, result = run_one(case, method, rep, panel,
                              protocol["per_process_wall_limit_seconds"], cpu)
        rows.append(row)
        if result is not None:
            reduced = {key: result[key] for key in ("basis", "geometric_marks")}
            if case in supports and supports[case] != reduced:
                raise ArithmeticError("Exact output disagreement on "+case)
            supports[case] = reduced
        payload.update({"runs": rows, "supports": supports,
                        "summary": summarize(rows), "completed": len(rows)==len(schedule)})
        output.write_text(json.dumps(payload, indent=2)+"\n")
        print(case, method, rep, row["status"], round(row["process_wall_seconds"],3), flush=True)
    checkpoint["completed_utc"] = datetime.now(timezone.utc).isoformat()
    checkpoint["dependencies_unchanged"] = dependencies()==checkpoint["research_source_sha256"]
    output.write_text(json.dumps(payload, indent=2)+"\n")
    if not checkpoint["dependencies_unchanged"]:
        raise ArithmeticError("A research dependency changed during the resumed run.")
    print(json.dumps(summarize(rows), indent=2), flush=True)


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--output", type=Path, default=HERE/"quartic-results.json")
    parser.add_argument("--cpu", type=int, default=8)
    args = parser.parse_args()
    main(args.output, args.cpu)
