"""Compare reported pre/post coheating outcomes with independent interval subtraction."""

import json
from pathlib import Path
from app import solve

ROOT = Path(__file__).parent


def replay():
    source = json.loads((ROOT / "field.json").read_text())
    example = json.loads((ROOT / "example.json").read_text())
    out = []
    for h in source["homes"]:
        b, be = h["before"]
        a, ae = h["after"]
        different = b - a
        example["total_W_K"] = [b - be, b + be]
        case = solve(example)
        out.append(
            {
                "home": h["id"],
                "measured_change_W_K": different,
                "independent_endpoint_change_W_K": [
                    different - be - ae,
                    different + be + ae,
                ],
                "published_change_and_uncertainty_W_K": h["published_change"],
                "pretest_status": case["status"],
                "pretest_next": case["next_test"],
            }
        )
    return out


if __name__ == "__main__":
    print(json.dumps(replay(), indent=2))
