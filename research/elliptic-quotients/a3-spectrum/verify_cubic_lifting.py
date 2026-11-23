"""Cross-check the four geometric cubic classes with the original Part XXX solver."""
import json
from pathlib import Path
import sys

from map_classes import count_classes


HERE = Path(__file__).resolve().parent
ENGINE = HERE.parents[2] / "applications/atlas/original/vendor/crest/scan_completeness_xxx/decision"
sys.path.insert(0, str(ENGINE))
from decide import decide


def main():
    source = json.loads((HERE.parent/"primitive-degree-four/ciani.json").read_text(encoding="utf-8"))
    result = decide(source, 11, "geometric", max_precision=3)
    # Normalize integer-key dictionaries and tuples exactly as the archived replay does.
    actual = json.loads(json.dumps(result))
    saved = json.loads((HERE / "cubic-lifting.json").read_text(encoding="utf-8"))
    if actual != saved:
        raise ArithmeticError("The archived Part XXX run did not reproduce")
    upper = [row["upper"] for row in result["levels"]]
    primitive = count_classes(3)["geometrically_primitive_classes"]
    if upper != [8, 8, 4] or upper[-1] != primitive:
        raise ArithmeticError("Local upper and global lower bounds disagree")
    print("Verified local upper bounds 8, 8, 4 and four distinct geometric quotient directions.")


if __name__ == "__main__":
    main()
