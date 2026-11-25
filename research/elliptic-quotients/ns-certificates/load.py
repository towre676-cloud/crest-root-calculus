"""Load certified maps (Proposition 6.1 records) from the period-reconstruction replay."""
import json
from pathlib import Path
import sys

from ffcore import gauss, poly_gauss

PERIOD = Path(__file__).resolve().parent.parent / "period-reconstruction"
sys.path.insert(0, str(PERIOD))
from certificate import certify_record

REPLAY = PERIOD / "replay.json"


def load_curve(prefix):
    d = json.load(open(REPLAY))["curves"]
    name = next(n for n in d if n.startswith(prefix))
    rec = d[name]
    maps = []
    for m in rec["maps"]:
        degree, _, _, _ = certify_record(rec["f"], m)
        diff = m["differential"]
        # "(A) dx/y + (B) x dx/y"
        A = diff.split(") dx/y")[0].lstrip("(")
        B = diff.split("+ (")[1].split(") x dx/y")[0]
        maps.append({"N": poly_gauss(m["N"]), "D": poly_gauss(m["D"]), "a": gauss(A), "b": gauss(B),
                     "g2": gauss(m["g2"]), "g3": gauss(m["g3"]), "degree": degree,
                     "gaussian": "w" in m["N"] + m["D"] + diff + m["g2"] + m["g3"]})
    return name, rec["f"], maps
