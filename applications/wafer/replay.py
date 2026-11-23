"""Reproduce a finite-library optical inspection exercise. Python 3.9+, standard library."""
from bisect import bisect_left, bisect_right
import csv
import gzip
from io import StringIO
from decimal import Decimal
from hashlib import sha256
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
FILES = {
    "Imean_L100P300_exp.csv": "d60502b0cc12e76c37e2df491e2e21363d2a9eb55087e6af8c6795a15af1829b",
    "Pi_L100P300_sim_All.csv": "d0482e4e0ab0a053fbf9bd7bfe58df1339ac9d4bfb54ce5c31b958534b8ac818",
    "Phi_L100P300_sim_All.csv": "d8d08a758c3691de0a9084b710f2fa57e935649b22907eb2db1636f7fe4c0f2d",
}
SCALE = 10**13
MAX_PAIR_GAP = 2 * SCALE // 1000

def rows(name):
    path = ROOT / "data" / (name + ".gz")
    raw = gzip.decompress(path.read_bytes())
    if sha256(raw).hexdigest() != FILES[name]:
        raise ValueError("Unexpected input bytes in " + name)
    return list(csv.reader(StringIO(raw.decode("ascii"))))

def fixed(value):
    scaled = Decimal(value) * SCALE
    if scaled != int(scaled):
        raise ValueError("Nonrepresentable intensity " + value)
    return int(scaled)

def main():
    geometry = rows("Pi_L100P300_sim_All.csv")
    profiles = [[fixed(x) for x in row] for row in rows("Phi_L100P300_sim_All.csv")]
    assert len(geometry) == len(profiles) == 2566 and all(len(x) == 84 for x in profiles)
    in_spec = [Decimal(row[0]) <= 135 for row in geometry]
    # Zero-based channel 10 is normal incidence with X polarization.
    opposite = sorted((p[10], i) for i, p in enumerate(profiles) if not in_spec[i])
    values = [x for x, _ in opposite]
    pairs = []
    for i, p in enumerate(profiles):
        if in_spec[i]:
            lo = bisect_left(values, p[10] - MAX_PAIR_GAP)
            hi = bisect_right(values, p[10] + MAX_PAIR_GAP)
            pairs.extend((i, j) for _, j in opposite[lo:hi])
    recipes = {"normal_X": [], "normal_X_and_nearby_X": [11],
               "normal_X_and_normal_Y": [31], "all_84": list(range(84))}
    report = {}
    for name, extra in recipes.items():
        overlaps = [(i, j) for i, j in pairs
                    if all(abs(profiles[i][k] - profiles[j][k]) <= MAX_PAIR_GAP for k in extra)]
        ambiguous = {p for pair in overlaps for p in pair}
        report[name] = {"decidable_profiles": len(profiles)-len(ambiguous),
                        "opposite_class_pairs_with_overlapping_readings": len(overlaps)}
    residuals = []
    measured = rows("Imean_L100P300_exp.csv")
    assert len(measured) == 9 and all(len(row) == 84 for row in measured)
    for row in measured:
        observation = [fixed(x) for x in row]
        best = min(max(abs(x-y) for x, y in zip(observation, profile)) for profile in profiles)
        residuals.append(str(Decimal(best)/SCALE))
    report["nine_measured_dies"] = {
        "compatible_at_0.001_direct_table_comparison": sum(Decimal(r) <= Decimal(".001") for r in residuals),
        "smallest_max_channel_discrepancy_per_die": residuals,
    }
    print(json.dumps(report, indent=2))
    return report

if __name__ == "__main__":
    actual = main()
    expected = ROOT / "results.json"
    if expected.exists() and actual != json.loads(expected.read_text()):
        raise SystemExit("Results differ from the recorded calculation.")
