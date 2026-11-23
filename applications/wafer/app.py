"""Decide a bottom-width limit or select another NIST optical reading."""
import argparse
import csv
from bisect import bisect_left, bisect_right
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from functools import lru_cache
import gzip
from hashlib import sha256
from io import StringIO
import json
from pathlib import Path

ROOT = Path(__file__).resolve().parent
SCALE = 10 ** 13
HASHES = {
    "Pi_L100P300_sim_All": "d0482e4e0ab0a053fbf9bd7bfe58df1339ac9d4bfb54ce5c31b958534b8ac818",
    "Phi_L100P300_sim_All": "d8d08a758c3691de0a9084b710f2fa57e935649b22907eb2db1636f7fe4c0f2d",
    "Imean_L100P300_exp": "d60502b0cc12e76c37e2df491e2e21363d2a9eb55087e6af8c6795a15af1829b",
}
ANGLES_X = [-36.65, -32.23, -28.07, -24.12, -20.34, -16.7, -13.19, -9.78, -6.45, -3.19, 0, 3.2, 6.45, 9.78, 13.2, 16.71, 20.35, 24.12, 28.07, 32.22, 36.65]
ANGLES_Y = [-36.65, -32.23, -28.07, -24.11, -20.34, -16.7, -13.2, -9.78, -6.45, -3.19, 0, 3.2, 6.46, 9.78, 13.2, 16.71, 20.35, 24.12, 28.07, 32.22, 36.65]

def number(value, name, minimum, maximum):
    try:
        n = Decimal(str(value))
    except (InvalidOperation, ValueError):
        raise ValueError(f"{name} must be a decimal number") from None
    if not n.is_finite() or not minimum <= n <= maximum or n.as_tuple().exponent < -13:
        raise ValueError(f"{name} must be between {minimum} and {maximum}, with at most 13 decimal places")
    return n

def scaled(value):
    return int(number(value, "Intensity or error", 0, 1) * SCALE)

def channel(value):
    if isinstance(value, bool) or str(value) not in {str(i) for i in range(1,85)}:
        raise ValueError("Channel must be a whole number from 1 to 84")
    return int(value)

def channels(value):
    if value == "all":
        return list(range(1, 85))
    if not isinstance(value, list):
        raise ValueError("Specify channels as a list of numbers or the word all")
    out = [channel(v) for v in value]
    if len(set(out)) != len(out):
        raise ValueError("Duplicate channels")
    return out

def label(k):
    group, index = divmod(k - 1, 21)
    pol = 'X' if group in (0,2) else 'Y'
    plane = 'X' if group < 2 else 'Y'
    angle = (ANGLES_X if plane == 'X' else ANGLES_Y)[index]
    return f"{pol} polarization, {angle:g}° in the {plane} plane"

def read_table(name):
    raw = gzip.decompress((ROOT/'data'/(name+'.csv.gz')).read_bytes())
    if sha256(raw).hexdigest() != HASHES[name]:
        raise ValueError("NIST source checksum mismatch: "+name)
    return list(csv.reader(StringIO(raw.decode('ascii'))))

@lru_cache(maxsize=1)
def library():
    sizes = [tuple(map(Decimal,r)) for r in read_table("Pi_L100P300_sim_All")]
    signals = [tuple(int(Decimal(v)*SCALE) for v in r) for r in read_table("Phi_L100P300_sim_All")]
    dies = [tuple(int(Decimal(v)*SCALE) for v in r) for r in read_table("Imean_L100P300_exp")]
    return sizes, signals, dies

def suggestions(good, bad, signals, options, costs, error):
    choices = []
    for k in options:
        high = sorted(signals[i][k-1] for i in bad)
        overlap = 0
        for i in good:
            centre = signals[i][k-1]
            overlap += bisect_right(high, centre+2*error)-bisect_left(high, centre-2*error)
        separated = len(good)*len(bad)-overlap
        cost = costs.get(k,Decimal(1))
        choices.append((Fraction(separated)/Fraction(cost), separated, -k,
                        {"channel":k, "setting":label(k), "cost":float(cost),
                         "opposite_decision_pairs_separated":separated,
                         "opposite_decision_pairs_left":overlap}))
    return [item[3] for item in sorted(choices, reverse=True)]

def solve(data):
    if not isinstance(data, dict):
        raise ValueError("Input must be an object")
    sizes, signals, dies = library()
    limit = number(data["width_limit_nm"], "Width limit", 0, 10000)
    error = scaled(data["error"])
    options = channels(data.get("available_channels", "all"))
    costs = {channel(k):number(v,"Cost",Decimal("0.000001"),1000000) for k,v in data.get("channel_costs",{}).items()}
    if "die" in data:
        if "observations" in data: raise ValueError("Choose die readings or supplied observations")
        die = data["die"]
        if type(die) is not int or not 1 <= die <= 9: raise ValueError("Die must be an integer from 1 to 9")
        used = channels(data.get("observed_channels","all"))
        obs = {k:dies[die-1][k-1] for k in used}
        source = f"NIST measured die {die}"
    else:
        obs = {channel(k):scaled(v) for k,v in data["observations"].items()}
        source = "Supplied readings"
    residuals = [max((abs(r[k-1]-v) for k,v in obs.items()),default=0) for r in signals]
    keep = [i for i,d in enumerate(residuals) if d<=error]
    good = [i for i in keep if sizes[i][0]<=limit]
    bad = [i for i in keep if sizes[i][0]>limit]
    if not keep:
        status = "outside_model"
        decision = "No library profile explains these readings within the specified error."
    elif bad and good:
        status = "ambiguous"
        decision = "Compatible profiles fall on both sides of the width limit."
    elif good:
        status = "within_limit"
        decision = "Every compatible library profile meets the bottom-width limit."
    else:
        status = "above_limit"
        decision = "Every compatible library profile exceeds the bottom-width limit."
    ranked = suggestions(good,bad,signals,sorted(set(options)-obs.keys()),costs,error) if good and bad else []
    bounds = [float(f(sizes[i][0] for i in keep)) for f in (min,max)] if keep else None
    examples = []
    for label_, subset in [('within',good),('above',bad)]:
        if subset:
            i = min(subset,key=lambda i:(residuals[i],i))
            examples.append({"class":label_,"row":i+1,"bottom_middle_top_nm":list(map(float,sizes[i])),
                             "largest_residual":residuals[i]/SCALE})
    return {"status":status,"decision":decision,"source":source,"width_limit_nm":float(limit),
            "error":error/SCALE,"observed_channel_count":len(obs),"surviving_profiles":len(keep),
            "within":len(good),"above":len(bad),"bottom_width_range_nm":bounds,
            "smallest_maximum_intensity_discrepancy":min(residuals)/SCALE,"witnesses":examples,
            "opposite_decision_pairs":len(good)*len(bad),
            "next_reading":ranked[0] if ranked and ranked[0]["opposite_decision_pairs_separated"] else None,
            "ranked_readings":ranked,
            "scope":"Discrete NIST library and supplied absolute error on each channel."}

if __name__ == "__main__":
    cli=argparse.ArgumentParser(description=__doc__)
    cli.add_argument("input",nargs="?",type=Path,default=ROOT/"example.json")
    args=cli.parse_args()
    try:
        print(json.dumps(solve(json.loads(args.input.read_text(encoding="utf-8"))),indent=2,allow_nan=False))
    except (ValueError,KeyError,TypeError,AttributeError,OSError) as e:
        cli.exit(2,f"Input error: {e}\n")
