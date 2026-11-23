"""Verify that the frozen blind slice reproduces from Sutherland's public text file."""
from argparse import ArgumentParser
from hashlib import sha256
import json
from pathlib import Path
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parent

def main():
    ap=ArgumentParser(description=__doc__)
    ap.add_argument("source", type=Path, nargs="?", help="use an already downloaded full text file")
    args=ap.parse_args()
    ref=json.loads((ROOT/"example.json").read_text())
    spec=ref["selection"]
    if args.source:
        source=args.source.read_bytes()
    else:
        with urlopen(spec["url"], timeout=40) as response: source=response.read()
    if sha256(source).hexdigest()!=spec["sha256"]: raise SystemExit("Source hash differs")
    lines=source.decode("ascii").splitlines()
    pinned=[x for x in ref["curves"] if "row" in x]
    rows=sorted(range(1,len(lines)+1),key=lambda i:sha256(f'{spec["seed"]}|{i}'.encode()).digest())[:len(pinned)]
    if len(lines)!=spec['row_count'] or rows!=[x['row'] for x in pinned]: raise SystemExit("Sample differs")
    for i,x in zip(rows,pinned):
        if lines[i-1]!=f'{x["absolute_discriminant"]}:[{x["quartic"]}]': raise SystemExit(f"Record differs at row {i}")
    print(f"Verified {len(rows)} selected records out of {len(lines)}.")
if __name__=="__main__": main()
