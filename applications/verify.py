"""Verify the offline release, replay reference studies, and run the unit suites."""
import argparse
from hashlib import sha256
import json
from pathlib import Path
import subprocess
import sys

ROOT=Path(__file__).resolve().parent
CASES={'wafer':('replay.py','results.json'),'water':('field_replay.py','receipt.json'),
       'battery':('benchmark.py','benchmark.json'),'retrofit':('field_replay.py','receipt.json'),
       'materials':('app.py','receipt.json'),'orbits':('benchmark.py','benchmark.json'),
       'atlas':('app.py','receipt.json')}

def compare(folder,cmd,expected):
    p=subprocess.run([sys.executable,'-B',cmd],cwd=folder,stdout=subprocess.PIPE,text=True,encoding='utf-8',check=True)
    if json.loads(p.stdout)!=json.loads((folder/expected).read_text(encoding='utf-8')):
        raise SystemExit(f'{folder.name}/{cmd} differs from {expected}')

def main():
    parser=argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--exact',action='store_true',help='replay the original Part XXX solver with NumPy and SymPy')
    args=parser.parse_args()
    if not __debug__:raise SystemExit('Run without -O so numerical assertions execute.')
    for line in (ROOT/'SHA256SUMS.txt').read_text(encoding='utf-8').splitlines():
        expected,name=line.split('  ',1)
        if sha256((ROOT/name).read_bytes()).hexdigest()!=expected:
            raise SystemExit('Checksum mismatch: '+name)
    print('All release files match their SHA-256 checksums.',flush=True)
    for app,(cmd,receipt) in CASES.items():
        folder=ROOT/app
        subprocess.run([sys.executable,'-B','-m','unittest','discover','-v'],cwd=folder,check=True)
        compare(folder,cmd,receipt)
        print(f'{app}: unit tests and {receipt} passed.',flush=True)
    if args.exact:
        for script in ['replay.py','audit_powers.py']:
            subprocess.run([sys.executable,'-B',script],cwd=ROOT/'atlas/original',check=True)
    print('All requested checks passed.')

if __name__=='__main__':main()
