"""Rebuild the six raw powder peak windows from a pinned APS synchrotron run."""
import argparse
from decimal import Decimal as D
import hashlib
import json
import math
from pathlib import Path
from urllib.request import urlopen

ROOT=Path(__file__).resolve().parent
URL='https://raw.githubusercontent.com/AdvancedPhotonSource/GSAS-II-tutorials/494c2700cde7c92591e328fdce1515c662cf11ab/FitPeaks/data/11bmb_3844.fxye'

def extract(source):
    lines=source.decode('ascii').splitlines()
    start=next(i for i,x in enumerate(lines) if x.startswith('BANK'))
    wavelength=float(next(s for s in lines[:start] if 'Mean wavelength is:' in s).split(':')[-1])
    xy=[(float(v[0])/100,float(v[1])) for v in (s.split() for s in lines[start+1:] if s.strip())]
    # Strongest separated maxima in a fixed low-angle window, selected without indexing.
    peaks=[]
    for angle,count in sorted(((a,y) for a,y in xy if 5<=a<14),key=lambda v:-v[1]):
        if all(abs(angle-prev)>.4 for prev in peaks):
            peaks.append(angle)
        if len(peaks)==6: break
    peaks.sort()
    def spacing(angle): return wavelength/(2*math.sin(math.radians(angle/2)))
    # ±0.01° allowance chosen for this replay; rounded out to 1e-6 angstrom.
    bounds=[[format(math.floor(spacing(a+.01)*1e6)*D('0.000001'),'f'),
             format(math.ceil(spacing(a-.01)*1e6)*D('0.000001'),'f')] for a in peaks]
    return {'mode':'cubic_powder','cell_range_A':[3,9],'d_spacing_intervals_A':bounds,
            'source':{'material':'NIST SRM 660a LaB6, APS 11-BM run 3844','url':URL,
            'sha256':hashlib.sha256(source).hexdigest(),'wavelength_A':wavelength,
            'selected_two_theta_deg':[round(a,4) for a in peaks],'assumed_half_width_deg':0.01}}

if __name__=='__main__':
    p=argparse.ArgumentParser(description=__doc__);p.add_argument('file',nargs='?',type=Path);args=p.parse_args()
    if args.file: src=args.file.read_bytes()
    else:
        with urlopen(URL,timeout=35) as r: src=r.read()
    actual=extract(src)
    if actual!=json.loads((ROOT/'powder.json').read_text()): raise SystemExit('Source or extraction differs')
    print('Verified six measured peaks and source hash.')
