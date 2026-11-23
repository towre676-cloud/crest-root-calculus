"""Extract measured voltage checkpoints from Ji et al., Zenodo 18630889 (CC BY 4.0)."""
import argparse
import csv
import datetime as dt
import hashlib
import io
import re
from pathlib import Path
from urllib.request import urlopen
from xml.etree import ElementTree as ET
from zipfile import ZipFile

URL = 'https://zenodo.org/api/records/18630889/files/SLBs_LFP_charging_data.zip/content'
SHA256 = '16afff956da104fcfb0acae5ee86a53f4c5c5e030c06f6fa8e7dfd566c715d07'
STAGES = (60, 300, 600, 1200, 1800)
NS = '{http://schemas.openxmlformats.org/spreadsheetml/2006/main}'

def extract(raw):
    if hashlib.sha256(raw).hexdigest() != SHA256:
        raise ValueError('Source ZIP differs from the pinned version')
    results = []
    with ZipFile(io.BytesIO(raw)) as bundle:
        for name in sorted(bundle.namelist()):
            match = re.fullmatch(r'(\d+(?:\.\d+)?)Ah_charge.xlsx', name)
            if not match:
                raise ValueError('Unexpected source filename ' + name)
            with ZipFile(io.BytesIO(bundle.read(name))) as workbook:
                strings = [item.find(NS+'t').text for item in ET.fromstring(workbook.read('xl/sharedStrings.xml'))]
                sheet = ET.fromstring(workbook.read('xl/worksheets/sheet1.xml'))
                samples = {}
                origin = None
                for row in sheet.find(NS+'sheetData')[1:]:
                    fields = {re.match('[A-Z]+', c.get('r')).group(): c for c in row}
                    when = dt.datetime.fromisoformat(strings[int(fields['A'].find(NS+'v').text)])
                    origin = origin or when
                    elapsed = int((when-origin).total_seconds())
                    if elapsed in STAGES:
                        samples[elapsed] = round(float(fields['B'].find(NS+'v').text), 6)
                if set(samples) != set(STAGES):
                    raise ValueError('Missing staged reading in '+name)
            results.append([name, match.group(1), *[samples[t] for t in STAGES]])
    return results

if __name__ == '__main__':
    ap = argparse.ArgumentParser(description=__doc__)
    ap.add_argument('archive', nargs='?', help='Existing source ZIP; otherwise download DOI data')
    ap.add_argument('--out', type=Path, required=True)
    args = ap.parse_args()
    raw = Path(args.archive).read_bytes() if args.archive else urlopen(URL, timeout=90).read()
    with args.out.open('w', newline='', encoding='utf8') as f:
        writer = csv.writer(f, lineterminator='\n')
        writer.writerow(['file', 'capacity_ah', *['volts_'+str(t)+'s' for t in STAGES]])
        writer.writerows(extract(raw))
