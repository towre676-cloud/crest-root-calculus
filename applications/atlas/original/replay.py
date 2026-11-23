"""Replay Part XXX decisions against the recorded catalogue and two controls."""
import json
from hashlib import sha256
from pathlib import Path
import sys

HERE=Path(__file__).resolve().parent
sys.path.insert(0,str(HERE.parent))
from app import solve
from algebra import parse_form

def verify(case, records):
    if len(case['curves'])!=len(records):
        raise AssertionError("Catalogue and recorded case counts differ.")
    curves=solve(case)['curves']
    for source,row,expected in zip(case['curves'],curves,records):
        assert source['id']==row['id']==expected['id']
        original={tuple(m[:3]):m[3] for m in expected['terms']}
        assert parse_form(source['quartic'])==original
        if 'proof' in expected:
            if json.loads(json.dumps(row['proof']))!=expected['proof']:
                raise AssertionError('The original proof did not reproduce for '+row['id'])
        else:
            assert row['status']=='unresolved' and row.get('reason') and 'proof' not in row
    return {r['id']:r['status'] for r in curves}

if __name__=='__main__':
    if not __debug__:raise SystemExit('Run without -O to retain assertions.')
    def read(path):return json.loads(path.read_text(encoding='utf-8'))
    for k,digest in read(HERE/'source-hashes.json')['files'].items():
        assert sha256((HERE/'vendor/crest'/k).read_bytes()).hexdigest()==digest
    report=verify(read(HERE.parent/'original-slice.json'),read(HERE/'catalogue.json'))
    for n in ('C-star','nonlift'):
        report.update(verify(read(HERE.parent/('original-'+n+'.json')),[read(HERE/(n+'.json'))]))
    print(json.dumps(report,indent=2))
