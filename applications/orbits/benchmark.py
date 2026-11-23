"""Replay visibility-screened planning versus earliest slot; quantify the actual RV sampling alias."""
import cmath
import hashlib
import json
import math
from app import ROOT, evaluate

SOURCES={'UID_0043587_RVC_001.tbl':'85622e16449b95ac0033dac6fc7c5858dc38e3643bc5d5ed73ac5c9f973c5639',
         'UID_0043587_RVC_002.tbl':'28483eb05e8728a0fc0315cbe60496a3b5fa461f8819d1b6620c47f1a9436839'}

def bench():
    obs=[]
    for file, required in SOURCES.items():
        path=ROOT/'data'/file
        if hashlib.sha256(path.read_bytes()).hexdigest()!=required:
            raise ValueError('RV file checksum mismatch: '+file)
        for line in path.read_text().splitlines():
            if line.strip() and line.strip()[0].isdigit():
                t,rv,err=map(float,line.split())
                obs.append(t+2440000)  # NASA files use the reduced dates of Fischer et al.
    f=1/.736539-1/2.81705
    origin=obs[0]
    z=sum(cmath.exp(2j*math.pi*f*(t-origin)) for t in obs)/len(obs)
    case=json.loads((ROOT/'example.json').read_text())
    forecast=evaluate(case)
    first=min((w for w in case['windows'] if w['available']),key=lambda w:float(w['days'][0]))
    scores={w['name']:len(w['distinguished_pairs']) for w in forecast['windows']}
    observed=evaluate(json.loads((ROOT/'confirmed.json').read_text()))
    return {'measured_RV_timestamps':len(obs),'frequency_difference_per_day':round(f,8),
            'spectral_window_magnitude_at_alias_frequency':round(abs(z),8),
            'historical_schedule_example':{'first_available':first['name'],'first_pairs':scores[first['name']],
            'selected':forecast['recommended_window'],'selected_pairs':scores[forecast['recommended_window']]},
            'confirmation_survivors':observed['surviving_hypotheses']}
if __name__=='__main__': print(json.dumps(bench(),indent=2))
