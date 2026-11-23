"""Screen second-life capacity grades using measured voltage traces and bounded mismatch."""
import csv
import json
import sys
from decimal import Decimal, InvalidOperation
from fractions import Fraction
from pathlib import Path

ROOT = Path(__file__).resolve().parent
STAGES = (60, 300, 600, 1200, 1800)

def number(value, name):
    if isinstance(value, bool):
        raise ValueError(name + ' must be a finite number')
    try:
        d = Decimal(str(value))
    except (ValueError, InvalidOperation):
        raise ValueError(name + ' must be a finite number') from None
    if not d.is_finite():
        raise ValueError(name + ' must be finite')
    return d

def library():
    with (ROOT / 'data/measurements.csv').open(encoding='utf8') as stream:
        return [{'id': row['file'], 'capacity_ah': row['capacity_ah'],
                 'volts': {str(t): row['volts_'+str(t)+'s'] for t in STAGES}}
                for row in csv.DictReader(stream)]

def prepare(rows):
    out = []
    ids = set()
    for row in rows:
        ident = row['id']
        if not isinstance(ident, str) or not ident or ident in ids:
            raise ValueError('Reference identifiers must be distinct nonempty strings')
        ids.add(ident)
        cap = number(row['capacity_ah'], 'capacity_ah')
        trace = {t: number(row['volts'][str(t)], 'reference volts') for t in STAGES}
        if cap <= 0 or min(trace.values()) <= 0:
            raise ValueError('Reference capacities and voltages must be positive')
        out.append((ident, cap, trace))
    if not out:
        raise ValueError('A nonempty reference library is required')
    return out

def evaluate(data, refs=None):
    mode = data.get('policy','adaptive')
    if mode not in ('adaptive','fixed_300s'):
        raise ValueError('policy must be adaptive or fixed_300s')
    limit = number(data['min_capacity_ah'], 'min_capacity_ah')
    tol = number(data['voltage_tolerance_v'], 'voltage_tolerance_v')
    caperr = number(data['capacity_error_ah'], 'capacity_error_ah')
    if limit <= 0 or min(tol, caperr) < 0:
        raise ValueError('Positive capacity threshold and nonnegative error limits required')
    obs = {}
    for stage, value in data['observations'].items():
        if stage not in {str(t) for t in STAGES}:
            raise ValueError('Unknown checkpoint '+str(stage))
        obs[int(stage)] = number(value, 'observed voltage')
    if mode=='fixed_300s' and any(t>300 for t in obs):
        raise ValueError('fixed_300s uses only observations at 60 and 300 seconds')
    if obs and min(obs.values()) <= 0:
        raise ValueError('Observed voltages must be positive')
    if refs is None:
        refs = prepare(data['reference'] if 'reference' in data else library())
    survivors = [r for r in refs if all(abs(r[2][t]-v) <= tol for t, v in obs.items())]
    result = {'candidate_count':len(survivors), 'candidate_ids':[r[0] for r in survivors],
              'capacity_range_ah':None, 'next_checkpoint_s':None, 'checkpoint_options':[]}
    if not survivors:
        return dict(result, decision='reference_test', reason='No measured reference matches all checkpoints.')
    low = max(Decimal(0), min(r[1] for r in survivors)-caperr)
    high = max(r[1] for r in survivors)+caperr
    result['capacity_range_ah'] = [float(low), float(high)]
    if mode == 'fixed_300s':
        required = [t for t in (60,300) if t not in obs]
        if required:
            if max(obs,default=0)>required[0]:
                return dict(result,decision='reference_test',reason='A required five-minute-screen checkpoint was missed.')
            return dict(result,decision='continue_screen',next_checkpoint_s=required[0],
                        reason='The fixed five-minute policy uses readings from both 60 and 300 seconds.')
    if low >= limit or high < limit:
        return dict(result, decision='meets_capacity_grade' if low >= limit else 'below_capacity_grade',
                    reason='All surviving reference intervals agree on this grade.')
    if mode == 'fixed_300s':
        return dict(result,decision='reference_test',reason='The five-minute screening budget is spent and the grade remains unresolved.')
    now = max(obs, default=0)
    pairs = [(a,b) for i,a in enumerate(survivors) for b in survivors[i+1:]
             if min(a[1], b[1])-caperr < limit <= max(a[1], b[1])+caperr]
    ranked = []
    for t in STAGES:
        if t > now:
            separated = sum(abs(a[2][t]-b[2][t]) > 2*tol for a,b in pairs)
            ranked.append((Fraction(separated,t-now), separated, -t, t))
            result['checkpoint_options'].append({'elapsed_s':t, 'separated_cross_grade_pairs':separated,
                                                  'additional_s':t-now})
    if not ranked or (not pairs and len(survivors)==1):
        return dict(result, decision='reference_test', reason='The grade remains unresolved. Use a reference capacity test.')
    best = max(ranked)
    result['next_checkpoint_s'] = best[-1]
    if best[1] == 0:
        return dict(result, decision='continue_screen',
                    reason='No available checkpoint guarantees a cross-grade split. A particular reading may still discard rivals. Try the earliest remaining checkpoint, or use a full reference test when bench time is scarce.')
    return dict(result, decision='continue_screen',
                reason='Rank future checkpoints by separated cross-grade pairs per additional second; ties use the larger separation count, then earlier time.')

def main():
    if len(sys.argv)>2:
        sys.exit('Usage: python app.py [input.json]')
    try:
        path = Path(sys.argv[1]) if len(sys.argv)>1 else ROOT/'example.json'
        print(json.dumps(evaluate(json.loads(path.read_text(encoding='utf8'), parse_float=str)), indent=2))
    except (ValueError, KeyError, TypeError, AttributeError, OSError) as exc:
        sys.exit('Input error: '+str(exc))
if __name__ == '__main__':
    main()
