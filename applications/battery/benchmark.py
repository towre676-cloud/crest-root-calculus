"""Leave each physically distinct cell out and score short-screen decisions against published labels."""
from collections import Counter
import json
from app import STAGES, evaluate, library, prepare

PARAMETERS = {'min_capacity_ah':4.8, 'voltage_tolerance_v':0.020, 'capacity_error_ah':0.050}

def benchmark():
    rows = library()
    output = {}
    for strategy in ('adaptive', 'fixed_300s', 'fixed_1800s'):
        outcomes = Counter()
        samples = []
        for row in rows:
            refs = prepare([r for r in rows if r['id'] != row['id']])
            obs = {}
            if strategy == 'adaptive':
                stage = STAGES[0]
                while True:
                    obs[str(stage)] = row['volts'][str(stage)]
                    result = evaluate(dict(PARAMETERS, observations=obs), refs)
                    if result['decision'] != 'continue_screen':
                        break
                    stage = result['next_checkpoint_s']
            else:
                stage = 300 if strategy=='fixed_300s' else 1800
                obs = {str(t):row['volts'][str(t)] for t in STAGES if t<=stage}
                result = evaluate(dict(PARAMETERS, observations=obs), refs)
            action = result['decision']
            truth = float(row['capacity_ah'])>=PARAMETERS['min_capacity_ah']
            if action not in {'below_capacity_grade','meets_capacity_grade'}:
                outcome='deferred'
            elif (action=='meets_capacity_grade') == truth:
                outcome='correct'
            else:
                outcome='false_accept' if not truth else 'false_reject'
            outcomes[outcome]+=1
            samples.append({'held_out_cell':row['id'],'stop_s':stage,'outcome':outcome})
        output[strategy] = {'total_cells':len(rows), 'short_screen_seconds':sum(r['stop_s'] for r in samples),
                            **{k:outcomes[k] for k in ('correct','false_accept','false_reject','deferred')},
                            'held_out':samples}
    return {'parameters':PARAMETERS,'strategies':output}
if __name__=='__main__':
    print(json.dumps(benchmark(),indent=2))
