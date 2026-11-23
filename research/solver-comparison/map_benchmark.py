"""Compare constant-vertical maps against a frozen, verified exact reference."""
import argparse
from datetime import datetime, timezone
import json
import os
from pathlib import Path
import random
import shutil
import sys
import time

import benchmark
from benchmark import platform_record, run_one, sha256

HERE = Path(__file__).resolve().parent
METHODS = ('crest','singular','sympy-decomposition','corrected-dense-decomposition')


def wait_for_study(path):
    while True:
        try:
            if json.loads(path.read_text()).get('completed'):
                return
        except (FileNotFoundError,json.JSONDecodeError):
            pass
        time.sleep(5)


def prepare_reference(path,panel):
    """Certify every exact answer before starting any timed method."""
    sys.path.insert(0,str(HERE.parent/'bounded-correspondences'))
    from composition import recover_constant,verify_constant
    data = json.loads(panel.read_text())
    if path.exists():
        reference = json.loads(path.read_text())
        if reference['panel_sha256']!=sha256(panel):
            raise ArithmeticError('The existing reference belongs to another input panel.')
    else:
        records = {}
        for case in data['cases']:
            result = recover_constant(case['source'],case['target'])
            records[case['id']] = result
        reference = {'schema':1,'panel_sha256':sha256(panel),
                     'proof':'Complete normalized polynomial-decomposition theorem and exact finite-algebra certificates. Each certificate is verified before timing.',
                     'cases':records,'created_utc':datetime.now(timezone.utc).isoformat()}
        path.write_text(json.dumps(reference,indent=2)+'\n')
    if set(reference['cases'])!={case['id'] for case in data['cases']}:
        raise ArithmeticError('The reference omits or adds an input.')
    outputs = {}
    for case in data['cases']:
        record = reference['cases'][case['id']]
        if record['source']!=list(map(str,case['source'])) or record['target']!=list(map(str,case['target'])):
            raise ArithmeticError('Reference curves differ from the frozen input.')
        verify_constant(record)
        outputs[case['id']] = {'geometric_maps':record['geometric_maps'],
                              'rational_maps':sorted(record['rational_maps'],key=lambda value:json.dumps(value,sort_keys=True))}
    return outputs


def summary(rows):
    result = benchmark.summarize(rows)
    for method in METHODS:
        result[method]['wrong_answers'] = sum(row['method']==method and row['status']=='wrong_answer' for row in rows)
    return result


def main(output,seconds,cpu,wait_for):
    if wait_for is not None:
        wait_for_study(wait_for)
    benchmark.METHODS = METHODS
    panel = HERE/'map-panel.json'
    data = json.loads(panel.read_text())
    protocol_path = output.with_name(output.stem+'-protocol.json')
    if output.exists() or protocol_path.exists():
        raise ValueError('Choose a fresh output path to preserve the frozen study.')
    reference_path = output.with_name(output.stem+'-reference.json')
    outputs = prepare_reference(reference_path,panel)
    source_names = ['bounded-correspondences/'+name for name in ('composition.py','finite.py','maps.py')]
    source_names += ['solver-comparison/worker.py','solver-comparison/benchmark.py']
    source_hashes = {name:sha256(HERE.parent/name) for name in source_names}
    native_extractor = os.environ.get('CREST_SINGULAR') or shutil.which('Singular')
    protocol = {
        'scope':data['task'],'panel_sha256':sha256(panel),
        'worker_sha256':sha256(HERE/'map_worker.py'),'runner_sha256':sha256(Path(__file__)),
        'research_source_sha256':source_hashes,
        'reference_file':reference_path.name,'reference_sha256':sha256(reference_path),
        'methods':list(METHODS),'replicates':2,'per_process_wall_limit_seconds':seconds,
        'memory_limit_bytes':4*1024**3,'worker_cpu':cpu,
        'timing':'Serial fresh processes. Solver timer starts after imports and includes curve validation, coefficient construction, elimination, extraction of all rational maps and exact curve identities. Process wall includes interpreter startup. Frozen reference construction and verification are untimed.',
        'order':'Two seeded shuffles of every case/method pair, seed 76124.',
        'baseline':'Singular receives the full coefficient ideal G(U)=F*c^2 with inverse*u_d=1 and uses heuristic groebner in dp, radical and FGLM. Stock SymPy uses the fixed-degree dense routines underlying public decompose. The separately labelled corrected-dense-decomposition changes the right-component recurrence to i-(r+1)*j and retains the same dense division and residual equations. It is a project-side correction, not stock SymPy.',
        'correctness':'Every output is compared with the frozen, independently replayed exact reference. Incorrect outputs retain their actual values and runtimes with status wrong_answer. They do not count as completed correct answers and do not enter paired speed ratios.',
        'shared_extraction':{'algorithm':'The three decomposition methods use the same exact residual affine equations and rational-point engine.',
                             'engine':'Singular' if native_extractor else 'SymPy',
                             'native_executable':native_extractor},
        'platform':platform_record(),'created_utc':datetime.now(timezone.utc).isoformat(),
    }
    protocol_path.write_text(json.dumps(protocol,indent=2)+'\n')
    print('protocol',sha256(protocol_path),flush=True)
    rows = []
    rng = random.Random(76124)
    for replicate in (1,2):
        order = [(case['id'],method) for case in data['cases'] for method in METHODS]
        rng.shuffle(order)
        for case,method in order:
            row,result = run_one(case,method,replicate,panel,seconds,cpu,HERE/'map_worker.py')
            if result is not None:
                actual = {key:result[key] for key in ('geometric_maps','rational_maps')}
                if actual!=outputs[case]:
                    row['status'] = 'wrong_answer'
                    row['actual_output'] = actual
                    row['expected_output'] = outputs[case]
            rows.append(row)
            payload = {'protocol_sha256':sha256(protocol_path),'runs':rows,'outputs':outputs,
                       'summary':summary(rows),'completed':len(rows)==2*len(data['cases'])*len(METHODS)}
            output.write_text(json.dumps(payload,indent=2)+'\n')
            print(case,method,replicate,row['status'],round(row['process_wall_seconds'],3),flush=True)
    if source_hashes!={name:sha256(HERE.parent/name) for name in source_names}:
        raise ArithmeticError('Research source changed during timing. Preserve this run and start a fresh study.')
    if protocol['worker_sha256']!=sha256(HERE/'map_worker.py') or protocol['reference_sha256']!=sha256(reference_path):
        raise ArithmeticError('The timed adapter or exact reference changed.')
    payload['dependencies_unchanged'] = True
    payload['completed_utc'] = datetime.now(timezone.utc).isoformat()
    output.write_text(json.dumps(payload,indent=2)+'\n')
    print(json.dumps(summary(rows),indent=2),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'map-results.json')
    parser.add_argument('--seconds',type=float,default=10.0)
    parser.add_argument('--cpu',type=int)
    parser.add_argument('--wait-for',type=Path)
    args = parser.parse_args()
    main(args.output,args.seconds,args.cpu,args.wait_for)
