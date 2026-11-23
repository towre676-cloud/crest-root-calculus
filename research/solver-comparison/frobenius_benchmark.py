"""Compare complete local zeta numerators with native PARI/GP."""
import argparse
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import resource
import signal
import subprocess
import sys
import tempfile
import time

import sympy as sp

HERE = Path(__file__).resolve().parent
FIELD_BUDGET = 2_000_000
sys.path.insert(0,str(HERE.parent/'bounded-correspondences'))
from frobenius_filter import curve, local_data, reduction
from worker import cpu_time
from map_benchmark import wait_for_study
from benchmark import child_limit, sha256


def create_panel(path):
    if path.exists():
        raise ValueError('A frozen Frobenius panel already exists.')
    rng = random.Random(2026112323)
    cases = []
    choices = [(g,2*g+1,p) for g in (2,3,4) for p in (5,7,11)]
    choices += [(g,2*g+2,p) for g,p in ((2,5),(3,7),(4,11))]
    for genus,degree,prime in choices:
        while True:
            coefficients = [rng.randint(-4,4) for _ in range(degree)]+[rng.choice((1,2,3))]
            try:
                f = curve(coefficients)
            except ValueError:
                continue
            if reduction(f,prime) is not None:
                break
        cases.append({'id':f'g{genus}-d{degree}-p{prime}','genus':genus,
                      'source':coefficients,'prime':prime})
    small_ids = [case['id'] for case in cases]
    large_rng = random.Random(2026112324)
    large_ids = []
    for genus,prime in ((2,1009),(3,97),(4,31)):
        degree = 2*genus+1
        while True:
            coefficients = [large_rng.randint(-4,4) for _ in range(degree)]+[large_rng.choice((1,2,3))]
            try:
                f = curve(coefficients)
            except ValueError:
                continue
            if reduction(f,prime) is not None:
                break
        identifier = f'g{genus}-d{degree}-p{prime}'
        large_ids.append(identifier)
        cases.append({'id':identifier,'genus':genus,'source':coefficients,'prime':prime})
    data = {'schema':2,'seed':2026112323,'large_field_seed':2026112324,
            'scope':'complete local zeta numerator of a smooth hyperelliptic curve',
            'selection_rule':'Preserve the twelve original seeded small-prime inputs. Append three independently seeded smooth odd-degree inputs of genera2/3/4 at primes1009/97/31, respectively, before any new timing; reject bad reduction only. The two cohorts are declared before measurement.',
            'field_element_budget':FIELD_BUDGET,
            'cohorts':{'small-prime':small_ids,'large-field':large_ids},
            'cases':cases}
    path.write_text(json.dumps(data,indent=2)+'\n')


def solve(case,method):
    f = curve(case['source'])
    p = case['prime']
    if reduction(f,p) is None:
        raise ValueError('The panel contains bad reduction.')
    if method=='crest-direct':
        return local_data(f,p,FIELD_BUDGET)
    gp = os.environ.get('PARI_GP','gp')
    expression = '+'.join(f'({coefficient})*x^{i}' for i,coefficient in enumerate(case['source']))
    script = f'print(Vec(hyperellcharpoly(Mod(1,{p})*({expression}))));\nquit\n'
    completed = subprocess.run([gp,'-q','-f','-s','64000000'],input=script,text=True,
                               capture_output=True,check=True)
    lines = [line.strip() for line in completed.stdout.splitlines() if line.startswith('[')]
    if len(lines)!=1 or '***' in completed.stderr:
        raise ArithmeticError(completed.stdout+completed.stderr)
    numerator = json.loads(lines[0])
    g = case['genus']
    if len(numerator)!=2*g+1 or numerator[0]!=1 or numerator[-1]!=p**g:
        raise ArithmeticError('Unexpected Frobenius polynomial normalization.')
    traces = []
    for k in range(1,g+1):
        traces.append(-k*numerator[k]-sum(numerator[j]*traces[k-j-1] for j in range(1,k)))
    return {'counts':[p**k+1-traces[k-1] for k in range(1,g+1)],
            'zeta_numerator':numerator}


def worker(panel,case,method,output):
    entry = next(row for row in json.loads(panel.read_text())['cases'] if row['id']==case)
    wall,cpu = time.perf_counter(),cpu_time()
    result = solve(entry,method)
    output.write_text(json.dumps({'result':result,'solve_wall_seconds':time.perf_counter()-wall,
                                 'solve_cpu_seconds':cpu_time()-cpu},indent=2)+'\n')


def study(output,cpu,wait_for):
    if wait_for is not None:
        wait_for_study(wait_for)
    panel = HERE/'frobenius-panel.json'
    data = json.loads(panel.read_text())
    names = ['bounded-correspondences/frobenius_filter.py','cubic-tower-factors/field.py',
             'cubic-tower-factors/tower.py','bounded-correspondences/finite.py',
             'solver-comparison/worker.py','solver-comparison/benchmark.py']
    hashes = {name:sha256(HERE.parent/name) for name in names}
    protocol = {'scope':data['scope'],'panel_sha256':sha256(panel),
                'worker_sha256':sha256(Path(__file__)), 'research_source_sha256':hashes,
                'methods':['crest-direct','pari'],'replicates':2,'worker_cpu':cpu,
                'field_element_budget':FIELD_BUDGET,'cohorts':data['cohorts'],
                'cohort_definition':'The12originalsmall-primecasesand3larger-fieldcasesarefixedbeforetiming. Maximumextensionordersinthelatterare1009^2,97^3and31^4.',
                'per_process_wall_limit_seconds':10,'memory_limit_bytes':4*1024**3,
                'timing':'Serial fresh processes. Includes curve validation, all field construction and complete zeta numerator. PARI invocation and decoding are included. Imports excluded from solver timer and included in process wall.',
                'cpu_timer':'POSIX getrusage user plus system time for the Python worker and reaped children; returned in seconds with microsecond precision. The GP millisecond timer is not used.',
                'normalization':'The descending coefficients of PARI characteristic polynomial are the ascending coefficients of the zeta numerator. Extension counts recovered by Newton identities.',
                'pari':subprocess.check_output([os.environ.get('PARI_GP','gp'),'--version-short'],text=True).strip(),
                'python':sys.version.split()[0],'sympy':sp.__version__,'numpy':__import__('numpy').__version__,'platform':platform.platform(),
                'created_utc':datetime.now(timezone.utc).isoformat()}
    protocol_path = output.with_name(output.stem+'-protocol.json')
    if output.exists() or protocol_path.exists():
        raise ValueError('Choose a fresh output path.')
    protocol_path.write_text(json.dumps(protocol,indent=2)+'\n')
    print('protocol',sha256(protocol_path),flush=True)
    rows,outputs = [],{}
    rng = random.Random(76125)
    for replicate in (1,2):
        order = [(row['id'],method) for row in data['cases'] for method in protocol['methods']]
        rng.shuffle(order)
        for case,method in order:
            with tempfile.TemporaryDirectory(prefix='crest-frobenius-timing-') as directory:
                result_path = Path(directory)/'result.json'
                env = os.environ.copy();env.update({'OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1','PYTHONHASHSEED':'0'})
                command = [sys.executable,str(Path(__file__)),'--worker',str(panel),case,method,str(result_path)]
                start = time.perf_counter()
                process = subprocess.Popen(command,stdout=subprocess.PIPE,stderr=subprocess.PIPE,text=True,
                                           env=env,start_new_session=True,preexec_fn=child_limit)
                if cpu is not None:
                    os.sched_setaffinity(process.pid,{cpu})
                row = {'case':case,'method':method,'replicate':replicate}
                try:
                    stdout,stderr = process.communicate(timeout=10)
                    row['process_wall_seconds'] = time.perf_counter()-start
                    if process.returncode or not result_path.exists():
                        row.update({'status':'error','error':(stdout+stderr)[-1500:]})
                    else:
                        payload = json.loads(result_path.read_text());result = payload.pop('result')
                        row.update({'status':'complete',**payload})
                        row['result_sha256'] = hashlib.sha256(json.dumps(result,sort_keys=True).encode()).hexdigest()
                        if case in outputs and outputs[case]!=result:
                            raise ArithmeticError(f'Frobenius polynomial disagreement on {case}.')
                        outputs[case] = result
                except subprocess.TimeoutExpired:
                    os.killpg(process.pid,signal.SIGKILL);process.communicate()
                    row.update({'status':'timeout','process_wall_seconds':time.perf_counter()-start,'limit_seconds':10})
            rows.append(row)
            output.write_text(json.dumps({'protocol_sha256':sha256(protocol_path),'runs':rows,'outputs':outputs,
                                        'completed':len(rows)==4*len(data['cases'])},indent=2)+'\n')
            print(case,method,replicate,row['status'],round(row['process_wall_seconds'],3),flush=True)
    if hashes!={name:sha256(HERE.parent/name) for name in names}:
        raise ArithmeticError('Research source changed during timing.')
    payload = json.loads(output.read_text())
    payload['dependencies_unchanged'] = True
    payload['completed_utc'] = datetime.now(timezone.utc).isoformat()
    output.write_text(json.dumps(payload,indent=2)+'\n')


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--create-panel',action='store_true')
    parser.add_argument('--worker',nargs=4,metavar=('PANEL','CASE','METHOD','OUTPUT'))
    parser.add_argument('--output',type=Path,default=HERE/'frobenius-results.json')
    parser.add_argument('--cpu',type=int)
    parser.add_argument('--wait-for',type=Path)
    args = parser.parse_args()
    if args.create_panel:
        create_panel(HERE/'frobenius-panel.json')
    elif args.worker:
        panel,case,method,output = args.worker
        worker(Path(panel),case,method,Path(output))
    else:
        study(args.output,args.cpu,args.wait_for)
