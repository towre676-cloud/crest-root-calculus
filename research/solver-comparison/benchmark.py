"""Run frozen marked-quartic comparisons serially with wall limits and exact outputs."""
import argparse
from collections import defaultdict
from datetime import datetime, timezone
import hashlib
import json
import os
from pathlib import Path
import platform
import random
import resource
import signal
import statistics
import subprocess
import sys
import tempfile
import time

HERE = Path(__file__).resolve().parent
METHODS = ('crest', 'singular', 'sympy-buchberger')


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def platform_record():
    import sympy
    cpu = 'unavailable'
    cpuinfo = Path('/proc/cpuinfo')
    if cpuinfo.exists():
        cpu = next((line.split(':',1)[1].strip() for line in cpuinfo.read_text().splitlines()
                    if line.startswith('model name')),cpu)
    singular = os.environ.get('SINGULAR','Singular')
    version = subprocess.run([singular,'--version'],text=True,capture_output=True,check=True).stdout.splitlines()[0]
    return {'python':sys.version.split()[0], 'sympy':sympy.__version__, 'singular':version,
            'os':platform.platform(), 'cpu':cpu, 'cpu_affinity':sorted(os.sched_getaffinity(0))
            if hasattr(os,'sched_getaffinity') else None}


def child_limit():
    # Four GiB excludes unbounded runaway allocations while keeping identical limits.
    resource.setrlimit(resource.RLIMIT_AS,(4*1024**3,4*1024**3))


def run_one(case,method,replicate,panel,seconds,cpu,worker=None):
    with tempfile.TemporaryDirectory(prefix='crest-timing-') as directory:
        output = Path(directory)/'result.json'
        command = [sys.executable,str(worker or HERE/'worker.py'),str(panel),case,method,str(output)]
        env = os.environ.copy()
        env.update({'PYTHONHASHSEED':'0','OMP_NUM_THREADS':'1','OPENBLAS_NUM_THREADS':'1'})
        start = time.perf_counter()
        process = subprocess.Popen(command,text=True,stdout=subprocess.PIPE,stderr=subprocess.PIPE,
                                   env=env,start_new_session=True,preexec_fn=child_limit)
        if cpu is not None:
            os.sched_setaffinity(process.pid,{cpu})
        try:
            stdout,stderr = process.communicate(timeout=seconds)
        except subprocess.TimeoutExpired:
            os.killpg(process.pid,signal.SIGKILL)
            stdout,stderr = process.communicate()
            return {'case':case,'method':method,'replicate':replicate,'status':'timeout',
                    'process_wall_seconds':time.perf_counter()-start,'limit_seconds':seconds},None
        common = {'case':case,'method':method,'replicate':replicate,
                  'process_wall_seconds':time.perf_counter()-start}
        if process.returncode or not output.exists():
            return {**common,'status':'error','returncode':process.returncode,
                    'error':(stdout+stderr)[-2000:]},None
        result = json.loads(output.read_text())
        values = result.pop('result')
        return {**common,'status':'complete',**result,
                'geometric_marks':values['geometric_marks'],
                'support_sha256':hashlib.sha256(json.dumps(values['basis'],sort_keys=True).encode()).hexdigest()},values


def summarize(rows):
    groups = defaultdict(list)
    for row in rows:
        groups[row['method'],row['case']].append(row)
    summary = {}
    for method in METHODS:
        full = [vals for (m,c),vals in groups.items() if m==method]
        completed = [vals for vals in full if len(vals)==2 and all(v['status']=='complete' for v in vals)]
        summary[method] = {'cases_completed_twice':len(completed), 'total_cases':len(full),
                           'timeouts':sum(v['status']=='timeout' for vals in full for v in vals),
                           'errors':sum(v['status']=='error' for vals in full for v in vals)}
        if completed:
            summary[method]['median_cpu_seconds'] = statistics.median(
                statistics.median(v['solve_cpu_seconds'] for v in vals) for vals in completed)
    return summary


def main(output,seconds,cpu):
    panel = HERE/'quartic-panel.json'
    data = json.loads(panel.read_text())
    protocol = {
        'scope':data['task'], 'panel_sha256':sha256(panel), 'worker_sha256':sha256(HERE/'worker.py'),
        'methods':list(METHODS), 'replicates':2, 'per_process_wall_limit_seconds':seconds,
        'memory_limit_bytes':4*1024**3, 'timing':'Serial fresh processes. Solver timer starts after imports and includes input normalization, smoothness check, exact elimination, radical support and map reconstruction. Process wall includes interpreter startup.',
        'order':'Two seeded shuffles of every case/method pair, seed 76123.',
        'baseline':'Same h-open marked chart. General solvers receive q, Db-B and four original square equations with invertible h factors removed. Singular heuristic groebner in dp, radical, then FGLM. SymPy Buchberger in grevlex, then FGLM and exact fibre radical.',
        'correctness':'Complete outputs must have identical rational lex ideals in c,s and identical geometric counts. The common original-equation fibre reconstruction verifies the target identities.',
        'platform':platform_record(), 'created_utc':datetime.now(timezone.utc).isoformat(),
    }
    protocol_path = output.with_name(output.stem+'-protocol.json')
    if output.exists() or protocol_path.exists():
        raise ValueError('Choose a fresh output path to preserve the frozen study.')
    protocol_path.write_text(json.dumps(protocol,indent=2)+'\n')
    print('protocol',sha256(protocol_path),flush=True)
    rows, supports = [], {}
    rng = random.Random(76123)
    for replicate in (1,2):
        order = [(case['id'],method) for case in data['cases'] for method in METHODS]
        rng.shuffle(order)
        for case,method in order:
            row,result = run_one(case,method,replicate,panel,seconds,cpu)
            rows.append(row)
            if result is not None:
                previous = supports.get(case)
                if previous and previous['basis']!=result['basis']:
                    raise ArithmeticError(f'Geometric support disagreement on {case}.')
                if previous and previous['geometric_marks']!=result['geometric_marks']:
                    raise ArithmeticError(f'Geometric count disagreement on {case}.')
                supports[case] = {key:result[key] for key in ('basis','geometric_marks')}
            payload = {'protocol_sha256':sha256(protocol_path), 'runs':rows,
                       'supports':supports, 'summary':summarize(rows),
                       'completed':len(rows)==2*len(data['cases'])*len(METHODS)}
            output.write_text(json.dumps(payload,indent=2)+'\n')
            print(case,method,replicate,row['status'],round(row['process_wall_seconds'],3),flush=True)
    print(json.dumps(summarize(rows),indent=2),flush=True)


if __name__=='__main__':
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument('--output',type=Path,default=HERE/'quartic-results.json')
    parser.add_argument('--seconds',type=float,default=30.0)
    parser.add_argument('--cpu',type=int)
    args = parser.parse_args()
    main(args.output,args.seconds,args.cpu)
